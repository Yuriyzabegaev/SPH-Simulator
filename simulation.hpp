#pragma once
#include "boundary.hpp"
#include "grid.hpp"
#include "print.hpp"
#include "vec3.hpp"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <functional>
#include <thread>

double distL2(const vec3<double> &a, const vec3<double> &b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

constexpr int8_t OFFESTS_3X3[27][3] = {
    {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1},
    {-1, 1, -1},  {-1, 1, 0},  {-1, 1, 1},  {0, -1, -1}, {0, -1, 0}, {0, -1, 1},
    {0, 0, -1},   {0, 0, 0},   {0, 0, 1},   {0, 1, -1},  {0, 1, 0},  {0, 1, 1},
    {1, -1, -1},  {1, -1, 0},  {1, -1, 1},  {1, 0, -1},  {1, 0, 0},  {1, 0, 1},
    {1, 1, -1},   {1, 1, 0},   {1, 1, 1}};

double W(double r, double h) {
    double q = r / h;
    double sigma_3d = 8 / (M_PI * h * h * h);
    assert(q >= 0);
    if (q <= 0.5) {
        return sigma_3d * (6 * (q * q * q - q * q) + 1);
    } else if (q <= 1) {
        double tmp = 1 - q;
        return sigma_3d * (2 * tmp * tmp * tmp);
    } else {
        return 0;
    }
}

double gradW(double r, double h) {
    double q = r / h;
    double sigma_3d = 8 / (M_PI * h * h * h);
    assert(q >= 0);
    if (q <= 0.5) {
        return sigma_3d * 6 * q * (3 * q - 2) / (h * r);
    } else if (q <= 1) {
        double tmp = 1 - q;
        return -sigma_3d * 6 * tmp * tmp / (h * r);
    } else {
        return 0;
    }
}

double compute_pressure(const Particle *p) {
    const double RHO_0 = 1000;
    const double k = 100;
    double density_diff = std::max(p->density - RHO_0, 0.);
    return k * density_diff;
}

void for_pair_neighbor_cells(
    const Grid &grid, std::function<void(Particle *, Particle *)> inner) {
    for (size_t i = 0; i < grid.m_grid_dims.z; ++i) {
        for (size_t j = 0; j < grid.m_grid_dims.y; ++j) {
            for (size_t k = 0; k < grid.m_grid_dims.x; ++k) {
                auto &particles_cell_1 = grid.particles_in_cell({i, j, k});
                if (particles_cell_1.empty())
                    continue;
                for (auto [di, dj, dk] : OFFESTS_3X3) {
                    if (((i == 0) && (di < 0)) || ((j == 0) && (dj < 0)) ||
                        ((k == 0) && (dk < 0)))
                        continue;
                    size_t i1 = i + di, j1 = j + dj, k1 = k + dk;
                    if ((i1 >= grid.m_grid_dims.z) ||
                        (j1 >= grid.m_grid_dims.y) ||
                        (k1 >= grid.m_grid_dims.x))
                        continue;

                    auto &particles_cell_2 =
                        grid.particles_in_cell({i1, j1, k1});
                    for (auto &p2 : particles_cell_2) {
                        for (auto &p1 : particles_cell_1) {
                            if (p1 < p2) {
                                inner(p1, p2);
                            }
                        }
                    }
                }
            }
        }
    }
}

class Simulation {

    ExternalBoundaries m_boundaries;

    void update_density_pair_cells(Particle *p1, Particle *p2) const {
        vec3<double> r = p2->position - p1->position;
        double r_norm = std::sqrt(r.dot(r));
        double collision_dist = p1->radius + p2->radius;
        double W_ = W(r_norm, collision_dist);
        p1->density += p2->mass * W_;
        p2->density += p1->mass * W_;
    }

    void update_densities() {
        for (size_t i = 0; i < m_particles.size(); ++i) {
            Particle *p = m_particles[i].get();
            double W_ = W(0, p->radius * 2);
            p->density = p->mass * W_;
        }
        for_pair_neighbor_cells(m_grid, [&](Particle *p1, Particle *p2) {
            update_density_pair_cells(p1, p2);
        });
    }

    void update_pressure_force(Particle *p1, Particle *p2) const {
        vec3<double> r = p2->position - p1->position;
        double r_norm = std::sqrt(r.dot(r));
        const double tol = 1e-9;
        if (r_norm < tol) {
            return;
        }
        double collision_dist = p1->radius + p2->radius;
        double gradW_ = gradW(r_norm, collision_dist);
        double pressure_1 = compute_pressure(p1);
        double pressure_2 = compute_pressure(p2);
        double tmp = pressure_1 / (p1->density * p1->density) +
                     pressure_2 / (p1->density * p2->density);
        tmp *= gradW_;
        p1->force += p2->mass * tmp * r;
        p2->force -= p1->mass * tmp * r;
        assert(!isnan(p1->force));
        assert(!isnan(p2->force));
    }

    void update_gravity_force(Particle *p1) const {
        p1->force.y -= 9.8 * p1->mass;
    }

    void update_viscous_force(Particle *p1, Particle *p2) const {
        vec3<double> r = p2->position - p1->position;
        double r_norm = std::sqrt(r.dot(r));
        const double tol = 1e-9;
        if (r_norm < tol) {
            return;
        }
        double collision_dist = p1->radius + p2->radius;
        double gradW_ = gradW(r_norm, collision_dist);
        const double visc = 0.01;
        // F_i = m_i * visc * Sum_j (m_j / rho_j * (v_j - v_i) * 2 ||nabla W||
        // / ||r||)
        auto visc_force = (-1) * p1->mass * visc * p2->mass *
                          (p2->velocity - p1->velocity) * 2 * gradW_;
        p1->force += visc_force / p2->density;
        p2->force -= visc_force / p1->density;
    }

    void update_forces() {
        for (size_t i = 0; i < m_particles.size(); ++i) {
            update_gravity_force(m_particles[i].get());
        }

        for_pair_neighbor_cells(m_grid, [&](Particle *p1, Particle *p2) {
            update_pressure_force(p1, p2);
            update_viscous_force(p1, p2);
        });
    }

    void update_positions(double dt) {
        for (auto &particle : m_particles) {
            assert(!isnan(particle->force));

            auto old_grid_cell = m_grid.position_to_cell(particle->position);
            const auto old_position = particle->position;

            particle->velocity += particle->force / particle->mass * dt;

            // // Velocity limiter
            // double max_velocity =
            //     std::min(m_grid.m_domain_limits.x, m_grid.m_domain_limits.y)
            //     / (4.0 * dt);
            // particle.velocity.z =
            //     std::clamp(particle.velocity.z, -max_velocity, max_velocity);
            // particle.velocity.y =
            //     std::clamp(particle.velocity.y, -max_velocity, max_velocity);
            // particle.velocity.x =
            //     std::clamp(particle.velocity.x, -max_velocity, max_velocity);

            particle->position += particle->velocity * dt;

            // Reset force for the next iteration
            particle->force = {0, 0, 0};

            // Collisions
            m_boundaries.handle_collision(old_position, particle->position,
                                          particle->velocity);

            particle->position.z = m_grid.m_domain_limits.z / 2;
            particle->velocity.z = 0;

            assert(m_grid.within_domain_bounds(particle->position));

            // Update particle position in the grid
            auto new_grid_cell = m_grid.position_to_cell(particle->position);
            if (new_grid_cell != old_grid_cell) {
                m_grid.remove_particle(particle.get(), old_grid_cell);
                m_grid.add_particle(particle.get());
            }
        }
    }

  public:
    std::vector<std::unique_ptr<Particle>> m_particles;
    Grid m_grid;

    Simulation() = delete;
    Simulation(Grid grid, std::vector<std::unique_ptr<Particle>> particles)
        : m_boundaries(grid.m_domain_limits), m_particles(std::move(particles)),
          m_grid(std::move(grid)) {
        for (auto &particle : m_particles) {
            m_grid.add_particle(particle.get());
        }
    }

    void add_particle(const Particle particle) {
        m_particles.emplace_back(
            std::make_unique<Particle>(std::move(particle)));
        m_grid.add_particle(m_particles.back().get());
    }

    void apply_external_force(vec3<double> position, vec3<double> acceleration,
                              double radius) {
        for (auto &particle : m_particles) {
            auto r = particle->position - position;
            auto distance = std::sqrt(r.dot(r));
            if (distance >= radius) {
                continue;
            }
            particle->force += particle->mass * acceleration;
        }
    }

    void update(double dt) {
        update_densities();
        update_forces();
        update_positions(dt);
    }
};
