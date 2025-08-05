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
#include <limits>
#include <thread>

constexpr int8_t OFFESTS_3X3[27][3] = {
    {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1},
    {-1, 1, -1},  {-1, 1, 0},  {-1, 1, 1},  {0, -1, -1}, {0, -1, 0}, {0, -1, 1},
    {0, 0, -1},   {0, 0, 0},   {0, 0, 1},   {0, 1, -1},  {0, 1, 0},  {0, 1, 1},
    {1, -1, -1},  {1, -1, 0},  {1, -1, 1},  {1, 0, -1},  {1, 0, 0},  {1, 0, 1},
    {1, 1, -1},   {1, 1, 0},   {1, 1, 1}};

static double W(double r, double h) {
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

static double gradW(double r, double h) {
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

static inline vec3<double> a(const Particle *p) { return p->force / p->mass; }

static void clamp_velosity(Particle *p, double vmax) {
    p->velocity.x = std::clamp(p->velocity.x, -vmax, vmax);
    p->velocity.y = std::clamp(p->velocity.y, -vmax, vmax);
    p->velocity.z = std::clamp(p->velocity.z, -vmax, vmax);
}

static void integrate_semi_implicit_euler(Particle *p, double dt, double vmax) {
    p->velocity += a(p) * dt;

    // CFL limiter
    clamp_velosity(p, vmax);

    p->position += p->velocity * dt;
}

static void integrateRK4(Particle *p, double dt, double vmax) {
    auto k1_v = a(p);
    auto k1_x = p->velocity;

    auto k2_v = a(p); // force assumed constant over dt
    auto k2_x = p->velocity + 0.5 * dt * k1_v;

    auto k3_v = a(p);
    auto k3_x = p->velocity + 0.5 * dt * k2_v;

    auto k4_v = a(p);
    auto k4_x = p->velocity + dt * k3_v;

    p->velocity += (dt / 6.0) * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v);

    // CFL limiter
    clamp_velosity(p, vmax);

    p->position += (dt / 6.0) * (k1_x + 2.0 * k2_x + 2.0 * k3_x + k4_x);
}

/**
 * @brief Calls the inner function for each pair of particles within the same
 * and neighboring cells. Each pair is guaranteed to be called exactly once.
 */
static void
for_pair_neighbor_cells(const Grid &grid,
                        std::function<void(Particle *, Particle *)> inner) {
    for (size_t i = 0; i < grid.grid_dims_.z; ++i) {
        for (size_t j = 0; j < grid.grid_dims_.y; ++j) {
            for (size_t k = 0; k < grid.grid_dims_.x; ++k) {
                auto &particles_cell_1 = grid.particles_in_cell({i, j, k});
                if (particles_cell_1.empty())
                    continue;
                for (auto [di, dj, dk] : OFFESTS_3X3) {
                    if (((i == 0) && (di < 0)) || ((j == 0) && (dj < 0)) ||
                        ((k == 0) && (dk < 0)))
                        continue;
                    size_t i1 = i + di, j1 = j + dj, k1 = k + dk;
                    if ((i1 >= grid.grid_dims_.z) ||
                        (j1 >= grid.grid_dims_.y) || (k1 >= grid.grid_dims_.x))
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

struct SimulationParameters {
    double viscosity = 0.01;
    double gravity = 10;
    double specific_volume = 2e6;
};
const double ARTIFICIAL_VISC = 0;
const double RHO_0 = -1.;

static double compute_pressure(const Particle *p, double k) {
    double density_diff = std::max(p->density - p->initial_density, 0.) / p->initial_density;
    // double density_diff = (p->density - p->initial_density) / p->initial_density;

    return k * density_diff;
}

class Simulation {

    ExternalBoundaries boundaries_;
    SimulationParameters parameters_;

    void update_density_pair_cells(Particle *p1, Particle *p2) const {
        vec3<double> r = p2->position - p1->position;
        double r_norm = std::sqrt(r.dot(r));
        double W_ = W(r_norm, 4 * PARTICLE_RADIUS);
        p1->density += p2->mass * W_;
        p2->density += p1->mass * W_;
    }

    void apply_artificial_viscosity_pair_cells(Particle *p1,
                                               Particle *p2) const {
        // artificial viscosity to reduce oscillations
        auto r = p2->position - p1->position;
        double r_norm = std::sqrt(r.dot(r));
        double W_ = W(r_norm, 4 * PARTICLE_RADIUS);

        auto delta_v = p2->velocity - p1->velocity;
        auto rho = (p1->density + p2->density) / 2;

        p1->velocity += ARTIFICIAL_VISC * p2->mass / rho * delta_v * W_;
        p2->velocity -= ARTIFICIAL_VISC * p1->mass / rho * delta_v * W_;
    }

    void update_densities() {
        // First compute density influence of each particle to itself.
        for (size_t i = 0; i < particles_.size(); ++i) {
            Particle *p = particles_[i].get();
            double W_ = W(0, PARTICLE_RADIUS * 2);
            p->density = p->mass * W_;
        }
        // Compute influence of nearby particles.
        for_pair_neighbor_cells(grid_, [&](Particle *p1, Particle *p2) {
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
        double gradW_ = gradW(r_norm, 4 * PARTICLE_RADIUS);
        double pressure_1 = compute_pressure(p1, parameters_.specific_volume);
        double pressure_2 = compute_pressure(p2, parameters_.specific_volume);
        double tmp = pressure_1 / (p1->density * p1->density) +
                     pressure_2 / (p1->density * p2->density);
        tmp *= gradW_;
        p1->force += p2->mass * tmp * r;
        p2->force -= p1->mass * tmp * r;
        assert(!isnan(p1->force));
        assert(!isnan(p2->force));
    }

    void update_gravity_force(Particle *p1) const {
        p1->force.y -= parameters_.gravity * p1->mass;
    }

    void update_viscous_force(Particle *p1, Particle *p2) const {
        vec3<double> r = p2->position - p1->position;
        double r_norm = std::sqrt(r.dot(r));
        const double tol = 1e-9;
        if (r_norm < tol) {
            return;
        }
        double gradW_ = gradW(r_norm, 4 * PARTICLE_RADIUS);
        const double visc = parameters_.viscosity;
        // F_i = m_i * visc * Sum_j (m_j / rho_j * (v_j - v_i) * 2 ||nabla W||
        // / ||r||)
        auto visc_force = (-1) * p1->mass * visc * p2->mass *
                          (p2->velocity - p1->velocity) * 2 * gradW_;
        p1->force += visc_force / p2->density;
        p2->force -= visc_force / p1->density;
    }

    void update_positions(double dt) {
        double min_cell_size =
            std::min(grid_.grid_cell_size_.x, grid_.grid_cell_size_.y);
        min_cell_size = std::min(min_cell_size, grid_.grid_cell_size_.z);
        double vmax = 0.4 * min_cell_size / dt;

        for (auto &particle : particles_) {
            assert(!isnan(particle->force));

            auto old_grid_cell = grid_.position_to_cell(particle->predicted_position);
            const auto old_position = particle->position;

            // Update velocity using forces
            integrate_semi_implicit_euler(particle.get(), dt, vmax);

            // Reset force for the next iteration
            particle->force = {0, 0, 0};

            // Handling 2D simulation
            particle->position.z = grid_.domain_limits_.z / 2;
            particle->velocity.z = 0;

            // Collisions
            boundaries_.handle_collision(old_position, particle->position,
                                         particle->velocity);

            assert(grid_.within_domain_bounds(particle->position));

            // Update particle position in the grid
            update_particle_position_in_grid(particle.get(), old_grid_cell);
        }
    }

    void update_particle_position_in_grid(Particle *p,
                                          const vec3<size_t> &old_grid_cell) {
        auto new_grid_cell = grid_.position_to_cell(p->position);
        if (new_grid_cell != old_grid_cell) {
            grid_.remove_particle(p, old_grid_cell);
            grid_.add_particle(p, p->position);
        }
    }

  public:
    std::vector<std::unique_ptr<Particle>> particles_;
    Grid grid_;

    Simulation() = delete;
    Simulation(Grid grid, std::vector<std::unique_ptr<Particle>> particles)
        : boundaries_(grid.domain_limits_), particles_(std::move(particles)),
          grid_(std::move(grid)) {
        for (auto &particle : particles_) {
            grid_.add_particle(particle.get(), particle->position);
        }
    }

    void add_particle(const vec3<double> position, double density) {
        particles_.push_back(
            std::make_unique<Particle>(make_particle(position, density)));
        particles_.back()->velocity = {
            0,
            0,
            static_cast<double>(rand()) / RAND_MAX * 1e-1,
        };

        grid_.add_particle(particles_.back().get(), particles_.back()->position);
    }

    void remove_particle_at(const vec3<double> position) {
        const auto cell = grid_.position_to_cell(position);
        const auto &particles_in_cell = grid_.particles_in_cell(cell);

        // search for closest particle.
        Particle *particle_to_remove = nullptr;
        double minimal_dist = std::numeric_limits<double>::infinity();
        for (const auto p : particles_in_cell) {
            auto dst = p->position - position;
            auto dst_norm = dst.dot(dst);
            if (dst_norm < minimal_dist) {
                minimal_dist = dst_norm;
                particle_to_remove = p;
            }
        }
        if (particle_to_remove == nullptr) {
            return;
        }

        // removing from grid.
        grid_.remove_particle(particle_to_remove, cell);

        // removing from particles.
        auto it = std::find_if(particles_.begin(), particles_.end(),
                               [&particle_to_remove](const auto &p1) {
                                   return p1.get() == particle_to_remove;
                               });
        if (it != particles_.end()) {
            particles_.erase(it);
        }
    }

    void apply_central_force(vec3<double> center, double acceleration,
                             double radius) {
        for (auto &particle : particles_) {
            auto r = center - particle->position;
            auto distance = std::sqrt(r.dot(r));
            if (distance >= radius) {
                continue;
            }
            particle->force += particle->mass * acceleration * r;
        }
    }

    void set_viscosity(double viscosity) { parameters_.viscosity = viscosity; }

    void set_gravity(double gravity) { parameters_.gravity = gravity; }

    void set_specific_volume(double specific_volume) {
        parameters_.specific_volume = specific_volume;
    }

    void update(double dt) {
        for (auto &p : particles_) {
            auto old_grid_cell = grid_.position_to_cell(p->position);
            update_gravity_force(p.get());
            p->velocity += dt * a(p.get());
            double artificial_dt = 1 / 60.;
            
            // p->predicted_position = p->position;
            p->predicted_position = p->position + p->velocity * artificial_dt;
            
            p->force = {0, 0, 0};
            
            auto predicted_grid_cell = grid_.position_to_cell(p->predicted_position);
            if (predicted_grid_cell != old_grid_cell) {
                grid_.remove_particle(p.get(), old_grid_cell);
                grid_.add_particle(p.get(), p->predicted_position);
            }
        }

        update_densities();
        for_pair_neighbor_cells(grid_, [&](Particle *p1, Particle *p2) {
            update_pressure_force(p1, p2);
            update_viscous_force(p1, p2);
        });
        // for_pair_neighbor_cells(grid_, [&](Particle *p1, Particle *p2) {
            // apply_artificial_viscosity_pair_cells(p1, p2);
        // });
        update_positions(dt);
    }

    std::vector<Particle *> get_particles_raw() {
        std::vector<Particle *> ptrs;
        ptrs.reserve(particles_.size());
        for (auto &uptr : particles_) {
            ptrs.push_back(uptr.get());
        }
        return ptrs;
    }

    vec3<double> get_domain_limits() { return grid_.domain_limits_; }
};
