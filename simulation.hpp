#pragma once
#include "boundary.hpp"
#include "grid.hpp"
#include "materials.hpp"
#include "print.hpp"
#include <algorithm>
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

void for_particle_neighbor_cells(const Grid &grid, Particle &particle,
                                 std::function<void(Particle &)> inner) {
    // Iterate over the neighboring cells of the particle
    // and apply the inner function to each particle in those cells, including
    // the same particle.
    auto grid_pos = grid.position_to_cell(particle.position);
    size_t xmin = grid_pos.x > 0 ? grid_pos.x - 1 : grid_pos.x;
    size_t ymin = grid_pos.y > 0 ? grid_pos.y - 1 : grid_pos.y;
    size_t zmin = grid_pos.z > 0 ? grid_pos.z - 1 : grid_pos.z;
    size_t xmax = std::min(grid_pos.x + 2, grid.m_grid_dims.x);
    size_t ymax = std::min(grid_pos.y + 2, grid.m_grid_dims.y);
    size_t zmax = std::min(grid_pos.z + 2, grid.m_grid_dims.z);
    for (size_t i = zmin; i < zmax; ++i) {
        for (size_t j = ymin; j < ymax; ++j) {
            for (size_t k = xmin; k < xmax; ++k) {
                for (auto &neighbor : grid.particles_in_cell({i, j, k})) {
                    inner(*neighbor);
                }
            }
        }
    }
}

void for_particle_pairs_neighbor_cells(
    const Grid &grid, std::function<void(Particle *, Particle *)> inner) {
    for (size_t pos_i = 0; pos_i < grid.m_grid_dims.z; ++pos_i) {
        for (size_t pos_j = 0; pos_j < grid.m_grid_dims.y; ++pos_j) {
            for (size_t pos_k = 0; pos_k < grid.m_grid_dims.x; ++pos_k) {
                auto grid_pos = grid.idx_3d(grid.idx_1d(pos_i, pos_j, pos_k));
                auto &particles_this_cell = grid.grid()[grid.idx_1d(grid_pos)];
                if (particles_this_cell.empty()) {
                    // No particles in this cell, nothing to do.
                    continue;
                }

                // Process the current grid cell
                for (auto it1 = particles_this_cell.begin();
                     it1 != particles_this_cell.end(); ++it1) {
                    auto it2 = it1;
                    ++it2;
                    for (; it2 != particles_this_cell.end(); ++it2) {
                        inner(*it1, *it2);
                    }
                }
                // Update particles in neighboring cells.
                size_t xmin = grid_pos.x > 0 ? grid_pos.x - 1 : grid_pos.x;
                size_t ymin = grid_pos.y > 0 ? grid_pos.y - 1 : grid_pos.y;
                size_t zmin = grid_pos.z > 0 ? grid_pos.z - 1 : grid_pos.z;
                for (size_t i = zmin; i <= grid_pos.z; ++i) {
                    for (size_t j = ymin; j <= grid_pos.y; ++j) {
                        for (size_t k = xmin; k <= grid_pos.x; ++k) {
                            if (i == grid_pos.z && j == grid_pos.y &&
                                k == grid_pos.x) {
                                // Skip the same cell, already processed.
                                continue;
                            }
                            for (auto *p2 : grid.particles_in_cell({i, j, k})) {
                                for (auto *p1 : particles_this_cell) {
                                    inner(p1, p2);
                                }
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

  public:
    std::vector<Particle> m_particles;
    Grid m_grid;

    Simulation() = delete;
    Simulation(Grid grid, std::vector<Particle> particles)
        : m_boundaries(grid.m_domain_limits), m_particles(std::move(particles)),
          m_grid(std::move(grid)) {
        for (auto &particle : m_particles) {
            m_grid.add_particle(&particle);
        }
    }

    // void add_particle(const Particle particle) {
    //     m_particles.push_back(std::move(particle));
    //     m_grid.add_particle(&m_particles.back());
    // }

    void update_pressure_force(Particle *p1, Particle *p2) const {
        vec3<double> force = {0, 0, 0};
        double tol = 1e-6;
        force.z = (p2->position.z - p1->position.z);
        force.y = (p2->position.y - p1->position.y);
        force.x = (p2->position.x - p1->position.x);

        double dist = distL2(p1->position, p2->position);
        double y0 = 3;
        double y1 = 0;
        double x0 = 0;
        double x1 = 1.25 * CELL_SIZE;
        if (dist < tol) {
            // Randomly adjust velocities
            force = {0.01 * (rand() % 100) / 100.0,
                     0.01 * (rand() % 100) / 100.0,
                     0.01 * (rand() % 100) / 100.0};
            force.z = 0;
            p1->velocity += force;
            p2->velocity -= force;
        } else if (dist < x1) {
            double v = y0 + (y1 - y0) / (x1 - x0) * (dist - x0);
            p1->velocity -= v * force;
            p2->velocity += v * force;
        }
    }

    void update_gravity_force(Particle *p1, Particle *p2) const {
        p1->force.y -= 0.00001;
        p2->force.y -= 0.00001;
    }

    void update_viscous_force(Particle *p1, Particle *p2) const {
        p1->velocity *= 0.99;
        p2->velocity *= 0.99;
    }

    void update_particles_pair_forces(Particle *p1, Particle *p2) const {
        update_pressure_force(p1, p2);
        update_gravity_force(p1, p2);
        update_viscous_force(p1, p2);
    }

    void update_forces() const {
        for_particle_pairs_neighbor_cells(
            m_grid, [&](Particle *p1, Particle *p2) {
                update_particles_pair_forces(p1, p2);
            });
    }

    void update_positions(double dt) {
        for (auto &particle : m_particles) {
            auto old_grid_cell = m_grid.position_to_cell(particle.position);
            const auto old_position = particle.position;

            particle.velocity += particle.force / particle.mass * dt;
            // Velocity limiter
            double max_velocity = std::min(std::min(m_grid.m_domain_limits.x,
                                                    m_grid.m_domain_limits.y),
                                           m_grid.m_domain_limits.z) /
                                  (4.0 * dt);
            particle.velocity.z =
                std::clamp(particle.velocity.z, -max_velocity, max_velocity);
            particle.velocity.y =
                std::clamp(particle.velocity.y, -max_velocity, max_velocity);
            particle.velocity.x =
                std::clamp(particle.velocity.x, -max_velocity, max_velocity);

            particle.position += particle.velocity * dt;

            // Reset force for the next iteration
            particle.force = {0, 0, 0};

            // Collisions
            m_boundaries.handle_collision(old_position, particle.position,
                                          particle.velocity);

            particle.position.z = m_grid.m_domain_limits.z / 2;
            particle.velocity.z = 0;

            if (!m_grid.within_domain_bounds(particle.position)) {
                throw std::out_of_range(
                    "Particle position is out of domain bounds");
            }

            // Update particle position in the grid
            auto new_grid_cell = m_grid.position_to_cell(particle.position);
            if (new_grid_cell != old_grid_cell) {
                m_grid.remove_particle(&particle, old_grid_cell);
                m_grid.add_particle(&particle);
            }
        }
    }

    void update(double dt) {
        update_forces();
        update_positions(dt);
    }

};
