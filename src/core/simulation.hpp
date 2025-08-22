#pragma once
#include "boundary.hpp"
#include "grid.hpp"
#include "kernel_function.hpp"
#include "pressure_solver.hpp"
#include "print.hpp"
#include "simulation_parameters.hpp"
#include "time_integration.hpp"
#include "vec3.hpp"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <thread>

class Simulation {

    ExternalBoundaries boundaries_;
    SimulationParameters parameters_;
    std::unique_ptr<PressureSolver> pressure_solver_ =
        std::make_unique<PressureSolverSESPH>(&parameters_, &grid_);

    void update_density_pair_cells(Particle *p1, Particle *p2) const {
        vec3<double> r = p2->position - p1->position;
        double r_norm = std::sqrt(r.dot(r));
        double W_ = W(r_norm, 4 * PARTICLE_RADIUS);
        p1->density += p2->mass * W_;
        p2->density += p1->mass * W_;
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

            auto old_grid_cell =
                grid_.position_to_cell(particle->predicted_position);
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

        grid_.add_particle(particles_.back().get(),
                           particles_.back()->position);
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
            auto dstSqr = r.dot(r);
            if (dstSqr >= (radius * radius)) {
                continue;
            }
            particle->force += particle->mass * acceleration * r;
        }
    }

    void apply_rotor_force(vec3<double> center, double acceleration,
                           double radius) {

        for (auto &particle : particles_) {
            auto r = center - particle->position;
            auto dstSqr = r.dot(r);
            if (dstSqr >= (radius * radius)) {
                continue;
            }

            vec3<double> rotor_dir{0, r.x, -r.y};
            rotor_dir /= std::sqrt(rotor_dir.dot(rotor_dir));
            particle->force += particle->mass * acceleration * rotor_dir;
        }
    }

    void set_viscosity(double viscosity) { parameters_.viscosity = viscosity; }

    void set_gravity(double gravity) { parameters_.gravity = gravity; }

    void set_specific_volume(double specific_volume) {
        parameters_.specific_volume = specific_volume;
    }

    void set_target_density(double target_density) {
        parameters_.target_density = target_density;
    }

    void set_pressure_solver_type(int type) {
        if (type == 0)
            pressure_solver_.reset(
                new PressureSolverSESPH(&parameters_, &grid_));
        else if (type == 1) {
            pressure_solver_.reset(
                new PressureSolverIISPH(&parameters_, &grid_, &particles_));
        } else {
            throw std::out_of_range("Use type 0 for PressureSolverSESPH and "
                                    "type 1 for PressureSolverIISPH.");
        }
    }

    void update(double dt) {
        update_densities();
        for_pair_neighbor_cells(grid_, [&](Particle *p1, Particle *p2) {
            update_viscous_force(p1, p2);
        });
        for (auto &p : particles_) {
            update_gravity_force(p.get());

            auto old_grid_cell = grid_.position_to_cell(p->position);
            p->velocity += dt * a(p.get());
            double artificial_dt = 1 / 60.;

            p->predicted_position = p->position + p->velocity * artificial_dt;

            p->force = {0, 0, 0};

            auto predicted_grid_cell =
                grid_.position_to_cell(p->predicted_position);
            if (predicted_grid_cell != old_grid_cell) {
                grid_.remove_particle(p.get(), old_grid_cell);
                grid_.add_particle(p.get(), p->predicted_position);
            }
        }

        pressure_solver_->update_pressure(dt);
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
