#pragma once
#include "grid.hpp"
#include "kernel_function.hpp"
#include "particle.hpp"
#include "simulation_parameters.hpp"
#include "vec3.hpp"
#include <cassert>

static vec3<double> compute_pressure_acceleration(Particle *p1, Particle *p2,
                                                  double pressure_1,
                                                  double pressure_2) {
    vec3<double> r = p2->position - p1->position;
    double r_norm = std::sqrt(r.dot(r));
    const double tol = 1e-9;
    if (r_norm < tol) {
        return {0., 0., 0.};
    }
    double gradW_ = gradW(r_norm, 4 * PARTICLE_RADIUS);
    double tmp = pressure_1 / (p1->density * p1->density) +
                 pressure_2 / (p1->density * p2->density);
    tmp *= gradW_;
    auto result = tmp * r;
    assert(!isnan(result));
    return result;
}

class PressureSolver {

  public:
    virtual ~PressureSolver() {}
    virtual void update_pressure(double dt) = 0;
};

class PressureSolverSESPH : public PressureSolver {

    SimulationParameters *parameters_;
    Grid *grid_;

  public:
    PressureSolverSESPH() = delete;

    explicit PressureSolverSESPH(SimulationParameters *params, Grid *grid)
        : parameters_(params), grid_(grid) {}

    void update_pressure(double dt) override {
        for_pair_neighbor_cells(*grid_, [&](Particle *p1, Particle *p2) {
            update_pressure_force(p1, p2);
        });
    }

  private:
    double compute_pressure(const Particle *p, double k,
                            double target_density) const {
        double density_diff =
            std::max(p->density - target_density, 0.) / target_density;
        // double density_diff = (p->density - target_density) / target_density;
        return k * density_diff;
    }

    void update_pressure_force(Particle *p1, Particle *p2) const {
        double pressure_1 = compute_pressure(p1, parameters_->specific_volume,
                                             parameters_->target_density);
        double pressure_2 = compute_pressure(p2, parameters_->specific_volume,
                                             parameters_->target_density);
        auto acc =
            compute_pressure_acceleration(p1, p2, pressure_1, pressure_2);
        p1->force += p2->mass * acc;
        p2->force -= p1->mass * acc;
    }
};

class PressureSolverIISPH : public PressureSolver {

    SimulationParameters *parameters_;
    Grid *grid_;
    std::vector<std::unique_ptr<Particle>> *particles_;
    std::vector<double> a11_;
    std::vector<double> rhs_;
    std::vector<double> A_dot_p_;

  public:
    PressureSolverIISPH() = delete;

    explicit PressureSolverIISPH(
        SimulationParameters *params, Grid *grid,
        std::vector<std::unique_ptr<Particle>> *particles)
        : parameters_(params), grid_(grid), particles_(particles) {}

    void update_pressure(double dt) override {
        update_a11_rhs(dt);
        for (int i = 0; i < 10; ++i) {
            update_pressure_accelerations();
            update_matvec_operator(dt);
            auto residual = do_jacobi_step();
            // print(residual);
        }
        update_pressure_forces();
    }

  private:
    void update_a11_rhs(double dt) {
        a11_.resize(particles_->size(), 0);
        rhs_.resize(particles_->size(), 0);

        double dt_square = dt * dt;

        for (size_t i = 0; i < particles_->size(); ++i) {
            auto &p1 = (*particles_)[i];
            double rho_m2 = 1 / (p1->density * p1->density);
            rhs_[i] = p1->initial_density - p1->density;

            // Setting initial guess for pressure.
            p1->pressure = 0;

            for_particle_neighbors(*grid_, p1.get(), [&](Particle *p2) {
                vec3<double> r = p1->position - p2->position;
                double r_norm = std::sqrt(r.dot(r));
                double gradW_ = gradW(r_norm, 4 * PARTICLE_RADIUS);

                // Formula: a_ii = -Δt^2 Σ_j (m_j m_i / ρ^2 ∇W_ij · ∇W_ij)
                a11_[i] -=
                    dt_square * p2->mass * p1->mass * rho_m2 * gradW_ * gradW_;

                // Formula: s_i = ρ0 - ρ_i - Δt Σ_j (m_j (v_i - v_j) ·
                // ∇W_ij)
                rhs_[i] -=
                    dt * p2->mass *
                    (p1->velocity - p2->velocity).dot(r / r_norm * gradW_);
            });
        }
    }

    void update_pressure_accelerations() {
        for (auto &p : *particles_) {
            p->acceleration = {0, 0, 0};
        }

        for_pair_neighbor_cells(*grid_, [&](Particle *p1, Particle *p2) {
            auto acc = compute_pressure_acceleration(p1, p2, p1->pressure,
                                                     p2->pressure);
            // Using force to temporarily store accelerations.
            p1->acceleration += acc;
            p2->acceleration -= acc;
        });
    }

    void update_matvec_operator(double dt) {
        A_dot_p_.assign(particles_->size(), 0);
        for (size_t i = 0; i < particles_->size(); ++i) {
            auto &p1 = (*particles_)[i];
            for_particle_neighbors(*grid_, p1.get(), [&](Particle *p2) {
                vec3<double> r = p1->position - p2->position;
                double r_norm = std::sqrt(r.dot(r));
                double gradW_ = gradW(r_norm, 4 * PARTICLE_RADIUS);
                // Formula: (A·p)_i = Δt^2 Σ_j (a_i - a_j) · ∇W_ij
                A_dot_p_[i] += (p1->acceleration - p2->acceleration).dot(r / r_norm * gradW_);
            });
            A_dot_p_[i] *= dt * dt;
        }
    }

    double do_jacobi_step() {
        const double w = 0.5;
        double max_residual = 0;
        for (size_t i = 0; i < particles_->size(); ++i) {
            auto &p1 = (*particles_)[i];
            double residual = rhs_[i] - A_dot_p_[i];
            p1->pressure += w / a11_[i] * residual;
            max_residual = std::max(std::abs(residual), max_residual);
        }
        return max_residual;
    }

    void update_pressure_forces() const {
        for (auto &p : *particles_) {
            p->force = {0, 0, 0};
        }

        for_pair_neighbor_cells(*grid_, [&](Particle *p1, Particle *p2) {
            auto acc = compute_pressure_acceleration(
                p1, p2, std::max(p1->pressure, 0.), std::max(p2->pressure, 0.));
            p1->force += p2->mass * acc;
            p2->force -= p1->mass * acc;
        });
    }
};