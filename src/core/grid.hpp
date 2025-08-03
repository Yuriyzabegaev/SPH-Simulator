#pragma once
#include "particle.hpp"
#include "vec3.hpp"
#include <cassert>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#include "../core/print.hpp"

class Grid {
    std::vector<std::unordered_set<Particle *>> grid_;

    inline size_t idx_1d(const vec3<size_t> &idx_3d) const {
        return idx_1d(idx_3d.z, idx_3d.y, idx_3d.x);
    }

    inline size_t idx_1d(size_t z, size_t y, size_t x) const {
        return (z * grid_dims_.y + y) * grid_dims_.x + x;
    }

  public:
    const vec3<double> domain_limits_;
    const vec3<size_t> grid_dims_;
    const vec3<double> grid_cell_size_;
    Grid(vec3<size_t> grid_dims, vec3<double> domain_limits)
        : grid_([&grid_dims]() {
              size_t size_1d = grid_dims.z * grid_dims.y * grid_dims.x;
              return std::vector<std::unordered_set<Particle *>>(size_1d);
          }()),
          domain_limits_(std::move(domain_limits)),
          grid_dims_(std::move(grid_dims)), grid_cell_size_({
                                                 domain_limits.z / grid_dims.z,
                                                 domain_limits.y / grid_dims.y,
                                                 domain_limits.x / grid_dims.x,
                                             }) {}

    vec3<size_t> position_to_cell(const vec3<double> &position) const {
        return {
            static_cast<size_t>(position.z / grid_cell_size_.z),
            static_cast<size_t>(position.y / grid_cell_size_.y),
            static_cast<size_t>(position.x / grid_cell_size_.x),
        };
    }

    const std::unordered_set<Particle *> &
    particles_in_cell(const vec3<size_t> &cell) const {
        return grid_[idx_1d(cell)];
    }

    void add_particle(Particle *particle) {
        assert(within_domain_bounds(particle->position));
        grid_[idx_1d(position_to_cell(particle->position))].insert(particle);
    }

    void remove_particle(Particle *particle, const vec3<size_t> &grid_cell) {
        grid_[idx_1d(grid_cell)].erase(particle);
    }

    inline bool within_domain_bounds(const vec3<double> &pos) const {
        double tol = 1e-6; // Tolerance for floating-point comparison
        return (pos.z >= -tol && pos.z <= (domain_limits_.z + tol) &&
                pos.y >= -tol && pos.y <= (domain_limits_.y + tol) &&
                pos.x >= -tol && pos.x <= (domain_limits_.x + tol));
    }
};