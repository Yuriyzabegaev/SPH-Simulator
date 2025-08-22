#pragma once
#include "../core/print.hpp"
#include "particle.hpp"
#include "vec3.hpp"
#include <cassert>
#include <functional>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class Grid {
    std::vector<std::unordered_set<Particle *>> grid_;

    inline size_t idx_1d(const vec3<size_t> &idx_3d) const {
        return idx_1d(idx_3d.z, idx_3d.y, idx_3d.x);
    }

    inline size_t idx_1d(size_t z, size_t y, size_t x) const {
        z = std::clamp(z, static_cast<size_t>(0), grid_dims_.z);
        y = std::clamp(y, static_cast<size_t>(0), grid_dims_.y);
        x = std::clamp(x, static_cast<size_t>(0), grid_dims_.x);
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
        const double z = std::clamp(position.z, 0.,
                                    domain_limits_.z - grid_cell_size_.z / 4);
        const double y = std::clamp(position.y, 0.,
                                    domain_limits_.y - grid_cell_size_.y / 4);
        const double x = std::clamp(position.x, 0.,
                                    domain_limits_.x - grid_cell_size_.x / 4);
        return {
            static_cast<size_t>(z / grid_cell_size_.z),
            static_cast<size_t>(y / grid_cell_size_.y),
            static_cast<size_t>(x / grid_cell_size_.x),
        };
    }

    const std::unordered_set<Particle *> &
    particles_in_cell(const vec3<size_t> &cell) const {
        return grid_[idx_1d(cell)];
    }

    void add_particle(Particle *particle, const vec3<double> &position) {
        // assert(within_domain_bounds(position));
        size_t idx = idx_1d(position_to_cell(position));
        if (idx >= grid_.size()) {
            // Handle error: out-of-bounds particle position
            throw std::out_of_range("Particle position outside grid domain");
        }
        grid_[idx].insert(particle);
    }

    void remove_particle(Particle *particle, const vec3<size_t> &grid_cell) {
        auto &data = grid_[idx_1d(grid_cell)];
        assert(data.contains(particle));
        data.erase(particle);
    }

    inline bool within_domain_bounds(const vec3<double> &pos) const {
        double tol = 1e-6; // Tolerance for floating-point comparison
        return (pos.z >= -tol && pos.z <= (domain_limits_.z + tol) &&
                pos.y >= -tol && pos.y <= (domain_limits_.y + tol) &&
                pos.x >= -tol && pos.x <= (domain_limits_.x + tol));
    }
};

static constexpr int8_t OFFESTS_3X3[27][3] = {
    {-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1},
    {-1, 1, -1},  {-1, 1, 0},  {-1, 1, 1},  {0, -1, -1}, {0, -1, 0}, {0, -1, 1},
    {0, 0, -1},   {0, 0, 0},   {0, 0, 1},   {0, 1, -1},  {0, 1, 0},  {0, 1, 1},
    {1, -1, -1},  {1, -1, 0},  {1, -1, 1},  {1, 0, -1},  {1, 0, 0},  {1, 0, 1},
    {1, 1, -1},   {1, 1, 0},   {1, 1, 1}};

/**
 * @brief Calls the inner function for each pair of particles within the same
 * and neighboring cells. Each pair is guaranteed to be called exactly once.
 */
void for_pair_neighbor_cells(
    const Grid &grid, std::function<void(Particle *, Particle *)> inner) {
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

/**
 * @brief Calls the inner function for each particle within the same and
 * neighboring cells to the provided particle p1, excluding itself. Each pair is
 * guaranteed to be called exactly once.
 */
void for_particle_neighbors(const Grid &grid, Particle *p1,
                            std::function<void(Particle *)> inner) {
    auto [i, j, k] = grid.position_to_cell(p1->position);
    for (auto [di, dj, dk] : OFFESTS_3X3) {
        if (((i == 0) && (di < 0)) || ((j == 0) && (dj < 0)) ||
            ((k == 0) && (dk < 0)))
            continue;
        size_t i1 = i + di, j1 = j + dj, k1 = k + dk;
        if ((i1 >= grid.grid_dims_.z) || (j1 >= grid.grid_dims_.y) ||
            (k1 >= grid.grid_dims_.x))
            continue;

        auto &particles_cell_2 = grid.particles_in_cell({i1, j1, k1});
        for (auto p2 : particles_cell_2) {
            if (p2 != p1)
                inner(p2);
        }
    }
}
