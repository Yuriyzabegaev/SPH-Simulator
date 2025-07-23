#pragma once
#include "particle.hpp"
#include "vec3.hpp"
#include <stdexcept>
#include <unordered_set>
#include <vector>

class Grid {
    std::vector<std::unordered_set<Particle *>> m_grid;

  public:
    const vec3<double> m_domain_limits;
    const vec3<size_t> m_grid_dims;
    const vec3<double> m_grid_cell_size;
    Grid(vec3<size_t> grid_dims, vec3<double> domain_limits)

        : m_grid([&grid_dims]() {
              size_t size_1d = grid_dims.z * grid_dims.y * grid_dims.x;
              return std::vector<std::unordered_set<Particle *>>(size_1d);
          }()),
          m_domain_limits(std::move(domain_limits)),
          m_grid_dims(std::move(grid_dims)), m_grid_cell_size({
                                                 domain_limits.z / grid_dims.z,
                                                 domain_limits.y / grid_dims.y,
                                                 domain_limits.x / grid_dims.x,
                                             }) {}

    vec3<size_t> position_to_cell(const vec3<double> &position) const {
        return {
            static_cast<size_t>(position.z / m_grid_cell_size.z),
            static_cast<size_t>(position.y / m_grid_cell_size.y),
            static_cast<size_t>(position.x / m_grid_cell_size.x),
        };
    }

    inline vec3<size_t> idx_3d(size_t idx_1d) const {
        return {
            idx_1d / (m_grid_dims.y * m_grid_dims.x),
            (idx_1d / m_grid_dims.x) % m_grid_dims.y,
            idx_1d % m_grid_dims.x,
        };
    }

    inline size_t idx_1d(const vec3<size_t> &idx_3d) const {
        return idx_1d(idx_3d.z, idx_3d.y, idx_3d.x);
    }

    inline size_t idx_1d(size_t z, size_t y, size_t x) const {
        return (z * m_grid_dims.y + y) * m_grid_dims.x + x;
    }

    const std::unordered_set<Particle *> &
    particles_in_cell(const vec3<size_t> &cell) const {
        return m_grid[idx_1d(cell)];
    }

    void add_particle(Particle *particle) {
        if (!within_domain_bounds(particle->position)) {
            throw std::out_of_range(
                "Particle position is out of domain bounds");
        }
        m_grid[idx_1d(position_to_cell(particle->position))].insert(particle);
    }

    void remove_particle(Particle *particle, const vec3<size_t> &grid_cell) {
        m_grid[idx_1d(grid_cell)].erase(particle);
    }

    inline bool within_domain_bounds(const vec3<double> &pos) const {
        double tol = 1e-6; // Tolerance for floating-point comparison
        return (pos.z >= -tol && pos.z <= (m_domain_limits.z + tol) &&
                pos.y >= -tol && pos.y <= (m_domain_limits.y + tol) &&
                pos.x >= -tol && pos.x <= (m_domain_limits.x + tol));
    }

    inline const std::vector<std::unordered_set<Particle *>> &grid() const {
        return m_grid;
    }
};