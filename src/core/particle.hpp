#pragma once
#include "vec3.hpp"
#include <cstdint>
#include <memory>
#include <utility>

constexpr double PARTICLE_RADIUS = 0.1;

struct Particle {
    vec3<double> position;
    double initial_density;
    double mass; // kg
    double density = 0;
    vec3<double> force = {0, 0, 0};
    vec3<double> velocity = {0, 0, 0};

    Particle(vec3<double> position, double mass, double initial_density)
        : position(std::move(position)), mass(mass),
          initial_density(initial_density) {}
};

Particle make_particle(vec3<double> position, double density) {
    double mass = density * 4 / 3 * M_PI * PARTICLE_RADIUS * PARTICLE_RADIUS *
                  PARTICLE_RADIUS;
    return Particle(std::move(position), mass, density);
}