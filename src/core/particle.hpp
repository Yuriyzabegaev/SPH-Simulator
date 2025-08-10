#pragma once
#include "vec3.hpp"
#include <cstdint>
#include <memory>
#include <utility>

constexpr double PARTICLE_RADIUS = 0.3;
constexpr double MASS = 56.6;
constexpr double TARGET_DENSITY = 300;
struct Particle {
    vec3<double> position;
    vec3<double> predicted_position = {0, 0, 0};
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

    // return Particle(std::move(position), mass, density);
    return Particle(std::move(position), MASS, TARGET_DENSITY);

}