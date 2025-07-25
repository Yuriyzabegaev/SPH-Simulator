#pragma once
#include "vec3.hpp"
#include <cstdint>
#include <utility>

struct Particle {
    vec3<double> position;
    double mass;   // kg
    double radius; // m
    double density = 0;
    vec3<double> force = {0, 0, 0};
    vec3<double> velocity = {0, 0, 0};

    Particle(vec3<double> position, double radius, double density)
        : position(std::move(position)), mass([radius, density]() {
              return density * 4 / 3 * M_PI * radius * radius * radius;
          }()),
          radius(radius) {}
};