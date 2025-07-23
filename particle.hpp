#pragma once
#include "vec3.hpp"
#include <utility>

struct Particle {
    vec3<double> position;
    double mass = 0.001; // 1000 for water, particle occupies 0.01 m^3
    vec3<double> force = {0, 0, 0};
    vec3<double> velocity = {0, 0, 0};

    Particle(vec3<double> position) : position(std::move(position)) {}
};