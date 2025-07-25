#pragma once
#include "vec3.hpp"
#include <utility>
#include <cstdint>

struct Particle {
    vec3<double> position;
    double mass;   // kg
    double radius; // m
    double density = 0;
    vec3<double> force = {0, 0, 0};
    vec3<double> velocity = {0, 0, 0};
    uint8_t color[3] = {148, 160, 255};
    uint8_t opacity = 100;
    int generation = 0;

    Particle(vec3<double> position, double radius, double mass)
        : position(std::move(position)), mass(mass), radius(radius) {}
};