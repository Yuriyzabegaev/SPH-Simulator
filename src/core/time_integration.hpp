#pragma once
#include "particle.hpp"
#include "vec3.hpp"

static inline vec3<double> a(const Particle *p) { return p->force / p->mass; }

static void clamp_velosity(Particle *p, double vmax) {
    p->velocity.x = std::clamp(p->velocity.x, -vmax, vmax);
    p->velocity.y = std::clamp(p->velocity.y, -vmax, vmax);
    p->velocity.z = std::clamp(p->velocity.z, -vmax, vmax);
}

void integrate_semi_implicit_euler(Particle *p, double dt, double vmax) {
    p->velocity += a(p) * dt;

    // CFL limiter
    clamp_velosity(p, vmax);

    p->position += p->velocity * dt;
}

void integrateRK4(Particle *p, double dt, double vmax) {
    auto k1_v = a(p);
    auto k1_x = p->velocity;

    auto k2_v = a(p); // force assumed constant over dt
    auto k2_x = p->velocity + 0.5 * dt * k1_v;

    auto k3_v = a(p);
    auto k3_x = p->velocity + 0.5 * dt * k2_v;

    auto k4_v = a(p);
    auto k4_x = p->velocity + dt * k3_v;

    p->velocity += (dt / 6.0) * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v);

    // CFL limiter
    clamp_velosity(p, vmax);

    p->position += (dt / 6.0) * (k1_x + 2.0 * k2_x + 2.0 * k3_x + k4_x);
}