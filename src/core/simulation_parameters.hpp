#pragma once

struct SimulationParameters {
    double viscosity = 0.01;
    double gravity = 10;
    double specific_volume = 2e6;
    double target_density = 300;
};
const double RHO_0 = -1.;