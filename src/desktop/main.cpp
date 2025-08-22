#include "../core/grid.hpp"
#include "../core/renderer.hpp"
#include "../core/simulation.hpp"
#include "../core/initialize_simulation.hpp"
#include <algorithm>
#include <cmath>
#include <memory>

int main() {
    auto simulation = initialize_simulation();
    simulation->set_pressure_solver_type(1);
    SFMLRenderer renderer(std::move(simulation));
    renderer.run_until_complete();
}
