#include "../core/grid.hpp"
#include "../core/renderer.hpp"
#include "../core/simulation.hpp"
#include "../core/initialize_simulation.hpp"
#include <algorithm>
#include <cmath>
#include <memory>

int main() {
    auto simulation = initialize_simulation();
    SFMLRenderer renderer(std::move(simulation));
    renderer.run_until_complete();
}
