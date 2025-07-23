#pragma once
#include "grid.hpp"
#include "simulation.hpp"
#include <SFML/Graphics.hpp>
#include <iostream>
#include <memory>
#include <vector>

static constexpr int target_frame_ms = 1000 / 60; // ~16ms per frame for 60 FPS

class SFMLRenderer {

  public:
    Simulation m_sim;
    sf::RenderWindow m_window{sf::VideoMode(800, 600), "Particles"};
    SFMLRenderer(Simulation sim) : m_sim(std::move(sim)) {}

    bool render() {
        m_window.clear(sf::Color::White);

        for (const auto &particle : m_sim.m_particles) {
            // Convert particle position to window coordinates
            float x =
                particle.position.x * 800 / m_sim.m_grid.m_domain_limits.x;
            float y = 600 - (particle.position.y * 600 /
                             m_sim.m_grid.m_domain_limits.y);
            x = std::clamp(x, 0.f, 800.f);
            y = std::clamp(y, 0.f, 600.f);
            sf::CircleShape shape(50);
            shape.setOrigin(50, 50);
            shape.setPosition(x, y);
            shape.setOutlineColor(sf::Color::Blue);
            shape.setFillColor(sf::Color(148, 160, 255, 100));
            m_window.draw(shape);
        }

        m_window.display();
        sf::Event event;
        while (m_window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                m_window.close();
                return false;
            } else if (event.type == sf::Event::MouseButtonReleased) {
                // auto pos = m_window.mapPixelToCoords(
                //     {event.mouseButton.x, event.mouseButton.y});
                // // Reverse transformation: window coordinates to simulation
                // // coordinates
                // float sim_x = pos.x * m_sim.m_grid.m_domain_limits.x / 800.f;
                // float sim_y =
                //     (600.f - pos.y) * m_sim.m_grid.m_domain_limits.y / 600.f;
                // m_sim.add_particle({{0.1, sim_y, sim_x}});
            }
        }
        return true;
    }

    void run_until_complete() {
        constexpr double dt =
            static_cast<double>(target_frame_ms) / 1000; // Time step
        for (;;) {
            if (!render()) {
                break;
            }
            m_sim.update(dt);
        }
    }
};