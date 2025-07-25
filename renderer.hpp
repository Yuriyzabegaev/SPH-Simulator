#pragma once
#include "grid.hpp"
#include "simulation.hpp"
#include <SFML/Graphics.hpp>
#include <chrono>
#include <iostream>
#include <memory>
#include <thread>
#include <vector>

static constexpr int target_frame_ms = 1000 / 60; // ~16ms per frame for 60 FPS

class SFMLRenderer {
    float m_scale_x;
    float m_scale_y;
    float m_fps = -1;

  public:
    Simulation m_sim;
    sf::RenderWindow m_window{sf::VideoMode(800, 600), "Particles"};
    SFMLRenderer(Simulation sim) : m_sim(std::move(sim)) {
        m_scale_x = 800. / m_sim.m_grid.m_domain_limits.x;
        m_scale_y = 600. / m_sim.m_grid.m_domain_limits.y;
    }

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

            float radius = particle.radius * std::max(m_scale_y, m_scale_x);

            sf::CircleShape shape(radius);
            shape.setOrigin(radius, radius);
            shape.setPosition(x, y);
            shape.setOutlineColor(sf::Color::Blue);
            shape.setFillColor(sf::Color(particle.color[0], particle.color[1],
                                         particle.color[2], particle.opacity));
            m_window.draw(shape);
        }

        // static sf::Font font;
        // static bool fontLoaded = false;
        // if (!fontLoaded) {
        //     if (!font.loadFromFile("arial.ttf")) {
        //     std::cerr << "Failed to load font\n";
        //     } else {
        //     fontLoaded = true;
        //     }
        // }
        // if (fontLoaded) {
        // sf::Text fpsText;
        // fpsText.setFont(font);
        // fpsText.setString("FPS: " + std::to_string(static_cast<int>(m_fps)));
        // fpsText.setCharacterSize(40);
        // fpsText.setFillColor(sf::Color::Red);
        // fpsText.setPosition(400.f, 400.f);
        // m_window.draw(fpsText);
        // }
        std::cout << "FPS: " << m_fps << std::endl;

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
            auto frame_start = std::chrono::steady_clock::now();
            if (!render()) {
                break;
            }
            m_sim.update(dt);
            auto frame_end = std::chrono::steady_clock::now();
            auto frame_duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    frame_end - frame_start)
                    .count();
            m_fps = 1000. / frame_duration;
            if (frame_duration < target_frame_ms) {
                std::this_thread::sleep_for(std::chrono::milliseconds(
                    target_frame_ms - frame_duration));
            }
        }
    }
};