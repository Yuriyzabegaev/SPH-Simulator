#pragma once
#include "grid.hpp"
#include "simulation.hpp"
#include <SFML/Graphics.hpp>
#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <thread>
#include <vector>

static constexpr int target_frame_ms = 1000 / 60; // ~16ms per frame for 60 FPS

std::array<int, 3> blue_to_red_gradient(float value, float max_value) {
    float t =
        value <= 0 ? 0.0f : (value >= max_value ? 1.0f : value / max_value);
    int r = static_cast<int>(255 * t);
    int g = 0;
    int b = static_cast<int>(255 * (1.0f - t));
    return {r, g, b};
}

class SFMLRenderer {
    float m_scale_sim_to_window;
    bool m_is_window_closed = false;
    bool m_is_dragging = false;
    int m_mouse_x = -1;
    int m_mouse_y = -1;
    const double m_dt =
        static_cast<double>(target_frame_ms) / 1000; // Time step
    float m_drag_radius = 100.f;
    Simulation m_sim;
    sf::RenderWindow m_window{sf::VideoMode(800, 600), "Particles"};

    inline std::array<float, 2> tranform_to_simulation_coords(int x, int y) {
        float sim_x = x * m_sim.m_grid.m_domain_limits.x / 800.f;
        float sim_y = (600.f - y) * m_sim.m_grid.m_domain_limits.y / 600.f;
        return {sim_x, sim_y};
    }

    void handle_start_drag(const sf::Event &event) {
        m_is_dragging = true;
        auto pos = m_window.mapPixelToCoords(
            {event.mouseButton.x, event.mouseButton.y});
        m_mouse_x = pos.x;
        m_mouse_y = pos.y;
    }

    void handle_end_drag(const sf::Event &event) { m_is_dragging = false; }

    void handle_move_drag(const sf::Event &event) {
        assert(m_is_dragging);
        auto pos =
            m_window.mapPixelToCoords({event.mouseMove.x, event.mouseMove.y});
        auto [sim_x, sim_y] = tranform_to_simulation_coords(pos.x, pos.y);
        auto [prev_x, prev_y] =
            tranform_to_simulation_coords(m_mouse_x, m_mouse_y);
        vec3<double> v{0, sim_y - prev_y, sim_x - prev_x};
        v /= m_dt * m_dt * 10;
        m_sim.apply_external_force(
            {m_sim.m_grid.m_domain_limits.z / 2, sim_y, sim_x}, v,
            m_drag_radius / m_scale_sim_to_window);

        m_mouse_x = event.mouseMove.x;
        m_mouse_y = event.mouseMove.y;
    }

    void handle_add_particle(const sf::Event &event) {
        auto pos = m_window.mapPixelToCoords(
            {event.mouseButton.x, event.mouseButton.y});
        // Reverse transformation: window coordinates to simulation
        // coordinates
        auto [sim_x, sim_y] = tranform_to_simulation_coords(pos.x, pos.y);
        m_sim.add_particle(Particle({0.25, sim_y, sim_x}, 0.1, 1000));
    }

    void handle_events() {
        sf::Event event;
        while (m_window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                m_window.close();
                m_is_window_closed = true;
                break;
            } else if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Button::Right)
                    handle_add_particle(event);
                else if (event.mouseButton.button == sf::Mouse::Button::Left)
                    handle_start_drag(event);
            } else if ((event.type == sf::Event::MouseMoved) && m_is_dragging) {
                handle_move_drag(event);
            } else if ((event.type == sf::Event::MouseButtonReleased) &&
                       m_is_dragging &&
                       event.mouseButton.button == sf::Mouse::Button::Left) {
                handle_end_drag(event);
            } else if (event.type == sf::Event::LostFocus) {
                handle_end_drag(event);
            }
        }
    }

    bool render() {
        m_window.clear(sf::Color::White);

        for (const auto &particle : m_sim.m_particles) {
            // Convert particle position to window coordinates
            float x =
                particle->position.x * 800 / m_sim.m_grid.m_domain_limits.x;
            float y = 600 - (particle->position.y * 600 /
                             m_sim.m_grid.m_domain_limits.y);
            x = std::clamp(x, 0.f, 800.f);
            y = std::clamp(y, 0.f, 600.f);

            float radius = particle->radius * m_scale_sim_to_window;

            sf::CircleShape shape(radius / 2);
            shape.setOrigin(radius, radius);
            shape.setPosition(x, y);
            auto color = blue_to_red_gradient(
                std::sqrt(particle->velocity.dot(particle->velocity)), 10.);
            // auto color = blue_to_red_gradient(particle->density - 2000, 500);

            shape.setFillColor(sf::Color(color[0], color[1], color[2], 100));
            m_window.draw(shape);
        }

        if (m_is_dragging) {
            sf::CircleShape drag_circle(m_drag_radius);
            drag_circle.setOrigin(m_drag_radius, m_drag_radius);
            drag_circle.setPosition({m_mouse_x, m_mouse_y});
            drag_circle.setOutlineColor(sf::Color::Black);
            drag_circle.setFillColor(sf::Color::Transparent);
            drag_circle.setOutlineThickness(1);
            m_window.draw(drag_circle);
        }

        handle_events();
        m_window.display();
    }

  public:
    SFMLRenderer(Simulation sim) : m_sim(std::move(sim)) {
        auto scale_x = 800. / m_sim.m_grid.m_domain_limits.x;
        auto scale_y = 600. / m_sim.m_grid.m_domain_limits.y;
        m_scale_sim_to_window = std::max(scale_y, scale_x);
    }

    void run_until_complete() {

        while (!m_is_window_closed) {
            auto frame_start = std::chrono::steady_clock::now();
            render();
            m_sim.update(m_dt);
            auto frame_end = std::chrono::steady_clock::now();
            auto frame_duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    frame_end - frame_start)
                    .count();
            std::cout << "Frame duration: " << frame_duration << std::endl;
            if (frame_duration < target_frame_ms) {
                std::this_thread::sleep_for(std::chrono::milliseconds(
                    target_frame_ms - frame_duration));
            }
        }
    }
};