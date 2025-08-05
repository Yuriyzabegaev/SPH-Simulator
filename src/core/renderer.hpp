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
    float scale_sim_to_window_;
    bool is_window_closed_ = false;
    bool is_dragging_ = false;
    int mouse_x_ = -1;
    int mouse_y_ = -1;
    const double dt_ =
        static_cast<double>(target_frame_ms) / 1000; // Time step
    float drag_radius_ = 100.f;
    std::shared_ptr<Simulation> sim_;
    sf::RenderWindow window_{sf::VideoMode(800, 600), "Particles"};

    inline std::array<float, 2> tranform_to_simulation_coords(int x, int y) {
        float sim_x = x * sim_->grid_.domain_limits_.x / 800.f;
        float sim_y = (600.f - y) * sim_->grid_.domain_limits_.y / 600.f;
        return {sim_x, sim_y};
    }

    void handle_start_drag(const sf::Event &event) {
        is_dragging_ = true;
        auto pos = window_.mapPixelToCoords(
            {event.mouseButton.x, event.mouseButton.y});
        mouse_x_ = pos.x;
        mouse_y_ = pos.y;
    }

    void handle_end_drag(const sf::Event &event) {
        (void)event;
        is_dragging_ = false;
    }

    void handle_move_drag(const sf::Event &event) {
        assert(is_dragging_);
        auto pos =
            window_.mapPixelToCoords({event.mouseMove.x, event.mouseMove.y});
        mouse_x_ = pos.x;
        mouse_y_ = pos.y;
    }

    void handle_add_particle(const sf::Event &event) {
        auto pos = window_.mapPixelToCoords(
            {event.mouseButton.x, event.mouseButton.y});
        // Reverse transformation: window coordinates to simulation
        // coordinates
        auto [sim_x, sim_y] = tranform_to_simulation_coords(pos.x, pos.y);
        sim_->add_particle({sim_->grid_.domain_limits_.z / 2, sim_y, sim_x}, RHO_0);
    }

    void handle_events() {
        sf::Event event;
        while (window_.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window_.close();
                is_window_closed_ = true;
                break;
            } else if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Button::Right)
                    handle_add_particle(event);
                else if (event.mouseButton.button == sf::Mouse::Button::Left)
                    handle_start_drag(event);
            } else if ((event.type == sf::Event::MouseMoved) && is_dragging_) {
                handle_move_drag(event);
            } else if ((event.type == sf::Event::MouseButtonReleased) &&
                       is_dragging_ &&
                       event.mouseButton.button == sf::Mouse::Button::Left) {
                handle_end_drag(event);
            } else if (event.type == sf::Event::LostFocus) {
                handle_end_drag(event);
            }
        }
    }

    void apply_central_force(double value) {
        auto [sim_x, sim_y] =
            tranform_to_simulation_coords(mouse_x_, mouse_y_);
        sim_->apply_central_force(
            {sim_->grid_.domain_limits_.z / 2, sim_y, sim_x}, value,
            drag_radius_ / scale_sim_to_window_);
    }

    void render() {
        window_.clear(sf::Color::White);

        for (const auto &particle : sim_->particles_) {
            // Convert particle position to window coordinates
            float x =
                particle->position.x * 800 / sim_->grid_.domain_limits_.x;
            float y = 600 - (particle->position.y * 600 /
                             sim_->grid_.domain_limits_.y);
            x = std::clamp(x, 0.f, 800.f);
            y = std::clamp(y, 0.f, 600.f);

            float radius = 5.f;

            sf::CircleShape shape(radius);
            shape.setOrigin(radius, radius);
            shape.setPosition(x, y);
            auto color = blue_to_red_gradient(
                std::sqrt(particle->velocity.dot(particle->velocity)), 10.);
            // auto color = blue_to_red_gradient(particle->density - 2000, 500);

            shape.setFillColor(sf::Color(color[0], color[1], color[2], 100));
            window_.draw(shape);
        }

        if (is_dragging_) {
            apply_central_force(30.);
            sf::CircleShape drag_circle(drag_radius_);
            drag_circle.setOrigin(drag_radius_, drag_radius_);
            drag_circle.setPosition(
                {static_cast<float>(mouse_x_), static_cast<float>(mouse_y_)});
            drag_circle.setOutlineColor(sf::Color::Black);
            drag_circle.setFillColor(sf::Color::Transparent);
            drag_circle.setOutlineThickness(1);
            window_.draw(drag_circle);
        }

        handle_events();
        window_.display();
    }

  public:
    SFMLRenderer(std::shared_ptr<Simulation> sim) : sim_(sim) {
        auto scale_x = 800. / sim->grid_.domain_limits_.x;
        auto scale_y = 600. / sim->grid_.domain_limits_.y;
        scale_sim_to_window_ = std::max(scale_y, scale_x);
    }

    void run_until_complete() {

        while (window_.isOpen()) {
            auto frame_start = std::chrono::steady_clock::now();
            render();
            sim_->update(dt_);
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