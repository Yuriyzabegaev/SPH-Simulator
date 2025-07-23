#pragma once
#include "vec3.hpp"
#include <cmath>

struct BoundaryPlane {
    vec3<double> normal;
    vec3<double> point;

    double signed_distance(const vec3<double> &pos) const {
        return normal.dot(pos - point);
    }
};

BoundaryPlane create_plane_3_points(const vec3<double> &p1,
                                    const vec3<double> &p2,
                                    const vec3<double> &p3) {
    vec3<double> v1 = p2 - p1;
    vec3<double> v2 = p3 - p1;
    vec3<double> normal = {v1.x * v2.y - v1.y * v2.x, v1.z * v2.x - v1.x * v2.z,
                           v1.y * v2.z - v1.z * v2.y};
    normal /= std::sqrt(normal.dot(normal));
    return {normal, p1};
}

class ExternalBoundaries {
    BoundaryPlane east;
    BoundaryPlane west;
    BoundaryPlane north;
    BoundaryPlane south;
    BoundaryPlane top;
    BoundaryPlane bottom;

  public:
    ExternalBoundaries(const vec3<double> &domain_limits)
        : east(create_plane_3_points({0, 0, domain_limits.x},
                                     {0, 1, domain_limits.x},
                                     {1, 0, domain_limits.x})),
          west(create_plane_3_points({0, 0, 0}, {1, 0, 0}, {0, 1, 0})),
          north(create_plane_3_points({0, domain_limits.y, 0},
                                      {1, domain_limits.y, 0},
                                      {0, domain_limits.y, 1})),
          south(create_plane_3_points({1, 0, 0}, {0, 0, 0}, {0, 0, 1})),
          top(create_plane_3_points({domain_limits.z, 0, 0},
                                    {domain_limits.z, 0, 1},
                                    {domain_limits.z, 1, 0})),
          bottom(create_plane_3_points({0, 0, 0}, {0, 1, 0}, {0, 0, 1})) {}

    void handle_collision(const vec3<double> &old_pos, vec3<double> &new_pos,
                          vec3<double> &vel) const {
        const BoundaryPlane *planes[] = {&east,  &west, &north,
                                         &south, &top,  &bottom};
        const double tol = 1e-6;

        static int count = 0;

        for (const BoundaryPlane *plane : planes) {
            double old_dst = plane->signed_distance(old_pos);
            double new_dst = plane->signed_distance(new_pos);
            if ((old_dst * new_dst) < 0) {
                count ++;
                new_pos -= plane->normal * new_dst;
                new_pos -= plane->normal * tol;
                // vel *= 0;
                vec3<double> vel_n = vel.dot(plane->normal) * plane->normal;
                vel -= 1.1 * vel_n;
                count++;
            }
        }
    }
};