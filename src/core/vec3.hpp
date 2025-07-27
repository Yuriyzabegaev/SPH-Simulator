#pragma once

#include <cmath>

template <typename T> struct vec3 {
    T z, y, x;

    vec3() : z(0), y(0), x(0) {};
    vec3(T z, T y, T x) : z(z), y(y), x(x) {};

    double dot(const vec3 &other) const {
        return (z * other.z + y * other.y + x * other.x);
    }

    friend inline void operator*=(vec3<T> &a, const vec3<T> &b) {
        a.z *= b.z;
        a.y *= b.y;
        a.x *= b.x;
    }

    friend inline void operator*=(vec3<T> &a, const T &b) {
        a.z *= b;
        a.y *= b;
        a.x *= b;
    }

    friend inline void operator+=(vec3<T> &a, const vec3<T> &b) {
        a.z += b.z;
        a.y += b.y;
        a.x += b.x;
    }

    friend inline void operator+=(vec3<T> &a, const T &b) {
        a.z += b;
        a.y += b;
        a.x += b;
    }

    friend inline void operator-=(vec3<T> &a, const vec3<T> &b) {
        a.z -= b.z;
        a.y -= b.y;
        a.x -= b.x;
    }

    friend inline void operator-=(vec3<T> &a, const T &b) {
        a.z -= b;
        a.y -= b;
        a.x -= b;
    }

    friend inline void operator/=(vec3<T> &a, const vec3<T> &b) {
        a.z /= b.z;
        a.y /= b.y;
        a.x /= b.x;
    }

    friend inline void operator/=(vec3<T> &a, const T &b) {
        a.z /= b;
        a.y /= b;
        a.x /= b;
    }

    // Left operators: vec3 <op> T
    friend inline vec3 operator+(const vec3 &a, const T &b) {
        return vec3(a.z + b, a.y + b, a.x + b);
    }

    friend inline vec3 operator-(const vec3 &a, const T &b) {
        return vec3(a.z - b, a.y - b, a.x - b);
    }

    friend inline vec3 operator*(const vec3 &a, const T &b) {
        return vec3(a.z * b, a.y * b, a.x * b);
    }

    friend inline vec3 operator/(const vec3 &a, const T &b) {
        return vec3(a.z / b, a.y / b, a.x / b);
    }

    // Right operators: T <op> vec3
    friend inline vec3 operator+(const T &a, const vec3 &b) {
        return vec3(a + b.z, a + b.y, a + b.x);
    }

    friend inline vec3 operator-(const T &a, const vec3 &b) {
        return vec3(a - b.z, a - b.y, a - b.x);
    }

    friend inline vec3 operator*(const T &a, const vec3 &b) {
        return vec3(a * b.z, a * b.y, a * b.x);
    }

    friend inline vec3 operator/(const T &a, const vec3 &b) {
        return vec3(a / b.z, a / b.y, a / b.x);
    }

    // vec3 <op> vec3
    friend inline vec3 operator+(const vec3 &a, const vec3 &b) {
        return vec3(a.z + b.z, a.y + b.y, a.x + b.x);
    }

    friend inline vec3 operator-(const vec3 &a, const vec3 &b) {
        return vec3(a.z - b.z, a.y - b.y, a.x - b.x);
    }

    friend inline vec3 operator*(const vec3 &a, const vec3 &b) {
        return vec3(a.z * b.z, a.y * b.y, a.x * b.x);
    }

    friend inline vec3 operator/(const vec3 &a, const vec3 &b) {
        return vec3(a.z / b.z, a.y / b.y, a.x / b.x);
    }

    friend inline bool operator==(const vec3 &a, const vec3 &b) {
        return (a.z == b.z && a.y == b.y && a.x == b.x);
    }
};

template <typename T> bool isnan(const vec3<T> &x) {
    return std::isnan(x.x) || std::isnan(x.y) || std::isnan(x.z);
}