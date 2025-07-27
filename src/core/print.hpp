#pragma once
#include <iostream>
#include <vector>

template <typename T> void print(const T &input) {
    std::cout << input << std::endl;
}

template <typename T> void print(const std::vector<T> &input) {
    for (const auto &x : input)
        print(x);
    std::cout << std::endl;
}
