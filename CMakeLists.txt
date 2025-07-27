cmake_minimum_required(VERSION 3.10)
project(MyProject)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Always use g++
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra")

# Shared C++ sources
file(GLOB_RECURSE CORE_SOURCES src/core/*.cpp src/core/*.hpp)

# ========== Web build ==========
if(EMSCRIPTEN)
    message(STATUS "Configuring for WebAssembly build")

    file(GLOB WEB_BINDINGS src/web/*.cpp)

    add_library(web_module SHARED ${CORE_SOURCES} ${WEB_BINDINGS})

    set_target_properties(web_module PROPERTIES
        OUTPUT_NAME "module"
        SUFFIX ".js"
        LINK_FLAGS "--bind -s MODULARIZE=1 -s EXPORT_NAME='Module' -s WASM=1"
    )

# ========== Desktop build ==========
else()
    message(STATUS "Configuring for Desktop build")

    add_executable(app src/desktop/main.cpp ${CORE_SOURCES})

    find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)
    target_link_libraries(app sfml-graphics sfml-window sfml-system)
endif()
