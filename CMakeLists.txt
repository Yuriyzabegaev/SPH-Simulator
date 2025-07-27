cmake_minimum_required(VERSION 3.10)
project(MyProject)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Debug for desktop, optimized for web
if(EMSCRIPTEN)
    set(CMAKE_BUILD_TYPE Release)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
else()
    set(CMAKE_BUILD_TYPE Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# Shared C++ sources
file(GLOB_RECURSE CORE_SOURCES src/core/*.cpp src/core/*.hpp)

# ========== Web build ==========
if(EMSCRIPTEN)
    message(STATUS "Configuring for WebAssembly build")

    set(WEB_DIST_DIR ${CMAKE_SOURCE_DIR}/dist)
    file(MAKE_DIRECTORY ${WEB_DIST_DIR})

    file(GLOB WEB_BINDINGS src/web/*.cpp)
    add_executable(web_module ${WEB_BINDINGS})

    target_include_directories(web_module PRIVATE
        ${CMAKE_SOURCE_DIR}/dependencies/emsdk/upstream/emscripten/system/include
    )

    # Link flags for emscripten
    set_target_properties(web_module PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${WEB_DIST_DIR}
        OUTPUT_NAME "sim_module"
        SUFFIX ".js"
        LINK_FLAGS "--bind -s MODULARIZE=1 -s EXPORT_NAME='SimModule' -s ALLOW_MEMORY_GROWTH=1"
    )

    # Copy HTML, CSS, user JS into dist/
    file(GLOB STATIC_FILES "web/*.html" "web/*.css" "web/*.js")
    foreach(f IN LISTS STATIC_FILES)
        add_custom_command(
            TARGET web_module POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${f} ${WEB_DIST_DIR}
        )
    endforeach()

# ========== Desktop build ==========
else()
    message(STATUS "Configuring for Desktop build")

    add_executable(app src/desktop/main.cpp ${CORE_SOURCES})

    find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)
    target_link_libraries(app sfml-graphics sfml-window sfml-system)
endif()
