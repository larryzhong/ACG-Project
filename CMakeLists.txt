cmake_minimum_required(VERSION 3.15)

project(CPUPathTracer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(pathtracer
    src/main.cpp
    src/io/gltf_loader.cpp
)
target_include_directories(pathtracer PRIVATE src)
