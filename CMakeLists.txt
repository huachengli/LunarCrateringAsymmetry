cmake_minimum_required(VERSION 3.18)
project(frequency_crater C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -std=c99 -O0")
add_executable(frequency_global main.c integrate_base.c integrate_model.c)
target_link_libraries(frequency_global -lm)
