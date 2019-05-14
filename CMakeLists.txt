cmake_minimum_required(VERSION 3.12)
project(delaunay)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")

add_executable(delaunay main.cpp)