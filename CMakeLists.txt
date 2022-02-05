cmake_minimum_required(VERSION 3.19)
project(dev_fdaPDE)

set(CMAKE_CXX_STANDARD 17)

add_executable(dev_fdaPDE space_time.cpp)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-O3")