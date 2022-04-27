cmake_minimum_required(VERSION 3.21)
project(network_reliability_unreliability)

set(CMAKE_CXX_STANDARD 17)

add_executable(network_reliability_unreliability reliability.hpp unreliability.hpp graph.hpp median_trick.hpp
        random.hpp main.cpp brute.hpp utils.hpp types.hpp graph_generator.hpp mincut.hpp digraph.hpp mincut_linear.hpp)
