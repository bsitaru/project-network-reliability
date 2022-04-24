#ifndef NETWORK_RELIABILITY_UNRELIABILITY_RANDOM_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_RANDOM_HPP

#include <bits/stdc++.h>
#include "types.hpp"

using namespace std;

namespace Random {
    static mt19937 gen;

    static void init_predictable(bool predictable) {
        auto seed = predictable ? 151917 : chrono::high_resolution_clock::now().time_since_epoch().count();
        gen = mt19937(seed);
    }

    static void init(int seed) {
        gen = mt19937(seed);
    }

    static bool get_ber(t_double p) {
        return bernoulli_distribution(p)(gen);
    }

    static int get_int(int l, int r) {
        return uniform_int_distribution<int>(l, r)(gen);
    }
};

#endif //NETWORK_RELIABILITY_UNRELIABILITY_RANDOM_HPP
