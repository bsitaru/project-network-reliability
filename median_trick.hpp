#ifndef NETWORK_RELIABILITY_UNRELIABILITY_MEDIAN_TRICK_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_MEDIAN_TRICK_HPP

#include <bits/stdc++.h>
#include "types.hpp"
using namespace std;

t_double median_trick(function<t_double()> fn, const int p, const t_double delta) {
    const t_double coef = 2.0 / (log2(p) - log2(p - 1));
    const int t = ceil(coef * log2(1.0 / delta));

    cout << "median tries: " << t << " -> ";
    cout.flush();
    long tot_time = 0.0;
    vector<t_double> results;
    for (int i = 0; i < t; i++) {
        auto start_t = chrono::high_resolution_clock::now();
        t_double r = fn();
        auto time_taken = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start_t);
        tot_time += time_taken.count();
        results.push_back(r);
        cout << ".";
        cout.flush();
    }
    tot_time /= t;
    cout << fixed << setprecision(3) << " average time per median trial: " << tot_time << " ms" << endl;
    sort(begin(results), end(results));
    return results[t / 2];
}

#endif //NETWORK_RELIABILITY_UNRELIABILITY_MEDIAN_TRICK_HPP