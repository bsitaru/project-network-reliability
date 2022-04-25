#ifndef NETWORK_RELIABILITY_UNRELIABILITY_UTILS_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_UTILS_HPP

#include "types.hpp"

t_double power(t_double x, int e, t_double ans = 1.0) {
    if(e == 0)  return ans;
    return power(x * x, e >> 1, ans * ( (e & 1) ? x : 1.0 ));
}

t_double n_root(t_double x, int n) {
    t_double l = 0.0, r = 1.0;
    int steps = 60;
    for(int i = 0; i < steps; i++) {
        t_double m = (l + r) / 2.0;
        if(power(m, n) <= x)    l = m;
        else    r = m;
    }
    return l;
}

struct UnionFind {
    unordered_map<t_node, t_node> uf;

    UnionFind(const vector<t_node> &nodes) {
        for (auto node: nodes) uf[node] = node;
    }

    t_node find(t_node x) {
        if (uf[x] == x) return x;
        return uf[x] = find(uf[x]);
    }

    void unite(t_node x, t_node y) {
        t_node fx = find(x), fy = find(y);
        if(fx != fy)
            uf[fy] = fx;
    }
};

void debug_measure_time(function<void()> fn, string text, bool print = true) {
    auto t_start = chrono::high_resolution_clock::now();
    fn();
    auto time_taken = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - t_start);
    if(print)   cout << text << ": " << time_taken.count() << " us" << endl;
}

void debug_measure_time_ms(function<void()> fn, string text, bool print = true) {
    auto t_start = chrono::high_resolution_clock::now();
    fn();
    auto time_taken = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t_start);
    if(print)   cout << text << ": " << time_taken.count() << " ms" << endl;
}

#endif //NETWORK_RELIABILITY_UNRELIABILITY_UTILS_HPP
