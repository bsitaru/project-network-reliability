#ifndef NETWORK_RELIABILITY_UNRELIABILITY_BRUTE_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_BRUTE_HPP

#include <bits/stdc++.h>
#include "graph.hpp"
#include "types.hpp"

using namespace std;

t_double brute(Graph g, bool reliability) {
    t_double ans_rel = 0.0, ans_unrel = 0.0;
    vector<int> ids;
    function<void(int, t_double)> bck = [&](int i, t_double w) {
        if (i >= g.edges.size()) {
            Graph gs = g.subgraph(ids);
            if (gs.connected())
                ans_rel += w;
            else
                ans_unrel += w;
            return;
        }

        t_double p = power(g.p, g.edges[i].cnt);

        bck(i + 1, w * p);

        ids.push_back(g.edges[i].id);
        bck(i + 1, w * (1.0 - p));
        ids.pop_back();
    };

    bck(0, 1.0);

    return reliability ? ans_rel : ans_unrel;
}

t_double brute_reliability(Graph g) {
    return brute(g, true);
}

t_double brute_unreliability(Graph g) {
    return brute(g, false);
}


#endif //NETWORK_RELIABILITY_UNRELIABILITY_BRUTE_HPP
