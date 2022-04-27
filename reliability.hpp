// Reliability - O( eps^(-2) (1-pmax)^(-3) m^2 n^3 )

#ifndef NETWORK_RELIABILITY_UNRELIABILITY_RELIABILITY_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_RELIABILITY_HPP

#include <bits/stdc++.h>
#include "graph.hpp"
#include "random.hpp"
#include "types.hpp"

using namespace std;
// TODO: repair
DiGraph sample_root_connected(DiGraph g) {
    int n = g.get_n();

    vector<t_edge> s_edges;
    for (auto e: g.edges)
        if (!Random::get_ber(g.p))
            s_edges.push_back(e);

    DiGraph gs = g.subgraph(s_edges);

    int steps = 0;
    while (true) {
        auto visited = gs.get_root_connected_nodes();
        if (visited.size() == n) break;

        vector<t_edge> new_edges;
        copy_if(begin(gs.edges), end(gs.edges), back_inserter(new_edges),
                [&](t_edge e) { return visited.find(e.from) != visited.end(); });

        vector<t_node> disconnected_nodes;
        copy_if(begin(gs.nodes), end(gs.nodes), back_inserter(disconnected_nodes),
                [&](t_node x) { return visited.find(x) == visited.end(); });

        for (auto node: disconnected_nodes)
            for (auto e: g.adj[node])
                if (!Random::get_ber(g.p))
                    new_edges.push_back(e);

        gs = g.subgraph(new_edges);
        steps++;
    }

    return gs;
}

t_double compute_r(DiGraph g, DiGraph gi, vector<t_edge> cntr_edges, const t_double eps) {
    // defined in proposition 9
    int s = ceil(5.0 * 1.0 / ((1.0 - g.p) * (1.0 - g.p) * eps * eps) * (g.get_n() - 1.0));

    int cnt = 0;
    for (int i = 0; i < s; i++) {
        //  sample gs connected in g contracted (gi)
        DiGraph gs = sample_root_connected(gi);

        vector<int> ids(gs.edges.size());
        transform(begin(gs.edges), end(gs.edges), begin(ids), [](t_edge e) { return e.id; });

        DiGraph gs_big = g.subgraph(ids);
        // Add edges independently
        for (auto e: cntr_edges)
            if (!Random::get_ber(g.p))
                gs_big.add_edge(e);

        // check if connected
        if (gs_big.is_root_connected())
            cnt++;
    }

    t_double r = t_double(cnt) / t_double(s);
    return r;
}

t_double compute_reliability(DiGraph g, const t_double eps) {
    t_double zreach = 1.0;
    int n = g.get_n();
    // we consider root = g.nodes[0] and always contract the last two nodes in g.
    for (int i = 1; i < n; i++) {
        auto[gi, cntr_edges] = unite(g, g.nodes[n - i - 1], g.nodes[n - i]);
        t_double r = compute_r(g, gi, cntr_edges, eps);
        zreach *= r;
        g = gi;
    }
    return zreach;
}

#endif //NETWORK_RELIABILITY_UNRELIABILITY_RELIABILITY_HPP

