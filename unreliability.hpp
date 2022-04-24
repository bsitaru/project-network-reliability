// Uneliability - O( n^2.78 ) - might be slower because of mincut

#ifndef NETWORK_RELIABILITY_UNRELIABILITY_UNRELIABILITY_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_UNRELIABILITY_HPP

#include <bits/stdc++.h>
#include "brute.hpp"
#include "graph.hpp"
#include "mincut.hpp"
#include "random.hpp"
#include "utils.hpp"
#include "types.hpp"

using namespace std;

const int NODES_BRUTE = 3;

t_double monte_carlo(Graph &g) {
    vector<t_edge> new_edges;
    for (auto &e: g.edges)
        if (Random::get_ber(1.0 - g.p))
            new_edges.push_back(e);
    Graph gs = g.subgraph(new_edges);
    if (gs.connected()) return 0.0;
    return 1.0;
}

Graph contract(Graph &g, t_double q) {
    UnionFind uf(g.nodes);

    for (auto &e: g.edges)
        if (Random::get_ber(1.0 - power(q, e.cnt))) {
            uf.unite(e.from, e.to);
        }

    vector<t_node> new_nodes;
    for (auto &node: g.nodes)
        if (uf.find(node) == node)
            new_nodes.push_back(node);

    vector<t_edge> new_edges;
    map<pair<t_node, t_node>, t_edge> edges;
    for (auto &e: g.edges) {
        auto f_from = uf.find(e.from), f_to = uf.find(e.to);
        if (f_from > f_to) swap(f_from, f_to);
        if (f_from != f_to) {
            if (!edges.count({f_from, f_to}))
                edges[{f_from, f_to}] = {f_from, f_to, e.id, e.cnt};
            else
                edges[{f_from, f_to}].cnt += e.cnt;
//            new_edges.push_back({f_from, f_to, e.id});
        }
    }

    for (const auto &x: edges)
        new_edges.push_back(x.second);

    return Graph(new_nodes, new_edges, g.p, false);
}

t_double unreliability(Graph &g, int depth = 1) {
    if (g.get_n() <= NODES_BRUTE)
        return brute_unreliability(g);

    int c;
    debug_measure_time([&]() {
//        c = KSMincut::fastmincut(g);
        c = EKMincut::mincut(g);
    }, "fastmincut", false);


    t_double pc = power(g.p, c);

    if (pc >= 0.5)
        return monte_carlo(g);

    t_double q = n_root(0.5, c);
    t_double ans = 0.0;
    for (int i = 0; i < 2; i++) {
        Graph gc = contract(g, q);
        gc.p = g.p / q;
        ans += unreliability(gc, depth + 1) / 2.0;

    }
    return ans;
}

t_double compute_unreliability(Graph &g, const t_double eps) {
    // t = 4 * r * eps^(-2)
    int t = ceil(4.0 / (eps * eps));

    t_double coef = 1.0 / t_double(t);
    t_double ans = 0.0;

    long tot_time = 0;
    for (int i = 0; i < t; i++) {
        auto time_s = chrono::high_resolution_clock::now();
        t_double ug = unreliability(g);
        auto time_taken_brute = chrono::duration_cast<chrono::milliseconds>(
                chrono::high_resolution_clock::now() - time_s);
//        cout << "sample time: " << time_taken_brute.count() << " ms" << endl;
        tot_time += time_taken_brute.count();

        ans += ug * coef;
    }

//    cout << "average sample time: " << tot_time / t << " ms" << endl;

    return ans;
}

#endif //NETWORK_RELIABILITY_UNRELIABILITY_UNRELIABILITY_HPP