#ifndef NETWORK_RELIABILITY_UNRELIABILITY_GRAPH_GENERATOR_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_GRAPH_GENERATOR_HPP

#include <bits/stdc++.h>

#include "random.hpp"
#include "types.hpp"
#include "graph.hpp"

using namespace std;
namespace Generator {
    Graph example() {
//        nodes = {0, 1, 2, 3};
//
//        edges = {{0, 1, 0},
//                 {1, 2, 1},
//                 {2, 3, 2}};

        vector<t_node> nodes = {0, 1, 2, 3, 4, 5};
        vector<t_edge> edges = {
                {0, 1, 0, 1},
                {0, 2, 1, 1},
                {1, 3, 2, 1},
                {1, 2, 3, 1},
                {2, 3, 4, 1},
                {2, 4, 5, 1},
                {3, 4, 6, 1},
                {3, 5, 7, 1},
                {4, 5, 8, 1}
        };
        t_double p = 0.01;
        bool directed = false;

        return Graph(nodes, edges, p, directed);
    }

    Graph erdos_renyi(int n, t_double p) {
        vector<t_node> nodes;
        for (int i = 1; i <= n; i++) nodes.push_back(i);

        int id = 0;
        vector<t_edge> edges;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (Random::get_ber(p)) {
                    edges.push_back({nodes[i], nodes[j], id++, 1});
                }

        t_double graph_p = 0.01;
        bool directed = false;

        return Graph(nodes, edges, graph_p, directed);
    }

    Graph complete_graph(int n) {
        return Generator::erdos_renyi(n, 1.0);
    }
}

#endif //NETWORK_RELIABILITY_UNRELIABILITY_GRAPH_GENERATOR_HPP
