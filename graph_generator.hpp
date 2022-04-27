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

        return Graph(nodes, edges, p);
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

        return Graph(nodes, edges, graph_p);
    }

    Graph complete_graph(int n) {
        return Generator::erdos_renyi(n, 1.0);
    }

    Graph random(int n, int m) {
        typedef pair<int, int> pii;

        vector<t_node> nodes;
        for (int i = 0; i < n; i++) nodes.push_back(i);

        assert(m >= n - 1);
        vector<pii> all;
        for(int i = 0; i < n; i++)
            for(int j = i + 1; j < n; j++)
                all.emplace_back(i, j);
        shuffle(begin(all), end(all), Random::gen);

        UnionFind uf(nodes);

        int id = 0;
        vector<t_edge> edges;
        for(int i = 0; i < all.size(); i++) {
            auto [x, y] = all[i];
            if(uf.find(x) == uf.find(y))    continue;
            uf.unite(x, y);
            edges.push_back({x, y, id++, 1});
            swap(all[i], all.back());
            all.pop_back();
            i--;
        }

        m -= n - 1;
        for(int i = 0; i < m; i++) {
            auto [x, y] = all[i];
            edges.push_back({x, y, id++, 1});
        }

        t_double graph_p = 0.01;

        return Graph(nodes, edges, graph_p);
    }

    Graph dodecahedron() {

        /*

         p = 0.1 -> u = 0.02286916405937333696

         */

        vector<int> nodes;
        for(int i = 1; i <= 20; i++)    nodes.push_back(i);

        vector< pair<int, int> > e = {
                {1, 2}, {1, 3}, {1, 4},
                {2, 5}, {2, 6},
                {3, 7}, {3, 8},
                {4, 9}, {4, 10},
                {5, 10}, {5, 11},
                {6, 7}, {6, 12},
                {7, 13},
                {8, 9}, {8, 14},
                {9, 15},
                {10, 16},
                {11, 12}, {11, 17},
                {12, 18},
                {13, 14}, {13, 18},
                {14, 19},
                {15, 16}, {15, 19},
                {16, 17},
                {17, 20},
                {18, 20},
                {19, 20}
        };

        vector<t_edge> edges;
        int id = 0;
        for(auto ee: e) {
            edges.push_back({ee.first, ee.second, id++, 1});
        }

        t_double graph_p = 0.01;

        return Graph(nodes, edges, graph_p);
    }

    Graph k4(int k) {
        int n = 2 * k;
        vector<t_node> nodes;
        for(int i = 0; i < n; i++)  nodes.push_back(i);

        vector<t_edge> edges;
        int id = 0;
        for(int i = 0; i < k; i++) {
            edges.push_back({i, i + k, id++, 1});

            if(i < k - 1) {
                edges.push_back({i, i + 1, id++, 1});
                edges.push_back({i, i + k + 1, id++, 1});
                edges.push_back({i + k, i + k + 1, id++, 1});
                edges.push_back({i + k, i + 1, id++, 1});
            }
        }

        t_double graph_p = 0.01;

        return Graph(nodes, edges, graph_p);
    }
}

#endif //NETWORK_RELIABILITY_UNRELIABILITY_GRAPH_GENERATOR_HPP
