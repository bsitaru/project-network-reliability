

#ifndef NETWORK_RELIABILITY_UNRELIABILITY_DIGRAPH_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_DIGRAPH_HPP

#include <bits/stdc++.h>

#include "graph.hpp"
#include "random.hpp"
#include "types.hpp"
#include "utils.hpp"

using namespace std;
// TODO: repair
struct DiGraph {

    // nodes of graph
    vector<t_node> nodes;

    // adjacent nodes for each node
    unordered_map<t_node, vector<t_edge>> adj, radj;

    // probability an edge fails
    t_double p;

    // edges of graph
    vector<t_edge> edges;

    // number of edges in graph
    int m;

    DiGraph() {
        p = 0.5;
        m = 0;
    }

    DiGraph(Graph g) {
        ;
    }

    DiGraph(vector<t_node> &nodes, vector<t_edge> &edges, t_double p) {
        this->nodes = nodes;
        this->edges = edges;
        this->p = p;
        adj.clear();
        radj.clear();
        m = 0;
        for (auto &e: edges) {
            adj[e.from].push_back(e);
            m += e.cnt;
            radj[e.to].push_back(e);
        }
    }

    void add_edge(t_edge e) {
        edges.push_back(e);
        adj[e.from].push_back(e);
        radj[e.to].push_back(e);
        m += e.cnt;
    }

    int get_n() const {
        return (int) nodes.size();
    }

    int get_m() const {
        return m;
    }

    DiGraph subgraph(vector<int> ids) {
        set<int> ids_set(begin(ids), end(ids));
        vector<t_edge> new_edges;
        for (auto e: edges)
            if (ids_set.find(e.id) != ids_set.end())
                new_edges.push_back(e);
        return DiGraph(nodes, new_edges, p);
    }

    DiGraph subgraph(vector<t_edge> new_edges) {
        return DiGraph(nodes, new_edges, p);
    }

    DiGraph to_directed() {
        assert(!directed);
        vector<t_edge> new_edges = edges;
        for (auto e: edges)
            new_edges.push_back({e.to, e.from, e.id + (int) edges.size()});
        return DiGraph(nodes, new_edges, p);
    }

    bool connected() {
        assert(!directed);
        unordered_set<t_node> visited;
        function<void(int)> DFS = [&](t_node node) {
            if (visited.find(node) != visited.end()) return;
            visited.insert(node);

            for (auto e: adj[node])
                DFS(e.to);
        };
        DFS(nodes[0]);
        return visited.size() == get_n();
    }

    unordered_set<t_node> get_root_connected_nodes() {
        assert(directed);
        unordered_set<t_node> visited;
        function<void(int)> DFS = [&](t_node node) {
            if (visited.find(node) != visited.end()) return;
            visited.insert(node);

            for (auto e: radj[node])
                DFS(e.from);
        };
        DFS(nodes[0]);
        return visited;
    }

    bool is_root_connected() {
        auto visited = get_root_connected_nodes();
        if (visited.size() != get_n()) return false;
        return true;
    }
};

static pair<DiGraph, vector<t_edge>> unite(DiGraph &g, t_node x, t_node y) {
    if (x == y) return {};
    vector<t_edge> new_edges;
    vector<t_edge> removed_edges;
    for (auto[from, to, id, cnt]: g.edges) {
        if ((from == x && to == y) || (from == y && to == x)) {
            removed_edges.push_back({from, to, id, cnt});
            continue;
        };

        if (from == y) from = x;
        if (to == y) to = x;
        new_edges.push_back({from, to, id, cnt});
    }
    auto new_nodes = g.nodes;
    auto it = find(begin(new_nodes), end(new_nodes), y);
    assert(it != end(new_nodes));
    new_nodes.erase(it);

    return {DiGraph(new_nodes, new_edges, g.p), removed_edges};
}
//
//Graph to_unit_edges(Graph &g) {
//    vector<t_edge> new_edges;
//    new_edges.reserve(g.get_m());
//    for(auto e: g.edges) {
//        int cnt = e.cnt;
//        e.cnt = 1;
//        for(int i = 0; i < cnt; i++)
//            new_edges.push_back(e);
//    }
//
//    return Graph(g.nodes, new_edges, g.p, g.directed);
//}

#endif //NETWORK_RELIABILITY_UNRELIABILITY_DIGRAPH_HPP
