

#ifndef NETWORK_RELIABILITY_UNRELIABILITY_GRAPH_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_GRAPH_HPP

#include <bits/stdc++.h>

#include "random.hpp"
#include "types.hpp"
#include "utils.hpp"

using namespace std;

struct Graph {

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

    // directed
    bool directed;

    // mincut
    int c;

    Graph() {
        p = 0.5;
        directed = false;
        c = -1;
        m = 0;
    }

    Graph(vector<t_node> &nodes, vector<t_edge> &edges, t_double p, bool directed) {
        this->nodes = nodes;
        this->edges = edges;
        this->p = p;
        this->directed = directed;
        adj.clear();
        radj.clear();
        m = 0;
        for (auto e: edges) {
            adj[e.from].push_back(e);
            m += e.cnt;
            if (directed)
                radj[e.to].push_back(e);
            else
                adj[e.to].push_back({e.to, e.from, e.id, e.cnt});
        }
        this->c = -1;
    }

    void add_edge(t_edge e) {
        edges.push_back(e);
        adj[e.from].push_back(e);
        if (directed)
            radj[e.to].push_back(e);
        else
            adj[e.to].push_back({e.to, e.from, e.id, e.cnt});
        m += e.cnt;
    }

    int get_n() const {
        return (int) nodes.size();
    }

    int get_m() const {
        return m;
    }

    Graph subgraph(vector<int> ids) {
        set<int> ids_set(begin(ids), end(ids));
        vector<t_edge> new_edges;
        for (auto e: edges)
            if (ids_set.find(e.id) != ids_set.end())
                new_edges.push_back(e);
        return Graph(nodes, new_edges, p, directed);
    }

    Graph subgraph(vector<t_edge> new_edges) {
        return Graph(nodes, new_edges, p, directed);
    }

    Graph to_directed() {
        assert(!directed);
        vector<t_edge> new_edges = edges;
        for (auto e: edges)
            new_edges.push_back({e.to, e.from, e.id + (int) edges.size()});
        return Graph(nodes, new_edges, p, true);
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

static pair<Graph, vector<t_edge>> unite(Graph &g, t_node x, t_node y) {
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

    return {Graph(new_nodes, new_edges, g.p, g.directed), removed_edges};
}

#endif //NETWORK_RELIABILITY_UNRELIABILITY_GRAPH_HPP
