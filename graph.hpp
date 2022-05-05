

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
    unordered_map<t_node, vector<int>> adj;

    // probability an edge fails
    t_double p;

    // edges of graph
    vector<t_edge> edges;

    // number of edges in graph
    int m;

    Graph() {
        p = 0.5;
        m = 0;
    }

    Graph(vector<t_node> &nodes, vector<t_edge> &edges, t_double p) {
        this->nodes = nodes;
        this->edges = edges;
        this->p = p;
        adj.clear();
        m = 0;
        for (int i = 0; i < edges.size(); i++) {
            auto e = edges[i];
            adj[e.from].push_back(i);
            adj[e.to].push_back(i);
            m += e.cnt;
        }
    }

    void add_edge(t_edge e) {
        edges.push_back(e);
        int id = (int)edges.size() - 1;
        adj[e.from].push_back(id);
        adj[e.to].push_back(id);
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
        return Graph(nodes, new_edges, p);
    }

    Graph subgraph(vector<t_edge> new_edges) {
        return Graph(nodes, new_edges, p);
    }

    bool connected() {
        unordered_set<t_node> visited;
        function<void(int)> DFS = [&](t_node node) {
            if (visited.find(node) != visited.end()) return;
            visited.insert(node);

            for (auto id: adj[node])
                DFS(edges[id].other(node));
        };
        DFS(nodes[0]);
        return visited.size() == get_n();
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

    return {Graph(new_nodes, new_edges, g.p), removed_edges};
}

Graph to_unit_edges(Graph &g) {
    vector<t_edge> new_edges;
    new_edges.reserve(g.get_m());
    for(auto e: g.edges) {
        int cnt = e.cnt;
        e.cnt = 1;
        for(int i = 0; i < cnt; i++)
            new_edges.push_back(e);
    }

    return Graph(g.nodes, new_edges, g.p);
}

void print_graph(Graph &g, ofstream &out, bool complete = false) {
    out << g.get_n() << " " << g.edges.size() << endl;
    for(auto &e: g.edges) {
        if(complete)
            out << e.from << " " << e.to << " " << e.id << " " << e.cnt << endl;
        else
            out << e.from << " " << e.to << endl;
    }
}

#endif //NETWORK_RELIABILITY_UNRELIABILITY_GRAPH_HPP
