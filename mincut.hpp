#ifndef NETWORK_RELIABILITY_UNRELIABILITY_MINCUT_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_MINCUT_HPP

#include <bits/stdc++.h>

#include "graph.hpp"
#include "random.hpp"
#include "types.hpp"
#include "utils.hpp"

#include "mincut_linear.hpp"

using namespace std;

struct SegmentTree {
    int n;
    vector<int> aint;

    SegmentTree(int n) {
        this->n = n;
        aint = vector<int>(4 * n, 0);
    }

    SegmentTree(int n, vector<int> &v) {
        this->n = n;
        aint = vector<int>(4 * n, 0);
        build(1, 0, n - 1, v);
    }

    void build(int nod, int st, int dr, vector<int> &v) {
        if (st == dr) {
            aint[nod] = v[st];
            return;
        }

        int m = (st + dr) / 2;
        build(nod * 2, st, m, v);
        build(nod * 2 + 1, m + 1, dr, v);
        aint[nod] = aint[nod * 2] + aint[nod * 2 + 1];
    }

    int sample() {
        int sum = aint[1];
        int x = Random::get_int(1, sum);

        int nod = 1, st = 0, dr = n - 1;
        while(st < dr) {
            int m = (st + dr) / 2;
            if(aint[nod * 2] >= x) {
                nod = nod * 2;
                dr = m;
            } else {
                sum -= aint[nod * 2];
                nod = nod * 2 + 1;
                st = m + 1;
            }
        }
        aint[nod] = 0;
        nod /= 2;
        while(nod > 0) {
            aint[nod] = aint[nod * 2] + aint[nod * 2 + 1];
            nod /= 2;
        }
        return st;
    }
};

namespace KSMincut {    // karger-stein

    constexpr static const t_double sqrt2 = 1.414213562373095;

    void weighted_shuffle(vector<t_edge> &v) {
        int n = (int) v.size();
        vector<int> weights(n);
        transform(begin(v), end(v), begin(weights), [](t_edge e) { return e.cnt; });

        SegmentTree aint(n, weights);

        vector<t_edge> shuffled(n);
        for (int i = 0; i < n; i++) {
            int id = aint.sample();
            shuffled[i] = v[id];
        }
        v = shuffled;
    }

    static Graph contract(Graph g, int t) {
        int n = g.get_n();
        auto edges = g.edges;
        UnionFind uf(g.nodes);
        shuffle(begin(edges), end(edges), Random::gen);
//        weighted_shuffle(edges);

        for (auto &e: edges) {
            if (uf.find(e.from) == uf.find(e.to)) continue;
            uf.unite(e.from, e.to);
            n--;
            if (n <= t) break;
        }

        vector<t_node> new_nodes;
        for (auto &node: g.nodes)
            if (uf.find(node) == node)
                new_nodes.push_back(node);

        vector<t_edge> new_edges;
        map<pair<t_node, t_node>, t_edge> edges_cnt;
        for (auto &e: g.edges) {
            auto f_from = uf.find(e.from), f_to = uf.find(e.to);
            if (f_from > f_to) swap(f_from, f_to);
            if (f_from != f_to) {
                if (!edges_cnt.count({f_from, f_to}))
                    edges_cnt[{f_from, f_to}] = {f_from, f_to, e.id, e.cnt};
                else
                    edges_cnt[{f_from, f_to}].cnt += e.cnt;
//            new_edges.push_back({f_from, f_to, e.id});
            }
        }

        for (const auto &x: edges_cnt)
            new_edges.push_back(x.second);

        return Graph(new_nodes, new_edges, g.p);
    }

    static int mincut(Graph &g) {
        int n = g.get_n();
        int t = ceil(t_double(n) * t_double(n - 1) * log(t_double(n)) / 2.0);
        int ans = g.get_m();
        for (int i = 0; i < t; i++) {
            Graph gc = contract(g, 2);
            ans = min(ans, gc.get_m());
        }
        return ans;
    }

    static int _fastmincut(Graph &g, int depth = 1) {
        int n = g.get_n();
        if (n <= 6) {
            return mincut(g);
        }
        int t = ceil(1.0 + t_double(n) / t_double(sqrt2));
        Graph g1 = contract(g, t);
        Graph g2 = contract(g, t);


        int mc1 = _fastmincut(g1, depth + 1);
        int mc2 = _fastmincut(g2, depth + 1);
        return min(mc1, mc2);
    }

    static int fastmincut(Graph &g) {
        Graph my_g = to_unit_edges(g);
        int n = my_g.get_n();
        int t = ceil(t_double(n) * log(n) / t_double(n - 1));
        int ans = my_g.get_m();
        for (int i = 0; i < t; i++) {
            int cut = _fastmincut(my_g);
            ans = min(ans, cut);
        }
        return ans;
    }
};

namespace EKFastMincut {    // edmond-karp
    typedef pair<int, int> pii;
    typedef pair<pii, int> p3i;

    const int NMAX = 35;

    class EdmondsKarp {
    public:

        int N, S, T;
        vector<p3i> edges;
        int cp[NMAX][NMAX];
        int flw[NMAX][NMAX];
        vector<int> edg[NMAX];
        int d[NMAX], f[NMAX];
        bool initGraph;

        EdmondsKarp() {
            N = S = T = 0;
            memset(cp, 0, sizeof(cp));
            memset(flw, 0, sizeof(flw));
            memset(d, 0, sizeof(d));
            memset(f, 0, sizeof(f));
            initGraph = false;
        }

        void clear() {
            edges.clear();
            for (int i = 0; i <= N; i++) edg[i].clear();
            N = S = T = 0;
            memset(cp, 0, sizeof(cp));
            memset(flw, 0, sizeof(flw));
            memset(d, 0, sizeof(d));
            memset(f, 0, sizeof(f));
            initGraph = false;
        }

        void addEdge(int x, int y, int c) {
            edges.push_back({{x, y}, c});
            initGraph = false;
        }

        void createGraph() {
            memset(cp, 0, sizeof(cp));

            N = 0;
            for (auto e: edges) {
                N = max(N, max(e.first.first, e.first.second));
                cp[e.first.first][e.first.second] += e.second;
            }

            for (int i = 1; i <= N; i++) edg[i].clear();
            for (int i = 1; i <= N; i++)
                for (int j = i + 1; j <= N; j++)
                    if (cp[i][j] || cp[j][i]) {
                        edg[i].push_back(j);
                        edg[j].push_back(i);
                    }

            initGraph = true;
        }

        void BFS(int start) {
            for (int i = 1; i <= N; i++) d[i] = 0;
            d[start] = 1;
            queue<int> q;
            q.push(start);
            while (!q.empty()) {
                int nod = q.front();
                q.pop();
                for (auto nxt: edg[nod])
                    if (!d[nxt] && cp[nod][nxt] - flw[nod][nxt]) {
                        d[nxt] = 1 + d[nod];
                        f[nxt] = nod;
                        q.push(nxt);
                    }
            }
        }

        int getFlow(int nod) {
            if (nod == S) return (1 << 30);
            int ans = min(cp[f[nod]][nod] - flw[f[nod]][nod], getFlow(f[nod]));
            return ans;
        }

        void pushFlow(int nod, int nxt, int flow) {
            flw[nod][nxt] += flow;
            flw[nxt][nod] -= flow;
            if (nod == S) return;
            pushFlow(f[nod], nod, flow);
        }

        int getMaxFlow(int _S, int _T) {
            S = _S;
            T = _T;

            if (!initGraph) createGraph();

            memset(flw, 0, sizeof(flw));
            int flow = 0;
            bool newFlow = true;
            while (newFlow) {
                newFlow = false;
                BFS(S);

                for (auto nod: edg[T])
                    if (cp[nod][T] - flw[nod][T] && d[nod] != 0) {
                        int addFlow = getFlow(nod);
                        addFlow = min(addFlow, cp[nod][T] - flw[nod][T]);
                        if (addFlow) {
                            pushFlow(nod, T, addFlow);
                            flow += addFlow;
                            newFlow = true;
                        }
                    }
            }

            return flow;
        }
    };

    static EdmondsKarp flow;

    int mincut(Graph &g) {
        flow.clear();
        for (auto e: g.edges) {
            flow.addEdge(e.from, e.to, e.cnt);
            flow.addEdge(e.to, e.from, e.cnt);
        }

        int s = g.nodes[0];
        int ans = g.get_m();
        for (auto t: g.nodes) {
            if (s == t) continue;
            int c = flow.getMaxFlow(s, t);
            ans = min(ans, c);
        }
        return ans;
    }
}

namespace EKMincut {    // edmond-karp
    struct FlowEdge {
        int from, to;
        int cap, flow = 0;
        FlowEdge(int from, int to, int cap) : from(from), to(to), cap(cap) {}
    };

    struct EdmondsKarp {
        int N, M, S, T;
        vector<FlowEdge> edges;
        vector< vector<int> > adj;
        vector<int> d, f;

        EdmondsKarp() {
            N = M = S = T = 0;
        }

        void init(int N, int exp_m) {
            this->N = N;
            M = 0;
            edges.clear();
            edges.reserve(exp_m);
            adj.resize(N);
            d.resize(N); f.resize(N);
            for(int i = 0; i < N; i++)  {
                adj[i].clear();
                d[i] = f[i] = 0;
            }
        }

        void add_edge(int x, int y, int c) {
            edges.emplace_back(x, y, c);
            edges.emplace_back(y, x, 0);
            adj[x].push_back(M);
            adj[y].push_back(M + 1);
            M += 2;
        }

        void BFS(int start) {
            for (int i = 0; i < N; i++) d[i] = 0;
            d[start] = 1;
            queue<int> q;
            q.push(start);
            while (!q.empty()) {
                int nod = q.front();
                q.pop();
                for (auto nxt_id: adj[nod]) {
                    int nxt = edges[nxt_id].to;
                    if (!d[nxt] && edges[nxt_id].cap - edges[nxt_id].flow) {
                        d[nxt] = 1 + d[nod];
                        f[nxt] = nxt_id;
                        q.push(nxt);
                    }
                }
            }
        }

        inline int get_flow(int nod) {
            int ans = 1 << 30;
            while(nod != S) {
                int id = f[nod];
                ans = min(ans, edges[id].cap - edges[id].flow);
                nod = edges[id].from;
            }
            return ans;
        }

        inline void push_flow(int nod, int id, int flow) {
            while(true) {
                edges[id].flow += flow;
                edges[id ^ 1].flow -= flow;
                if(nod == S)    break;
                id = f[nod];
                nod = edges[id].from;
            }
        }

        int get_max_flow(int _S, int _T) {
            S = _S;
            T = _T;

            for(auto &e: edges)  e.flow = 0;

            int flow = 0;
            bool newFlow = true;
            while (newFlow) {
                newFlow = false;
                BFS(S);

                for (auto id: adj[T]) {
                    int nod = edges[id].to;
                    id ^= 1;
                    if (nod != T && edges[id].cap - edges[id].flow && d[nod] != 0) {
                        int addFlow = get_flow(nod);
                        addFlow = min(addFlow, edges[id].cap - edges[id].flow);
                        if (addFlow) {
                            push_flow(nod, id, addFlow);
                            flow += addFlow;
                            newFlow = true;
                        }
                    }
                }
            }

            return flow;
        }

        void unite(int x, int y) {
            for(auto id: adj[y]) {
                if(edges[id].to != x)   adj[x].push_back(id);
                edges[id].from = x;
                edges[id ^ 1].to = x;
            }
        }
    };

    static EdmondsKarp flow;

    int mincut(Graph &g) {
        flow.init(g.get_n(), g.edges.size());
        unordered_map<t_node, int> mp;
        int n = 0;
        for (auto e: g.edges) {
            int x = e.from, y = e.to;
            if(!mp.count(x))    mp[x] = n++;
            if(!mp.count(y))    mp[y] = n++;
            x = mp[x]; y = mp[y];
            flow.add_edge(x, y, e.cnt);
            flow.add_edge(y, x, e.cnt);
        }

        int s = 0;
        int ans = g.get_m();
        for (int t = 1; t < n; t++) {
            int c = flow.get_max_flow(s, t);
            ans = min(ans, c);
        }
//        for (int i = 0; i < n - 1; i++) {
//            int t;
//            for(auto id: flow.adj[s])
//                if(flow.edges[id].to != s) {
//                    t = flow.edges[id].to;
//                    break;
//                }
//            int c = flow.get_max_flow(s, t);
//            ans = min(ans, c);
//            flow.unite(s, t);
//        }
        return ans;
    }
}

namespace DinMincut {   // dinic
    struct FlowEdge {
        int v, u;
        long long cap, flow = 0;
        FlowEdge(int v, int u, long long cap) : v(v), u(u), cap(cap) {}
    };

    struct Dinic {
        const long long flow_inf = 1e18;
        vector<FlowEdge> edges;
        vector<vector<int>> adj;
        int n, m = 0;
        int s, t;
        vector<int> level, ptr;
        queue<int> q;

        Dinic() {
            ;
        }

        void init(int n) {
            adj.resize(n);
            level.resize(n);
            ptr.resize(n);
            edges.clear();
            m = 0;
            for(int i = 0; i < n; i++) {
                level[i] = ptr[i] = 0;
                adj[i].clear();
            }
            while(!q.empty())   q.pop();
        }

        void add_edge(int v, int u, long long cap) {
            edges.emplace_back(v, u, cap);
            edges.emplace_back(u, v, 0);
            adj[v].push_back(m);
            adj[u].push_back(m + 1);
            m += 2;
        }

        bool bfs() {
            while (!q.empty()) {
                int v = q.front();
                q.pop();
                for (int id : adj[v]) {
                    if (edges[id].cap - edges[id].flow < 1)
                        continue;
                    if (level[edges[id].u] != -1)
                        continue;
                    level[edges[id].u] = level[v] + 1;
                    q.push(edges[id].u);
                }
            }
            return level[t] != -1;
        }

        long long dfs(int v, long long pushed) {
            if (pushed == 0)
                return 0;
            if (v == t)
                return pushed;
            for (int& cid = ptr[v]; cid < (int)adj[v].size(); cid++) {
                int id = adj[v][cid];
                int u = edges[id].u;
                if (level[v] + 1 != level[u] || edges[id].cap - edges[id].flow < 1)
                    continue;
                long long tr = dfs(u, min(pushed, edges[id].cap - edges[id].flow));
                if (tr == 0)
                    continue;
                edges[id].flow += tr;
                edges[id ^ 1].flow -= tr;
                return tr;
            }
            return 0;
        }

        long long flow(int s, int t){
            this->s = s; this->t = t;
            long long f = 0;
            for(auto &e: edges)
                e.flow = 0;
            while (true) {
                fill(level.begin(), level.end(), -1);
                level[s] = 0;
                q.push(s);
                if (!bfs())
                    break;
                fill(ptr.begin(), ptr.end(), 0);
                while (long long pushed = dfs(s, flow_inf)) {
                    f += pushed;
                }
            }
            return f;
        }
    };

    static Dinic flow;

    int mincut(Graph &g) {
        flow.init(g.get_n());
        unordered_map<t_node, int> mp;
        int n = 0;
        for (auto e: g.edges) {
            int x = e.from, y = e.to;
            if(!mp.count(x))    mp[x] = n++;
            if(!mp.count(y))    mp[y] = n++;
            x = mp[x]; y = mp[y];
            flow.add_edge(x, y, e.cnt);
            flow.add_edge(y, x, e.cnt);
        }

        int s = 0;
        int ans = g.get_m();
        for (int t = 1; t < n; t++) {
            int c = flow.flow(s, t);
            ans = min(ans, c);
        }
        return ans;
    }
};

namespace SWMincut {    // stoer-wagner
    int mincut_phase(Graph &ig) {
        Graph g = ig;

        int n = g.get_n();
        priority_queue< pair<int,int> > pq;
        unordered_map<t_node, bool> in_set;
        int node = g.nodes[0];
        in_set[node] = true;
        unordered_map<t_node, int> connect;
        for(auto id: g.adj[node])
            connect[ g.edges[id].other(node) ] += g.edges[id].cnt;
        for(auto oth: g.nodes) {
            if(node == oth) continue;
            pq.push({connect[oth], oth});
        }
        vector<int> order;
        order.reserve(n);
        order.push_back(node);
        for(int i = 1; i < n; i++) {
            auto [val, node] = pq.top(); pq.pop();
            while(connect[node] != val || in_set[node]) {
                tie(val, node) = pq.top();
                pq.pop();
            }
            order.push_back(node);
            in_set[node] = true;
            for(auto id: g.adj[node]) {
                auto e = g.edges[id];
                int to = e.other(node);
                if(!in_set[to]) {
                    connect[to] += e.cnt;
                    pq.push({connect[to], to});
                }
            }
        }

        node = order.back();
        int c = connect[node];
        order.pop_back();
        Graph g_red = unite(g, order.back(), node).first;
        ig = g_red;
        return c;
    }

    int mincut(Graph g) {
        int ans = g.get_m();

        while(g.get_n() > 1) {
            int c = mincut_phase(g);
            ans = min(ans, c);
        }
        return ans;
    }
}

namespace BenchmarkMincut {
    int mincut(Graph &g) {
        int ans;

        profiler.start("ks-mincut");
        ans = KSMincut::mincut(g);
        profiler.stop("ks-mincut");

//        profiler.start("kar-mincut");
//        ans = KargerLinearMincut::mincut(g);
//        profiler.stop("kar-mincut");

        profiler.start("din-mincut");
        ans = DinMincut::mincut(g);
        profiler.stop("din-mincut");

        profiler.start("ek-mincut");
        ans = EKMincut::mincut(g);
        profiler.stop("ek-mincut");

        return ans;
    }
}

#endif //NETWORK_RELIABILITY_UNRELIABILITY_MINCUT_HPP
