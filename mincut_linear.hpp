#include <bits/stdc++.h>

#include "graph.hpp"
#include "random.hpp"
#include "types.hpp"
#include "utils.hpp"

namespace KargerLinearMincut_GH {
    class edge {
    public:
        int u, v, idx;
        double w;

        edge(int _u, int _v, int _idx, double _w) { u = _u, v = _v, idx = _idx, w = _w; }
    };

    typedef edge *pedge;

    class graph {
    public:
        int n, m;
        vector<vector<pedge>> adj;
        vector<pedge> E;

        graph(int _n, int _m, vector<pedge> &_E) {
            n = _n, m = _m, E = _E;
            adj.clear();
            adj.resize(n);
            for (auto it: E) adj[it->u].push_back(it), adj[it->v].push_back(it);
        }
    };

    typedef graph *pgraph;

    class tree : public graph {
    public:
        tree(int _n, vector<pedge> &_E) : graph(_n, _n - 1, _E) {}
    };

    typedef tree *ptree;

    const double eps = 1 / 5.0;

    class packeredge {
    public:
        int u, v, idx;
        priority_queue<double> l;

        packeredge(int _u, int _v, int _idx, priority_queue<double> &_l) { u = _u, v = _v, idx = _idx, l = _l; }
    };

    typedef packeredge *ppackeredge;

    class packergraph {
    public:
        int n, m;
        vector<ppackeredge> E;

        packergraph(int _n, int _m, vector<ppackeredge> &_E) {
            n = _n, m = _m, E = _E;
        }
    };

    typedef packergraph *ppackergraph;

    inline int root(int x, vector<int> &P) { return (P[x] == x) ? x : (P[x] = root(P[x], P)); }

    void dsu(int x, int y, vector<int> &P, vector<int> &sz) {
        x = root(x, P), y = root(y, P);
        if (sz[x] > sz[y]) swap(x, y);
        P[x] = y;
        sz[y] += sz[x];
    }

// Algorithm 1: Given an unweighted, undirected graph, compute packing of weight at least .4c
    pair<double, vector<ptree>> packer(ppackergraph G) {
        double W = 0;
        vector<vector<ppackeredge>> packerEdges;

        vector<int> P(G->n), sz(G->n);
        vector<pair<int, int>> sorter(G->m);
        vector<ppackeredge> href(G->m);
        for (auto it: G->E) href[it->idx] = it;
        while (true) {
            for (auto it: G->E) sorter[it->idx] = {-it->l.top(), it->idx};
            sort(sorter.begin(), sorter.end());
            vector<ppackeredge> T;
            for (int i = 0; i < G->n; i++) P[i] = i, sz[i] = 1;
            for (auto idx: sorter) {
                auto it = href[idx.second];
                int a = it->u, b = it->v;
                if (root(a, P) != root(b, P)) {
                    dsu(a, b, P, sz);
                    T.emplace_back(it);
                    double newl = -it->l.top() + (eps * eps) / (3.0 * log(G->m));
                    if (newl > 1) {
                        pair<double, vector<ptree>> res;
                        res.first = W;
                        for (auto kt: packerEdges) {
                            vector<pedge> tmpedges;
                            for (auto gt: kt) {
                                pedge tmp = new edge(gt->u, gt->v, gt->idx, 0);
                                tmpedges.emplace_back(tmp);
                            }

                            ptree tmp = new tree(G->n, tmpedges);
                            res.second.emplace_back(tmp);
                        }

                        return res;
                    }

                    it->l.pop();
                    it->l.emplace(-newl);
                }
            }

            packerEdges.emplace_back(T);
            W += (eps * eps) / (3.0 * log(G->m));
        }
    }

    typedef long long int lli;

    const double inf = 1e12 + 5, eps1 = 1.0 / 100.0, eps2 = 1.0 / 6.0, eps3 = 1.0 / 5.0;
    const double f = 3 / 2.0 - ((2.0 + eps1) * (1.0 + eps2)) / ((2.0 - eps1) * (1.0 - eps3));
    const double cmpeps = 1e-11;
    default_random_engine generator(time(NULL));

    double compute_U(const pgraph G) {
        double res = inf;
        for (auto it: G->adj) {
            double tmp = 0;
            for (auto gt: it) tmp += gt->w;
            res = min(res, tmp);
        }
        return res;
    }

// inverse transform sampling in O(ceil) time
    lli binom(lli trials, long double p, lli ceil) {
        if (p > 1 - cmpeps)
            return trials;

        uniform_real_distribution<double> distribution(0.0, 1.0);
        double u = distribution(generator);

        long double prob = pow(1 - p, trials);
        long double cum_prob = prob;
        for (int i = 0; i <= ceil; ++i) {
            if (cum_prob >= u - cmpeps) {
                return i;
            }
            prob *= (long double) (trials - i) / (i + 1) * p / (1 - p);
            cum_prob += prob;
        }

        return ceil;
    }

    vector<ptree> sample(const vector<ptree> &packing, double d, int n) {
        uniform_int_distribution<int> distribution(0, int(packing.size()) - 1);
        vector<ptree> res;
        vector<bool> marker(int(packing.size()));
        int req = ceil(-d * log(n) / log(1 - f));
        for (int i = 0; i < req; i++) {
            int tmp = distribution(generator);
            marker[tmp] = 1;
            res.emplace_back(packing[tmp]);
        }

        for (int i = 0; i < int(packing.size()); i++) {
            if (!marker[i]) {
                for (auto it: packing[i]->E) delete (it);
                delete (packing[i]);
            }
        }

        return res;
    }

// Algorithm 2: Obtain spanning trees with probability of success 1 − 1/(G->n)^d
    vector<ptree> tworespectingtrees(double d, pgraph _G) {
        vector<ptree> res;
        vector<pedge> Edash;
        double normaliser = inf;
        for (auto it: _G->E) {
            if (it->w >= cmpeps) normaliser = min(normaliser, it->w);
        }

        normaliser = 1.0 / normaliser;
        for (auto it: _G->E) {
            auto tmp = new edge(it->u, it->v, it->idx, it->w * normaliser);
            Edash.emplace_back(tmp);
        }

        pgraph Gdash = new graph(_G->n, _G->m, Edash);

        double b = 3.0 * (d + 2.0) * log(Gdash->n) / (eps2 * eps2);

        for (auto it: Gdash->E) {
            it->w = round(it->w * (1 / eps1));
        }

        double c = compute_U(Gdash);

        bool lastrun = 0;
        vector<ppackeredge> HE;
        vector<ppackeredge> href(Gdash->m);

        priority_queue<double> tmpl;
        for (auto it: Gdash->E) {
            auto tmp = new packeredge(it->u, it->v, it->idx, tmpl);
            HE.emplace_back(tmp);
            href[it->idx] = HE.back();
        }

        ppackergraph H = new packergraph(Gdash->n, Gdash->m, HE);

        while (true) {
            //cout << c << " " << lastrun << " " << b/c << "\n";
            double p = min(1.0, b / c);
            for (auto it: Gdash->E) {
                double wt = binom(it->w, p, lli(ceil((1 + eps2) * 12.0 * b)));
                //if(lastrun) cout << it->idx << " " << wt << "\n";
                priority_queue<double> l;
                for (int i = 0; i < wt; i++) l.emplace(0);
                href[it->idx]->l = l;
            }

            H->E = HE;
            pair<double, vector<ptree>> packing = packer(H);

            if (lastrun || p >= 1.0 - cmpeps) {
                res = sample(packing.second, d, Gdash->n);
                break;
            } else {
                for (auto it: packing.second) {
                    for (auto gt: it->E) delete (gt);
                    delete (it);
                }
                if (packing.first >= ((1 - eps3) * b / (2.0 * (1 + eps2)))) {
                    c /= 6.0;
                    lastrun = 1;
                } else c /= 2.0;
            }
        }

        for (auto it: Edash) delete (it);
        for (auto it: HE) delete (it);
        delete (H);
        delete (Gdash);

        return res;
    }

    typedef long long int lli;

// Lemma 8 DS: Slow version
    class lemma8ds {
    public:
        const lli inf = lli(1e12) + 5;

        ptree T;
        int curst = 0, totchain = 0, maxlog;
        double totdelta = 0;
        vector<int> sz, H, inchain, inst, head;
        vector<double> st, lz, val;
        vector<vector<int>> jmp;

        lemma8ds(ptree _T) {
            T = _T;
            curst = 0;
            maxlog = ceil(log(T->n));
            st.clear();
            st.resize(4 * T->n + 10);
            lz.clear();
            lz.resize(4 * T->n + 10);
            sz.clear();
            sz.resize(T->n);
            H.clear();
            H.resize(T->n);
            inchain.clear();
            inchain.resize(T->n);
            val.clear();
            val.resize(T->n);
            inst.clear();
            inst.resize(T->n);
            head.clear();
            head.resize(T->n, -1);
            jmp.clear();
            jmp.resize(T->n, vector<int>(maxlog + 1, -1));

            dfs0(0, -1, 0);
            dfs1(0, -1, totchain++, inf);
            init();
            build(1, 0, curst - 1);
        }

        inline int left(int node) { return (node << 1); }

        inline int right(int node) { return (node << 1) + 1; }

        void build(int node, int L, int R) {
            if (L == R) st[node] = val[L];
            else {
                build(left(node), L, (L + R) / 2), build(right(node), (L + R) / 2 + 1, R);
                st[node] = min(st[left(node)], st[right(node)]);
            }
        }

        void shift(int node, int L, int R) {
            if (!lz[node] || L == R) {
                lz[node] = 0;
                return;
            }

            lz[left(node)] += lz[node], st[left(node)] += lz[node];
            lz[right(node)] += lz[node], st[right(node)] += lz[node];
            lz[node] = 0;
        }

        void upd(int node, int L, int R, int a, int b, double v) {
            if (a > R || b < L) return;
            else if (a <= L && R <= b) st[node] += v, lz[node] += v;
            else {
                shift(node, L, R);
                upd(left(node), L, (L + R) / 2, a, b, v), upd(right(node), (L + R) / 2 + 1, R, a, b, v);
                st[node] = min(st[left(node)], st[right(node)]);
            }
        }

        double qry(int node, int L, int R, int a, int b) {
            if (a > R || b < L) return 0;
            else if (a <= L && R <= b) return st[node];
            else {
                shift(node, L, R);
                return qry(left(node), L, (L + R) / 2, a, b) + qry(right(node), (L + R) / 2 + 1, R, a, b);
            }
        }

        void dfs0(int node, int par, int ht) {
            sz[node] = 1;
            H[node] = ht;
            jmp[node][0] = par;
            for (auto it: T->adj[node]) {
                int child = (it->u == node) ? it->v : it->u;
                if (child != par) {
                    dfs0(child, node, ht + 1);
                    sz[node] += sz[child];
                }
            }
        }

        void dfs1(int node, int par, int chain, double inc) {
            inchain[node] = chain;
            if (head[chain] == -1) head[chain] = node;
            val[curst] = inc;
            inst[node] = curst++;

            pair<int, pair<double, int>> largest = {-1, {-1, -1}};
            for (auto it: T->adj[node]) {
                int child = (it->u == node) ? it->v : it->u;
                if (child != par) largest = max(largest, {sz[child], {it->w, child}});
            }
            if (largest.second.second != -1) dfs1(largest.second.second, node, chain, largest.second.first);
            for (auto it: T->adj[node]) {
                int child = (it->u == node) ? it->v : it->u;
                if (child != par && child != largest.second.second) {
                    dfs1(child, node, totchain++, it->w);
                }
            }
        }


        void init() {
            for (int j = 1; j <= maxlog; j++) {
                for (int i = 0; i < T->n; i++) {
                    if (jmp[i][j - 1] != -1) jmp[i][j] = jmp[jmp[i][j - 1]][j - 1];
                }
            }
        }

        void pathupd(int x, int y, double delta) {
            if (H[x] > H[y]) swap(x, y);
            int origx = x, origy = y;
            for (int i = maxlog; i >= 0; i--) {
                if (H[y] - (1 << i) >= H[x]) y = jmp[y][i];
            }
            if (x == y) {
                int child = origy;
                for (int i = maxlog; i >= 0; i--) {
                    if (H[child] - (1 << i) > H[x]) child = jmp[child][i];
                }
                ancestryupd(origy, child, delta);
                return;
            }
            for (int i = maxlog; i >= 0; i--) {
                if (jmp[x][i] != jmp[y][i]) x = jmp[x][i], y = jmp[y][i];
            }

            ancestryupd(origx, x, delta);
            ancestryupd(origy, y, delta);
        }

        void ancestryupd(int node, int anc, double delta) {
            int cur = node;
            while (inchain[cur] != inchain[anc]) {
                upd(1, 0, curst - 1, inst[head[inchain[cur]]], inst[cur], delta);
                cur = jmp[head[inchain[cur]]][0];
            }
            upd(1, 0, curst - 1, inst[anc], inst[cur], delta);
        }

        void PathAdd(int u, int v, double x) {
            pathupd(u, v, x);
        }

        void NonPathAdd(int u, int v, double x) {
            totdelta += x;
            pathupd(u, v, -x);
        }

        double QueryMinimum() {
            return st[1] + totdelta;
        }

        double QueryEdge(int node) {
            return qry(1, 0, curst - 1, inst[node], inst[node]) + totdelta;
        }
    };

    typedef lemma8ds *plemma8ds;

    // Algorithm 4: initialise and call compute() to solve.
    class mincut {
    public:
        const double inf = 1e12 + 5;

        ptree T;
        pgraph G;
        plemma8ds D;
        int totchain = 0, maxlog;
        vector<int> sz, H, inchain, inst, head, st;
        vector<double> val;
        vector<vector<int>> jmp;

        mincut(pgraph _G, ptree _T) {
            G = _G, T = _T;
            D = new lemma8ds(T);
            maxlog = ceil(log(T->n));
            st.clear();
            sz.clear();
            sz.resize(T->n);
            H.clear();
            H.resize(T->n);
            inchain.clear();
            inchain.resize(T->n);
            val.clear();
            val.resize(T->n);
            inst.clear();
            inst.resize(T->n);
            head.clear();
            head.resize(T->n, -1);
            jmp.clear();
            jmp.resize(T->n, vector<int>(maxlog + 1, -1));

            dfs0(0, -1, 0);
            dfs1(0, -1, totchain++, inf);
            init();
        }

        void clear() {
            delete (D);
        }

        void dfs0(int node, int par, int ht) {
            sz[node] = 1;
            H[node] = ht;
            jmp[node][0] = par;
            for (auto it: T->adj[node]) {
                int child = (it->u == node) ? it->v : it->u;
                if (child != par) {
                    dfs0(child, node, ht + 1);
                    sz[node] += sz[child];
                }
            }
        }

        void dfs1(int node, int par, int chain, double inc) {
            inchain[node] = chain;
            if (head[chain] == -1) head[chain] = node;
            val[int(st.size())] = inc;
            inst[node] = int(st.size());
            st.emplace_back(node);

            pair<int, pair<double, int>> largest = {-1, {-1, -1}};
            for (auto it: T->adj[node]) {
                int child = (it->u == node) ? it->v : it->u;
                if (child != par) largest = max(largest, {sz[child], {it->w, child}});
            }
            if (largest.second.second != -1) dfs1(largest.second.second, node, chain, largest.second.first);
            for (auto it: T->adj[node]) {
                int child = (it->u == node) ? it->v : it->u;
                if (child != par && child != largest.second.second) {
                    dfs1(child, node, totchain++, it->w);
                }
            }
        }


        void init() {
            for (int j = 1; j <= maxlog; j++) {
                for (int i = 0; i < T->n; i++) {
                    if (jmp[i][j - 1] != -1) jmp[i][j] = jmp[jmp[i][j - 1]][j - 1];
                }
            }
        }

        vector<pair<int, int>> pathbreakdown(int x, int y) {
            if (H[x] > H[y]) swap(x, y);
            int origx = x, origy = y;
            for (int i = maxlog; i >= 0; i--) {
                if (H[y] - (1 << i) >= H[x]) y = jmp[y][i];
            }
            if (x == y) {
                int child = origy;
                for (int i = maxlog; i >= 0; i--) {
                    if (H[child] - (1 << i) > H[x]) child = jmp[child][i];
                }
                return ancestrybreakdown(origy, child);
            }
            for (int i = maxlog; i >= 0; i--) {
                if (jmp[x][i] != jmp[y][i]) x = jmp[x][i], y = jmp[y][i];
            }


            vector<pair<int, int>> res = ancestrybreakdown(origx, x);
            vector<pair<int, int>> tmp = ancestrybreakdown(origy, y);
            for (auto it: tmp) res.emplace_back(it);

            return res;
        }

        vector<pair<int, int>> ancestrybreakdown(int node, int anc) {
            vector<pair<int, int>> res;
            int cur = node;

            while (inchain[cur] != inchain[anc]) {
                res.push_back({inst[head[inchain[cur]]], inst[cur]});
                cur = jmp[head[inchain[cur]]][0];
            }

            res.push_back({inst[anc], inst[cur]});
            return res;
        }

        double compute() {
            double res = inf;
            set<int> intree;
            for (auto it: T->E) intree.insert(it->idx);

            //for(auto it: st) cout << it << " ";
            //cout << "\n";

            vector<vector<vector<pedge>>> events(int(st.size()), vector<vector<pedge>>(2));
            for (auto it: G->E) {
                if (intree.find(it->idx) != intree.end()) continue;


                vector<pair<int, int>> path = pathbreakdown(it->u, it->v);

                for (auto gt: path) {
                    int a = min(gt.first, gt.second), b = max(gt.first, gt.second);

                    events[a][0].emplace_back(it);
                    events[b][1].emplace_back(it);
                }

                D->PathAdd(it->u, it->v, it->w);
            }


            for (int i = 0; i < int(st.size()); i++) {
                for (auto it: events[i][0]) {
                    D->NonPathAdd(it->u, it->v, it->w);
                    D->PathAdd(it->u, it->v, -it->w);
                }

                double selfwt = D->QueryEdge(st[i]);
                D->PathAdd(st[i], st[i], inf);
                //cout << i << " " << st[i] << " " << selfwt << " " << D->QueryMinimum() << "\n";
                res = min(res, selfwt + D->QueryMinimum());
                D->PathAdd(st[i], st[i], -inf);

                for (auto it: events[i][1]) {
                    D->NonPathAdd(it->u, it->v, -it->w);
                    D->PathAdd(it->u, it->v, it->w);
                }
            }

            double singlecut = D->QueryMinimum();
            res = min(singlecut, res);

            return res;
        }
    };

    double tester(pgraph G) {
        double d = 0.2;
        vector<ptree> trees = tworespectingtrees(d, G);
//        printf("%d\n", int(trees.size()));

        // Fix wts using edge index
        for (auto it: trees) {
            for (auto gt: it->E) {
                gt->w = G->E[gt->idx]->w;
                //cout << gt->u << " " << gt->v << " " << gt->idx << " " << gt->w << "\n";
            }
            //cout << "\n";
        }

        //vector<ptree> trees;
        //trees.push_back(tree_reader());

        double res = double(1e9) + 5;
        for (auto it: trees) {
            mincut mc(G, it);
            res = min(res, mc.compute());
            mc.clear();
        }

//        for (auto it: trees) {
//            for (auto gt: it->E) delete (gt);
//            delete (it);
//        }
//        for (auto it: G->E) delete (it);
//        delete (G);

        return res;
    }

    pgraph graph_reader(Graph &g) {
        int n, m;
        n = g.get_n(); m = g.edges.size();
        vector<pedge> E;
        unordered_map<t_node , int> who;
        int w = 0;
        for (int i = 0; i < m; i++) {
            int u, v;
            double c;
//            scanf("%d%d%lf", &u, &v, &c);
            //u--, v--;
            u = g.edges[i].from; v = g.edges[i].to;
            c = g.edges[i].cnt;
            if(!who.count(u))   who[u] = w++;
            if(!who.count(v))   who[v] = w++;
            u = who[u]; v = who[v];
            pedge tmp = new edge(u, v, i, c);
            E.emplace_back(tmp);
        }

        pgraph res = new graph(n, m, E);
        return res;
    }
}

namespace KargerLinearMincut {
    using namespace KargerLinearMincut_GH;

    int mincut(Graph &g) {
        pgraph pg = graph_reader(g);
        auto ans = tester(pg);
        return (int)round(ans);
    }
}