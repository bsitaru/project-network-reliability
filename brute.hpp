#ifndef NETWORK_RELIABILITY_UNRELIABILITY_BRUTE_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_BRUTE_HPP

#include <bits/stdc++.h>
#include "graph.hpp"
#include "types.hpp"

using namespace std;

t_double brute_rel(Graph g) {
    // TODO: use UnionFindUndo to improve time
    t_double ans_rel = 0.0;

    UnionFind uf(g.nodes);

    function<void(int, t_double, UnionFind &)> bck = [&](int i, t_double w, UnionFind &uf) {
        if (uf.components == 1) {
            ans_rel += w;
            return;
        }

        while (i < g.edges.size() && uf.find(g.edges[i].from) == uf.find(g.edges[i].to))
            i++;

        if (i >= g.edges.size()) return;

        if (uf.components - 1 > (int) g.edges.size() - i)
            return;

        t_double p = power(g.p, g.edges[i].cnt);

        bck(i + 1, w * p, uf);

        UnionFind uf2 = uf;
        uf2.unite(g.edges[i].from, g.edges[i].to);
        bck(i + 1, w * (t_double(1.0) - p), uf2);
    };

    bck(0, 1.0, uf);

//    ans_unrel = 1.0 - ans_rel;

//    return reliability ? ans_rel : ans_unrel;
    return ans_rel;
}

t_double brute_unrel(Graph g) {
        int m = (int)g.edges.size();
        t_double ans_unrel = 0.0;

        UnionFindUndo uf(g.nodes);

        function<void(int, t_double)> bck = [&](int i, t_double w) {
            while (i < g.edges.size() && uf.find(g.edges[i].from) == uf.find(g.edges[i].to))
                i++;

            if(uf.components - 1 > m - i) {
                ans_unrel += w;
                return;
            }

            if (i >= m) {
                return;
            }

            t_double p = power(g.p, g.edges[i].cnt);

            bck(i + 1, w * p);

            if (uf.components > 2) {
                uf.unite(g.edges[i].from, g.edges[i].to);
//        bck(i + 1, w * (t_double(1.0) - p), uf2);
                bck(i + 1, w - w * p);
                uf.undo();
            }
        };

        bck(0, 1.0);

        return ans_unrel;
    }


t_double brute_reliability(Graph g) {
    return brute_rel(g);
}

t_double brute_unreliability(Graph g) {
    return brute_unrel(g);
}


#endif //NETWORK_RELIABILITY_UNRELIABILITY_BRUTE_HPP
