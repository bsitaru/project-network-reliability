#ifndef NETWORK_RELIABILITY_UNRELIABILITY_UTILS_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_UTILS_HPP

#include "types.hpp"

#include "libraries/json.hpp"

using json = nlohmann::json;

t_double power(t_double x, int e, t_double ans = 1.0) {
    if(e == 0)  return ans;
    return power(x * x, e >> 1, ans * ( (e & 1) ? x : 1.0 ));
}

t_double n_root(t_double x, int n) {
    t_double l = 0.0, r = 1.0;
    int steps = 100;
    for(int i = 0; i < steps; i++) {
        t_double m = (l + r) / 2.0;
        if(power(m, n) <= x)    l = m;
        else    r = m;
    }
    return l;
}

struct UnionFind {
    unordered_map<t_node, t_node> uf;
    int components;

    UnionFind(const vector<t_node> &nodes) {
        components = nodes.size();
        for (auto node: nodes) uf[node] = node;
    }

    t_node find(t_node x) {
        if (uf[x] == x) return x;
        return uf[x] = find(uf[x]);
    }

    void unite(t_node x, t_node y) {
        t_node fx = find(x), fy = find(y);
        if(fx != fy) {
            uf[fy] = fx;
            components--;
        }
    }
};

struct UnionFindUndo {
    struct Node {
        t_node father;
        int sz;
    };
    unordered_map<t_node, Node> uf;
    int components;
    vector< pair<t_node, t_node> > stack;

    UnionFindUndo(const vector<t_node> &nodes) {
        components = nodes.size();
        for (auto node: nodes) uf[node] = {node, 1};
        stack.reserve(nodes.size());
    }

    t_node find(t_node x) {
        if (uf[x].father == x) return x;
        return find(uf[x].father);
    }

    void unite(t_node x, t_node y) {
        t_node fx = find(x), fy = find(y);
        if(fx != fy) {
            if(uf[fy].sz > uf[fx].sz)
                swap(fx, fy);

            stack.push_back({fx, fy});

            uf[fy].father = fx;
            uf[fx].sz += uf[fy].sz;
            components--;
        }
    }

    void undo() {
        assert(!stack.empty());

        auto [fx, fy] = stack.back();
        stack.pop_back();

        uf[fy].father = fy;
        uf[fx].sz -= uf[fy].sz;
        components++;
    }
};

void debug_measure_time(function<void()> fn, string text, bool print = true) {
    auto t_start = chrono::high_resolution_clock::now();
    fn();
    auto time_taken = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - t_start);
    if(print)   cout << text << ": " << time_taken.count() << " us" << endl;
}

void debug_measure_time_ms(function<void()> fn, string text, bool print = true) {
    auto t_start = chrono::high_resolution_clock::now();
    fn();
    auto time_taken = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t_start);
    if(print)   cout << text << ": " << time_taken.count() << " ms" << endl;
}

struct Profiler {
    using time_stamp = chrono::time_point<chrono::system_clock>;

    unordered_map<string, long> time_spent;
    unordered_map<string, time_stamp> time_stamps;

    void reset() {
        time_spent.clear();
    }

    void print() {
        vector< pair<long, string> > prs;
        for(auto [x, y]: time_spent)
            prs.emplace_back(y, x);

        cout << "Profiling: " << endl;
        for(auto [time, who]: prs)
            cout << fixed << setprecision(3) << who << ": " << (double)time * 0.001 << " ms" << endl;
    }

    void profile(const string& who, function<void()> fn) {
        auto t_start = chrono::high_resolution_clock::now();
        fn();
        auto time_taken = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - t_start);
        time_spent[who] += time_taken.count();
    }

    void start(const string& who) {
        time_stamps[who] = chrono::high_resolution_clock::now();
    }

    void stop(const string& who) {
        auto now = chrono::high_resolution_clock::now();
        auto time_taken = chrono::duration_cast<chrono::microseconds>(now - time_stamps[who]);
        time_spent[who] += time_taken.count();
    }

    void update(Profiler p) {
        for(auto [name, time]: p.time_spent)
            time_spent[name] += time;
    }
};

static Profiler profiler;

#endif //NETWORK_RELIABILITY_UNRELIABILITY_UTILS_HPP
