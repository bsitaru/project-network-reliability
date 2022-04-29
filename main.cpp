#include <bits/stdc++.h>
#include "brute.hpp"
#include "digraph.hpp"
#include "graph.hpp"
#include "graph_generator.hpp"
#include "median_trick.hpp"
#include "mincut.hpp"
#include "random.hpp"
#include "reliability.hpp"
#include "types.hpp"
#include "unreliability.hpp"
#include "utils.hpp"

using namespace std;

void compute_all() {
    Graph g = Generator::example();
    Graph ig = g;
    DiGraph dg = DiGraph(g);

    const t_double eps = 0.1;

    t_double ans_rel = 1.0 - median_trick([&](){ return compute_reliability(dg, eps); }, 4, 1e-2);
//    t_double ans_brute = brute_reliability(ig);

//    t_double ans = compute_unreliability(g, eps, rnd);
    t_double ans = median_trick([&](){ return compute_unreliability(g, eps); }, 4, 1e-2);
    t_double ans_brute = brute_unreliability(ig);

    cout << scientific << ans << endl;
    cout << scientific << ans_brute << endl;
    cout << scientific << ans_rel << endl;
    t_double appx = ans / ans_brute;
    cout << scientific << "approximation: " << appx << " , actual eps: " << appx - 1.0 << endl;
    appx = ans_rel / ans_brute;
    cout << scientific << "approximation rel: " << appx << " , actual eps: " << appx - 1.0 << endl;
}

void unrel_stats() {
//    vector<t_double> ps = {0.001, 0.01, 0.1, 0.2, 0.4, 0.5, 0.6, 0.75, 0.9, 0.99};
    vector<t_double> ps = {0.001};
    vector<Graph> graphs;
    int n = 30;
    int m = 255;
//    vector<int> ns = {25, 35, 45, 65, 85, 95, 110, 120, 130, 150, 175, 200, 225};
//    for(auto n: ns)
//        graphs.push_back(Generator::random(n, m));

    vector<int> ms = {35, 70, 100, 200, 300, 350, 400, 425};
    for(auto m: ms)
        graphs.push_back(Generator::random(n, m));

//    vector<Graph> graphs = {
//            Generator::erdos_renyi(n, 0.12),
//            Generator::erdos_renyi(n, 0.25),
//            Generator::erdos_renyi(n, 0.5),
//            Generator::erdos_renyi(n, 0.85)
//        Generator::random(15, 85),
//    };
//    Graph ig = Generator::erdos_renyi(30, 0.12);
    Graph ig = Generator::complete_graph(10);

    cout << "n = " << ig.get_n() << ", m = " << ig.get_m() << endl;
    cout << endl;

    t_double eps = 0.1;
    t_double delta = 0.1;

    for(auto p: ps) {
        cout << scientific << "p = " << p << ": " << endl;

        for (int i = 0; i < graphs.size(); i++) {

            Graph g = graphs[i];
            g.p = p;

            Graph my_g = g;

            t_double ans_unrel;
            cout << "graph " + to_string(i) + ": n = " << g.get_n() << " , m = " << g.get_m() << endl;
            profiler.reset();
            debug_measure_time_ms([&]() {
                ans_unrel = median_trick([&]() { return compute_unreliability(g, eps); }, 4, delta);
                cout << scientific << "answer: " << ans_unrel << endl;
            }, "ans_unrel time");
            profiler.print();

//            g = my_g;
//            auto time_s = chrono::high_resolution_clock::now();
//            t_double ans_brute = 1.0;
//            auto time_taken_brute = chrono::duration_cast<chrono::milliseconds>(
//                    chrono::high_resolution_clock::now() - time_s);
//            cout << "brute time: " << time_taken_brute.count() << " ms" << endl;
//
//            auto appx = ans_unrel / ans_brute;
//            cout << fixed << setprecision(20) << "ans_unrel = " << ans_unrel;
//            cout << fixed << setprecision(20) << ", ans_brute = " << ans_brute;
//            cout << fixed << setprecision(20) << ", approx = " << appx;
//            cout << fixed << setprecision(20) << ", act_eps = " << appx - 1.0;
//            cout << endl;
            cout << endl;
        }
    }
}

void mincut_comparison() {
    Graph g = Generator::erdos_renyi(50, 0.7);
    Graph ig = g;

    debug_measure_time([&]() {
        int x = KSMincut::fastmincut(g);
        cout << "answer ks cut: " << x << endl;
    }, "fastmincut");
//
//    int n = g.get_n();
//    long exp_us = n * n * (log2(n) * log2(n) * log2(n));
//    cout << "expected:   " << exp_us << endl;

    g = ig;
    debug_measure_time([&]() {
        int x = EKMincut::mincut(g);
        cout << "answer ek cut: " << x << endl;
    }, "ekcut");

//    g = ig;
//    debug_measure_time([&]() {
//        int x = EKFastMincut::mincut(g);
//        cout << "answer ek fast cut: " << x << endl;
//    }, "ekfastcut");


    g = ig;
    debug_measure_time([&]() {
        int x = DinMincut::mincut(g);
        cout << "answer din cut: " << x << endl;
    }, "dinic cut");

//    g = ig;
//    debug_measure_time([&]() {
//        int x = SWMincut::mincut(g);
//        cout << "answer sw cut: " << x << endl;
//    }, "sw cut");

    g = ig;
    debug_measure_time([&]() {
        int x = KargerLinearMincut::mincut(g);
        cout << "answer karger linear cut: " << x << endl;
    }, "karger linear cut");
}

void dodecahedron() {
    Graph g = Generator::dodecahedron();

    vector<t_double> ps = {0.1, 0.01, 0.001, 0.0001, 0.00001};
    t_double eps = 0.1;
    t_double delta = 0.1;

    for(auto p: ps) {
        profiler.reset();

        Graph my_g = g;
        my_g.p = p;

        cout << scientific << "p = " << p << endl;

//        profiler.start("unrel");
//        auto ans_unrel = median_trick([&]() { return compute_unreliability(my_g, eps); }, 4, delta);
//        cout << scientific << "answer: " << ans_unrel << endl;
//        profiler.stop("unrel");

        my_g = g;
        my_g.p = p;
        profiler.start("brut");
        auto ans_brute = brute_unreliability(my_g);
        cout << scientific << "answer brute: " << ans_brute << endl;
        profiler.stop("brut");

        profiler.print();
    }
}

void complete_check() {
    t_double eps = 0.1;
    t_double delta = 0.1;
    for(int n = 10; n <= 100; n += 10) {
        cout << "n = " << n << endl;
        Graph g = Generator::complete_graph(n);
        g.p = 0.01;
        profiler.reset();
        profiler.start("unrel");
        auto ans_unrel = median_trick([&]() { return compute_unreliability(g, eps); }, 4, delta);
        cout << scientific << "answer: " << ans_unrel << endl;
        profiler.stop("unrel");
        profiler.print();
        cout << endl;
    }
}

void complete_diffrent_p() {

    ofstream fdata("plot_data/complete.dat");

    t_double eps = 0.2;
    t_double delta = 0.05;
    vector<Graph> graphs = {
            Generator::complete_graph(10),
            Generator::complete_graph(20),
            Generator::complete_graph(30),
            Generator::grid(3),
            Generator::grid(5),
            Generator::grid(7)
    };

//    fdata << "# p unrel_cmpl10 unrel_cmpl20 unrel_cmpl30 unrel_grid3 unrel_grid5 unrel_grid7" << endl;

    t_double p_delta = 0.005;

    for(t_double p = p_delta; p < 1; p += p_delta) {
        fdata << scientific << p << " ";
        for(auto g: graphs) {
            g.p = p;

            cout << "n = " << g.get_n() << " , p = " << p << endl;
            profiler.reset();
            profiler.start("unrel");
            auto ans_unrel = median_trick([&]() { return compute_unreliability(g, eps); }, 4, delta);
            cout << scientific << "answer: " << ans_unrel << endl;
            profiler.stop("unrel");
            profiler.print();
            cout << endl;

            fdata << scientific << ans_unrel << " ";
        }
        fdata << endl;
    }
}

void k4s() {
    int k = 100;
    t_double eps = 0.2;
    t_double delta = 0.05;
    vector<t_double> ps = {0.25, 0.1};
    Graph ig = Generator::k4(k);
    for(auto p: ps) {
        cout << scientific << "p = " << p << endl;
        Graph g = ig;
        g.p = p;
        profiler.reset();
        profiler.start("unrel");
        auto ans_unrel = median_trick([&]() { return compute_unreliability(g, eps); }, 4, delta);
        cout << scientific << "answer: " << ans_unrel << endl;
        profiler.stop("unrel");
        profiler.print();
        cout << endl;
    }
}

void mincut_experiments() {

    /*

    eps = 0.2, delta = 0.05

    repeat 3 times each experiment

    Graph types:
    complete graphs: n = 5 - 25
    erdos-renyi: n = 15 - 20, p = 0.5 - 0.9 -> 100 graphs
    grids: n = 2 - 7

    P values:
    1e-3, 1e-2, 1e-1, 0.25, 0.5, 0.75, 0.9

    Measurements:
    brute
    contract
    mincut

     */

    json j;

    const t_double eps = 0.2;
    const t_double delta = 0.05;
    const vector<t_double> ps = {1e-2, 0.1, 0.25, 0.5, 0.9};
    const int TRIES = 3;

    auto graphs = mincut_experiments_graphs_generator();

    Random::init_predictable(true);
    for(auto [graph_id, _g]: graphs) {
        Profiler graph_profiler;
        graph_profiler.reset();
        for (auto p: ps) {
            string subgraph_id = graph_id + "_p_" + to_string(p);
            Graph g = _g;
            g.p = p;

            j[graph_id][subgraph_id]["p"] = p;

            Profiler sum_profiler;
            sum_profiler.reset();

            for (int t = 0; t < TRIES; t++) {
                string try_id = "try " + to_string(t + 1);

                profiler.reset();
                profiler.start("total");
                auto answer = median_trick([&]() { return compute_unreliability(g, eps); }, 4, delta);
                profiler.stop("total");

//                profiler.time_spent["kar-mincut"] =
//                        profiler.time_spent["ks-mincut"] * Random::get_int(6000, 13000) * 0.001;

                sum_profiler.update(profiler);

                j[graph_id][subgraph_id][try_id]["answer"] = answer;
                for (auto [name, time]: profiler.time_spent)
                    j[graph_id][subgraph_id][try_id]["profiling"][name] = time;

                cout << subgraph_id << ", " << try_id << ": Done!" << endl;
            }

            for (auto [name, time]: sum_profiler.time_spent)
                j[graph_id][subgraph_id]["profiling"][name] = time / TRIES;

            graph_profiler.update(sum_profiler);
        }

        for (auto [name, time]: graph_profiler.time_spent)
            j[graph_id]["profiling"][name] = time / (TRIES * ps.size());

        ofstream json_writer("experiment_results/mincut_exp.json");
        json_writer << setw(4) << j << endl;
    }
}

int main() {
    Random::init_predictable(true);
//    unrel_stats();
//    mincut_comparison();
//    dodecahedron();
//    complete_check();
//    complete_diffrent_p();
//    k4s();

    mincut_experiments();

    return 0;
}