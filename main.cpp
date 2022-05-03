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
    t_double eps = 0.2;
    t_double delta = 0.05;

    for(auto p: ps) {
        profiler.reset();

        Graph my_g = g;
        my_g.p = p;

        cout << scientific << "p = " << p << endl;

        profiler.start("unrel");
        auto ans_unrel = median_trick([&]() { return compute_unreliability(my_g, eps); }, 4, delta);
        cout << scientific << "answer: " << ans_unrel << endl;
        profiler.stop("unrel");

//        my_g = g;
//        my_g.p = p;
//        profiler.start("brut");
//        auto ans_brute = brute_unreliability(my_g);
//        cout << scientific << "answer brute: " << ans_brute << endl;
//        profiler.stop("brut");

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

    /*

    eps = 0.2, delta = 0.05

    repeat 3 times each experiment

    Graph types:
    complete graphs: n = 5 - 25

    P values:
    0.001 - 0.999

    Measurements:
    total_time, answer

     */

    json j;

    vector<Graph> graphs;

//    for(int n = 10; n <= 200; n += 10) {
//        graphs.push_back(Generator::complete_graph(n));
//    }

//    for(int n = 2; n <= 10; n++)
//        graphs.push_back(Generator::grid(n));

    for(int i = 0; i < 10; i++)
        graphs.push_back(Generator::erdos_renyi(50, 0.05 * Random::get_int(10, 20)));

    const t_double eps = 0.2;
    const t_double delta = 0.05;
    const int TRIES = 1;

    Random::init_predictable(true);
    int i = 0;
    for(auto g: graphs) {
        int n = g.get_n();
        string num_id = to_string(i);
        num_id = string(3 - (int)num_id.size(), '0') + num_id;
        string graph_id = "er_" + num_id;
        i++;

        Profiler graph_profiler;
        graph_profiler.reset();

        t_double p_delta = 0.01;

        vector<t_double> ps;
        for (t_double p = p_delta; p < 1; p += 2 * p_delta) {
            ps.push_back(p);
            if(p > 0.25 && p < 0.9) p += p_delta;
        }

        for (auto p: ps) {
            string subgraph_id = graph_id + "_p_" + to_string(p);
            g.p = p;

            j[graph_id][subgraph_id]["p"] = p;

            Profiler sum_profiler;
            sum_profiler.reset();

            for(int t = 0; t < TRIES; t++) {
                string try_id = "try " + to_string(t + 1);

                profiler.reset();
                profiler.start("time_unrel");
                auto answer = median_trick([&]() { return compute_unreliability(g, eps); }, 4, delta);
                profiler.stop("time_unrel");

                sum_profiler.update(profiler);

                j[graph_id][subgraph_id][try_id]["answer"] = answer;
                j[graph_id][subgraph_id][try_id]["time_unrel"] = profiler.time_spent["time_unrel"];

                cout << subgraph_id << ", " << try_id << ": Done!" << endl;
            }

            j[graph_id][subgraph_id]["average_time"] = sum_profiler.time_spent["time_unrel"] / TRIES;

            graph_profiler.update(sum_profiler);
        }

        j[graph_id]["average_time"] = graph_profiler.time_spent["time_unrel"] / (TRIES * ps.size());

        ofstream json_writer("experiment_results/er_diff_p_exp.json");
        json_writer << setw(4) << j << endl;
    }
}

void experiment(string graph_id, Graph &ig, vector<t_double> ps) {
    cout << graph_id << endl;
    const t_double eps = 0.2;
    const t_double delta = 0.05;
    for(auto p: ps) {
        cout << scientific << "p = " << p << endl;
        Graph g = ig;
        g.p = p;
        profiler.reset();
        profiler.start("time_unrel");
        auto ans_unrel = median_trick([&]() { return compute_unreliability(g, eps); }, 4, delta);
        cout << scientific << "answer: " << ans_unrel << endl;
        profiler.stop("time_unrel");
        cout << "time: " << profiler.time_spent["time_unrel"] / 1000 << endl;
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
    1e-2, 1e-1, 0.25, 0.5, 0.9

    Measurements:
    brute
    contract
    mincut

     */

    json j;
    ifstream json_reader("experiment_results/mincut_exp.json");
    json_reader >> j;
    json_reader.close();

    const t_double eps = 0.2;
    const t_double delta = 0.05;
    const vector<t_double> ps = {1e-2, 0.1, 0.25, 0.5, 0.9};
    const int TRIES = 3;

    auto graphs = mincut_experiments_graphs_generator();

    Random::init_predictable(true);
    for(auto [graph_id, _g]: graphs) {
        if(j.contains(graph_id))    continue;

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

void table_format() {
    json j;
    ifstream json_reader("experiment_results/mincut_exp.json");
    json_reader >> j;

    const vector<string> order = {
            "brute", "contract",
            "ek-mincut", "din-mincut",
            "ks-mincut"
    };

    auto graphs_vec = mincut_experiments_graphs_generator();
    map<string, Graph> graphs;
    for(auto [id, g]: graphs_vec)
        graphs[id] = g;

    unordered_map<string, long> tot_time;
    const vector<string> procedures = {
            "brute", "contract",
            "ek-mincut",
            "din-mincut",
            "ks-mincut"
    };
    const int T = 100;
    for(int t = 0; t < T; t++) {
        string num_id = to_string(t);
        num_id = string(2 - num_id.size(), '0') + num_id;
        string graph_id = "er_" + num_id;

        cout << t << " & ";
        cout << graphs[graph_id].get_n() << " & ";
        cout << graphs[graph_id].get_m() << " & ";

        for(auto who: procedures)
            cout <<  (long)j[graph_id]["profiling"][who] / 1000 << " & ";
        cout << R"(0 \\ \hline)" << endl;
    }

//    for(auto who: procedures) {
//        cout << who << ": " << tot_time[who] / (T * 1000) << endl;
//    }

//    for(int n = 2; n <= 8; n++) {
//        string id = "grid_" + to_string(n);
//        cout << n;
//        for(auto procedure: order)
//            cout << " & " << (long)j[id]["profiling"][procedure] / 1000;
//
//        cout << " & 0 \\\\\n";
//    }
}

void epsilon_comparison() {
    json j;

    const vector<t_double> epses = {0.25, 0.2, 0.1, 0.05, 0.01};
    const vector<t_double> ps = {0.001, 0.01, 0.1, 0.2, 0.25, 0.5, 0.75, 0.9};
    const t_double delta = 0.05;

    Random::init(true);
    for(int t = 0; t < 1000; t++) {
        int n = Random::get_int(5, 10);
        int m = Random::get_int(n - 1, n * (n - 1) / 2);
        Graph g = Generator::random(n, m);

        cout << "Graph " << t << ": n = " << n << ", m = " << m << endl;

        string num_id = to_string(t);
        num_id = string(3 - num_id.size(), '0') + num_id;
        string graph_id = num_id;

        j[graph_id]["n"] = n;
        j[graph_id]["m"] = m;

        for(auto p: ps) {
            string p_id = to_string(p);
            g.p = p;

            profiler.reset();
            profiler.start("brute");
            auto ans_brute = brute_unreliability(g);
            profiler.stop("brute");

            j[graph_id][p_id]["brute"]["answer"] = ans_brute;
            j[graph_id][p_id]["brute"]["time"] = profiler.time_spent["brute"];
            cout << "Brute " << graph_id << ", p = " << p << ": Done!" << endl;


            for (auto eps: epses) {
                string now_id = to_string(eps);
                Graph gg = g;

                profiler.reset();
                profiler.start("unrel_time");
                auto answer = median_trick([&]() { return compute_unreliability(gg, eps); }, 4, delta);
                profiler.stop("unrel_time");

                j[graph_id][p_id][now_id]["answer"] = answer;
                j[graph_id][p_id][now_id]["time"] = profiler.time_spent["unrel_time"];
                t_double rap = answer / ans_brute;
                t_double act_eps = rap - 1.0;
                j[graph_id][p_id][now_id]["rap"] = rap;
                j[graph_id][p_id][now_id]["act_eps"] = act_eps;

                cout << "Unrel " << graph_id << ", p = " << p << ", eps = " << eps << ": Done!" << endl;
            }
        }

        ofstream json_writer("experiment_results/empirical_eps.json");
        json_writer << setw(4) << j << endl;

    }
}

void empiric_epsilon_stats() {
    freopen("logs/empirical_eps.log", "w", stdout);
    ifstream reader("experiment_results/empirical_eps.json");
    json j;
    reader >> j;
    reader.close();

    unordered_map<string, t_double> max_eps;
    unordered_map<string, string> max_graph_id;
    unordered_map<string, vector<t_double>> act_epses;
    const vector<t_double> ps = {0.001, 0.01, 0.1, 0.2, 0.25, 0.5, 0.75, 0.9};
    const vector<t_double> epses = {0.25, 0.2, 0.1, 0.05, 0.01};

    for (auto it = j.begin(); it != j.end(); it++) {
        string graph_id = it.key();
        for(auto p: ps) {
            string p_id = to_string(p);
            for(auto eps: epses) {
                string eps_id = to_string(eps);
                t_double act_eps = (t_double)j[graph_id][p_id][eps_id]["act_eps"];
                act_eps = fabs(act_eps);
                act_epses[eps_id].push_back(act_eps);
                if(max_eps[eps_id] < act_eps) {
                    max_eps[eps_id] = act_eps;
                    max_graph_id[eps_id] = graph_id;
                }
//                max_eps[eps_id] = max(max_eps[eps_id], fabs(act_eps));
            }
        }
    }

    for(auto [eps, _]: max_eps) {
        cout << "eps: " <<  eps << ":" << endl;
        t_double eps_val;
        for(auto e: epses)
            if(to_string(e) == eps)
                eps_val = e;

        auto epses_now = act_epses[eps];
        sort(begin(epses_now), end(epses_now));

        for(int i = 0; i < epses_now.size(); i++) {
            cout << scientific << i << ": " << epses_now[i] << " , rap: " << epses_now[i] / eps_val << endl;
        }

        cout << endl;
    }

    for(auto [eps, act_eps]: max_eps)
        cout << scientific << "eps = " << eps << ", max_act_eps = " << act_eps << ", graph_id = " << max_graph_id[eps] << endl;
}

void parallel_experiments() {
    Graph geant = read_graph("graphs/geant2009.in");
    Graph k50 = Generator::k4(50);
    Graph k100 = Generator::k4(100);
    vector<t_double> ps = {0.25, 0.1, 0.01, 0.001};

    experiment("geant", geant, ps);
    experiment("k50", k50, ps);
    experiment("k100", k100, ps);
}

void brute_comparison() {
    Random::init_predictable(true);
    for (int t = 0; t < 100000; t++) {
        int n = Random::get_int(5, 9);
        int m = Random::get_int(n - 1, n * (n - 1) / 2);
        Graph g = Generator::random(n, m);
        t_double p = t_double(0.001) * Random::get_int(1, 999);
        g.p = p;

        cout << "Graph " << t << ": n = " << n << ", m = " << m << ", p = " << p << endl;

        profiler.reset();
        profiler.start("unrel");
        auto ans_unrel = brute_unreliability(g);
        profiler.stop("unrel");

        profiler.start("rel");
        auto ans_rel = brute_reliability(g);
        profiler.stop("rel");

        auto tot = ans_rel + ans_unrel;
        assert( fabs(tot - 1.0) < 0.00000001 );

        cout << scientific << "ans_unrel: " << ans_unrel << endl;
        cout << scientific << "ans_rel: " << ans_rel << endl;
        cout << scientific << "tot: " << ans_rel + ans_unrel << endl;

        profiler.print();

        cout << endl;
    }
}

void rel_verify() {
    const t_double eps = 0.2;
    const t_double delta = 0.05;

    Random::init_predictable(true);
    for (int t = 0; t < 1000000; t++) {
        int n = Random::get_int(3, 6);
        int m = Random::get_int(n - 1, n * (n - 1) / 2);
        t_double p = t_double(0.001) * Random::get_int(1, 950);
        Graph g = Generator::random(n, m);
        g.p = p;
        DiGraph dg(g);

        cout << "Graph " << t << ": n = " << n << ", m = " << m << ", p = " << p << endl;

        profiler.reset();
        profiler.start("rel");
        auto ans_rel = median_trick([&]() { return compute_reliability(dg, eps); }, 4, delta);
        profiler.stop("rel");

        cout << scientific << "ans_rel: " << ans_rel << endl;

        profiler.start("brute");
        auto ans_brute = brute_reliability(g);
        profiler.stop("brute");

        cout << scientific << "ans_brute: " << ans_brute << endl;
        t_double act_eps = (ans_rel / ans_brute) - 1.0;
        cout << scientific << "epsilon: " << act_eps << endl;

        profiler.print();

        cout << endl;

        if(fabs(act_eps) > 2 * eps) {
            ofstream out("logs/rel_verify.log");
            print_graph(g, out);
            assert(0);
        }
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

//    mincut_experiments();

//    table_format();

//    epsilon_comparison();
//    empiric_epsilon_stats();

//    parallel_experiments();

//    brute_comparison();

    rel_verify();

    return 0;
}