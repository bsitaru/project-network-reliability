#ifndef NETWORK_RELIABILITY_UNRELIABILITY_TYPES_HPP
#define NETWORK_RELIABILITY_UNRELIABILITY_TYPES_HPP

typedef long double t_double;
typedef int t_node;
struct t_edge {
    t_node from, to;
    int id, cnt;

    int other(int x) const { return from ^ to ^ x; }
};

#endif