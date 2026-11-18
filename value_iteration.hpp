// value_iteration.hpp
// Weights on the states of Q, and exact bounds on the growth rate that they certify.
//
//   StateSet<P> Q = generateStates<P>();
//   Bounds b = valueIteration(Q, 1000, 10);       // growth rate in [b.lower, b.upper] / 10^6
//
// Everything is an integer: weights are at most SCALE = 10^8, and ratios are fractions with denominator
// RATIO = 10^6, rounded outward, so every bound and every comparison below is exact.
#pragma once
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "generate_states.hpp"

constexpr uint64_t RATIO = 1000000;  // ratios are fractions with denominator 10^6

// A ratio given in units of 1 / RATIO, as a decimal number: decimal(1306541) is "1.306541" (six digits,
// since RATIO = 10^6).
inline std::string decimal(uint64_t v) {
    const std::string fraction = std::to_string(v % RATIO);
    return std::to_string(v / RATIO) + "." + std::string(6 - fraction.size(), '0') + fraction;
}

// ---------------------------------------------------------------------------------------------------
// T_x(s): the smallest total weight of the successors that a list of three letters can leave,
// computed as in the article. A letter that completes a square contributes 0.
//   * Z_0(A), the last letter, is never allowed and contributes 0, so some minimizing triple contains it.
//     The two other letters are either Z_1(A) and an outside letter, or two outside letters on disjoint
//     lags (two new letters, X = 0, are allowed).
//   * If the state has no outside type at all no list is admissible: we return 0
//     (this never happens for a reachable state and this could be proven, but puting weight 0 "safe").
// ---------------------------------------------------------------------------------------------------
template <int P>
uint32_t cheapestList(const StateSet<P>& Q, const State<P>& s) {
    static thread_local std::vector<uint32_t> types;
    static thread_local std::vector<std::pair<uint32_t, uint32_t>> options;  // (weight of the successor, X)

    outsideTypes(s, types);
    if (types.empty()) return 0;
    
    options.clear();
    const uint32_t weightOfOne = Q.extensionWeight(s, s.lagsWith(State<P>::One));
    uint32_t minimumOutsideWeight = StateSet<P>::SCALE;
    for (uint32_t X : types) {
        const uint32_t weight = Q.extensionWeight(s, X);
        options.push_back({weight, X});
        minimumOutsideWeight = std::min(minimumOutsideWeight, weight);
    }
    if (weightOfOne <= minimumOutsideWeight)
        return weightOfOne + minimumOutsideWeight;
    std::sort(options.begin(), options.end());  // cheapest first

    // Z_0, Z_1 and the cheapest outside letter. Z_1 is never empty: lag 2 always has the label One.
    uint32_t best = weightOfOne + minimumOutsideWeight;    // Z_0 and two outside letters: for each i, the first disjoint j >= i is the best partner; stop when no
    // pair can win.
    for (size_t i = 0; i < options.size() && 2 * uint32_t(options[i].first) < best; ++i)
        for (size_t j = i; j < options.size() && uint32_t(options[i].first) + options[j].first < best; ++j) {
            const uint32_t X = options[i].second, Y = options[j].second;
            if ((X & Y) == 0 && (i != j || X == 0)) {
                best = uint32_t(options[i].first) + options[j].first;
                break;
            }
        }
    return uint32_t(best);  // at most 2 * SCALE < 2^32
}

// ---------------------------------------------------------------------------------------------------
// Value iteration: y = T_x, then x = y normalized so that the largest weight is SCALE. The growth rate lies
// between the smallest and the largest ratio y/x over the states with x > 0. Bounds holds them as fractions
// with denominator RATIO, the lower one rounded down and the upper one rounded up. When upper - lower is
// at most the tolerance (in units of 1 / RATIO), or at the last iteration, the iteration stops before
// normalizing: Q then holds the weights x the bounds refer to. Since T_x was computed exactly,
// RATIO * T_x(A) >= lower * x(A) for every state: with x > 0 at the start state, these weights are a
// certificate for alpha = lower / RATIO, and no separate check is needed.
// ---------------------------------------------------------------------------------------------------
struct Bounds {
    uint64_t lower = 0, upper = 0;  // in units of 1 / RATIO
    int iterations = 0;
};

template <int P>
Bounds valueIteration(StateSet<P>& Q, int maxIterations, uint64_t tolerance) {
    const size_t n = Q.size();
    Bounds b;
    for (b.iterations = 1; b.iterations <= maxIterations; ++b.iterations) {
// The following loop represents most of the computation time
// we use openmp to parallelize it, with dynamic scheduling and a chunk size of 1024
#pragma omp parallel for schedule(dynamic, 1024)
        for (size_t id = 0; id < n; ++id) 
            Q.setNextWeight(id, cheapestList(Q, Q[id]));

        // lower only has to be at most every ratio, so any starting value is valid; 3 * RATIO is above
        // the growth rate. The products y * RATIO < 3 * 10^14 fit in 64 bits.
        b.lower = 3 * RATIO;
        b.upper = 0;
        uint64_t maxY = 0;
        for (size_t id = 0; id < n; ++id) {
            const uint64_t x = Q.weight(id), y = Q.nextWeight(id);
            maxY = std::max(maxY, y);
            if (x == 0) continue;
            b.lower = std::min(b.lower, y * RATIO / x);            // rounded down
            b.upper = std::max(b.upper, (y * RATIO + x - 1) / x);  // rounded up
        }
        if (b.iterations % 10 == 0) 
            std::cout << "  iteration " << b.iterations << ": growth rate in [" << decimal(b.lower) << ", "
                      << decimal(b.upper) << "]" << std::endl;
        // Stop before normalizing, so that Q keeps the weights the bounds refer to (maxY = 0: every weight vanished).
        if (maxY == 0 || b.upper - b.lower <= tolerance || b.iterations == maxIterations) break;
        for (size_t id = 0; id < n; ++id)
            Q.setWeight(id, uint32_t(Q.nextWeight(id) * uint64_t(StateSet<P>::SCALE) / maxY));
    }
    return b;
}
