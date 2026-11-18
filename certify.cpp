// certify.cpp
// Generates the states, computes the weights, and outputs the growth rate they certify.
//
//   ./certify P [tolerance=1e-5] [maxIterations=1000]
//
// 1. Generates all states reachable from the start state (generate_states.hpp).
// 2. Iterates the weights until the bounds on the growth rate are within the tolerance (value_iteration.hpp).
//    The bounds are exact fractions with denominator 10^6. The last pass computes T_D exactly for the final
//    weights D, so alpha D(A) <= T_D(A) holds for every state with alpha = the lower bound.
// 3. Checks the remaining condition of the certificate: D > 0 at the start state, and outputs alpha.
//
// Everything is computed with integers. Build: g++ -O3 -march=native -fopenmp -std=c++17 -o certify certify.cpp
#include <sys/resource.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include "value_iteration.hpp"

constexpr int MIN_P = 3, MAX_P = 22;  // the values of P compiled into the program

// Peak memory used by the process so far, in MiB (Linux reports ru_maxrss in KiB; not tested on other OS).
double peakMemoryMiB() {
    rusage usage{};
    return getrusage(RUSAGE_SELF, &usage) == 0 ? usage.ru_maxrss / 1024.0 : -1;
}

// Wall-clock time since construction, in seconds.
class Stopwatch {
    std::chrono::steady_clock::time_point start_ = std::chrono::steady_clock::now();
public:
    double seconds() const { return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count(); }
};

// ---------------------------------------------------------------------------------------------------
template <int P>
int run(uint64_t tolerance, int maxIterations) {
    Stopwatch total;

    // 1. The states.
    Stopwatch generation;
    StateSet<P> Q = generateStates<P>();
    std::cout << "P=" << P << ": " << Q.size() << " states, generated in " << generation.seconds()
              << " s (peak memory " << peakMemoryMiB() << " MiB)" << std::endl;

    // 2. The weights.
    Stopwatch iteration;
    const Bounds b = valueIteration(Q, maxIterations, tolerance);
    std::cout << "value iteration: " << b.iterations << " iterations in " << iteration.seconds()
              << " s, growth rate in [" << decimal(b.lower) << ", " << decimal(b.upper) << "]\n";

    // 3. The certificate, with alpha = the lower bound: the inequalities hold by construction; it remains
    //    to check D > 0 at the start state (id 0).
    size_t positive = 0;
    for (size_t id = 0; id < Q.size(); ++id) positive += (Q.weight(id) > 0);
    std::cout << "D = " << Q.weight(0) << " at the start state, " << positive << " states of positive weight\n"
              << "peak memory " << peakMemoryMiB() << " MiB, total time " << total.seconds() << " s\n";
    const bool valid = Q.weight(0) > 0;
    if (valid)
        std::cout << "certified growth rate: alpha = " << b.lower << "/" << RATIO << " = " << decimal(b.lower) << std::endl;
    else
        std::cout << "no certificate: the start state has weight 0" << std::endl;
    return valid ? 0 : 1;
}

// run<P> for every compiled P, indexed by P - MIN_P.
template <int... I>
constexpr auto runTable(std::integer_sequence<int, I...>) {
    return std::array{&run<MIN_P + I>...};
}

int main(int argc, char** argv) {
    constexpr auto table = runTable(std::make_integer_sequence<int, MAX_P - MIN_P + 1>{});
    try {
        if (argc < 2) throw std::invalid_argument("missing P");
        const int p = std::stoi(argv[1]);
        if (p < MIN_P || p > MAX_P)
            throw std::invalid_argument("P must be between " + std::to_string(MIN_P) + " and " + std::to_string(MAX_P));
        const double tolerance = argc > 2 ? std::stod(argv[2]) : 1e-5;
        const int maxIterations = argc > 3 ? std::stoi(argv[3]) : 1000;
        if (argc > 4) throw std::invalid_argument("too many arguments");
        if (tolerance < 0) throw std::invalid_argument("the tolerance must be nonnegative");
        if (maxIterations < 1) throw std::invalid_argument("maxIterations must be at least 1");

        std::cout << std::fixed << std::setprecision(1);  // timings and memory with one decimal
        return table[p - MIN_P](uint64_t(std::llround(tolerance * RATIO)), maxIterations);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\nusage: " << argv[0]
                  << " P [tolerance=1e-5] [maxIterations=1000]\n";
        return 2;
    }
}
