// state.hpp
// The abstract states of the paper (Section "Finite states and deletion of long squares").
#pragma once
#include <cstdint>

// A state describes the last P positions of a word. 
// In the code, the position j of the paper (the letter
// at distance j from the end of the word) is called lag j: lag 1 is the last letter.
// A state holds:
//   * a label per lag: the last letter (Zero, label 0 of the paper), the letter before it (One, label 1),
//     or any other letter (Joker, label * of the paper);
//   * a run length per period q: r_q, the number of trailing letters equal to the letter q positions
//     earlier. For words without squares of period at most P, r_q <= q - 1, and r_1 = 0.
//
// Layout, 24 bytes:
//   labels_   2 bits per lag q = 1..P,        bits 2(q-1) .. 2(q-1)+1
//   runs_[0]  5 bits per period q = 1..12,    bits 5(q-1) .. 5(q-1)+4
//   runs_[1]  5 bits per period q = 13..24,   bits 5(q-13) .. 5(q-13)+4
// A run field never uses two memory words (12 fields x 5 bits = 60 bits per word).

template <int P>
class State {
    static_assert(2 <= P && P <= 24, "The period P must be in [2, 24].");
    uint64_t labels_ = 0;
    uint64_t runs_[2] = {0, 0};
public:
    enum Label : uint8_t { Zero = 0, One = 1, Joker = 2 };
    
    // The start state A_0 of the paper: the state of the prefix rho = c_P ... c_1 c_0 of P + 1 new letters
    // that starts every word. Lag 1 holds the last letter, lag 2 the letter before it, the other lags hold
    // a distinct letter each, and no run has started.
    static State start() {
        State s;
        for (int q = 2; q <= P; ++q) s.setLabel(q, q == 2 ? One : Joker);
        s.setLabel(1, Zero);
        return s;  // every run is 0
    }

    // Reads and writes the label at lag q.
    Label label(int q) const { return Label((labels_ >> (2 * (q - 1))) & 3u); }
    void setLabel(int q, Label v) {
        labels_ = (labels_ & ~(3ull << (2 * (q - 1)))) | (uint64_t(v) << (2 * (q - 1)));
    }

    // Reads and writes the run length r_q. (q - 1) / 12 selects the word, (q - 1) % 12 the field in that word.
    int run(int q) const { return int((runs_[(q - 1) / 12] >> (5 * ((q - 1) % 12))) & 31u); }
    void setRun(int q, int r) {
        uint64_t& w = runs_[(q - 1) / 12];
        w = (w & ~(31ull << (5 * ((q - 1) % 12)))) | (uint64_t(r) << (5 * ((q - 1) % 12)));
    }

    // The lags whose label is v, as a bit mask (bit q <=> lag q).
    uint32_t lagsWith(Label v) const {
        uint32_t lags = 0;
        for (int q = 1; q <= P; ++q)
            if (label(q) == v) lags |= 1u << q;
        return lags;
    }

    // ---- the successor Phi(A, X) of the paper ----
    // The appended letter is described by the set X of lags (bit q <=> lag q) where it occurs:
    //   the letter before the last one: X = lagsWith(One);
    //   an outside letter: X = an element of outsideTypes, a set of Joker lags (0 for a new letter).
    // The last letter, X = lagsWith(Zero), is never allowed: it would create the square aa.
    // Returns false if X is not allowed, that is, if the letter completes a square of period q <= P: 
    // it repeats the letter at lag q and r_q = q - 1.
    // Otherwise writes the successor to next and returns true.
    bool append(uint32_t X, State& next) const {
        State result;  // every run is 0 until set below
        // Labels: every lag moves one step further (lag q becomes lag q + 1; lag P drops out) and is
        // relabelled Zero -> One, One -> Joker, Joker -> Joker. On the two bits (hi, lo) of a label, this is
        // hi' = hi | lo and lo' = (hi == lo). Lag 1, the appended letter, gets Zero (bits 00) from the shift.
        constexpr uint64_t LOW = 0x5555555555555555ull;  // the low bit of every label
        const uint64_t lo = labels_ & LOW, hi = (labels_ >> 1) & LOW;
        result.labels_ = ((((hi | lo) << 1) | (~(hi ^ lo) & LOW)) << 2) & ((1ull << (2 * P)) - 1);
        // The lags q of X hold the appended letter: lag q + 1 becomes Zero, and the run of period q continues.
        for (uint32_t lags = X; lags; lags &= lags - 1) {  // removes the lowest bit of lags at each step
            const int q = __builtin_ctz(lags);             // the lowest lag of X not treated yet
            if (run(q) == q - 1) return false;             // the run would complete a square of period q
            result.setRun(q, run(q) + 1);
            if (q < P) result.setLabel(q + 1, Zero);
        }
        next = result;
        return true;
    }

    // ---- comparison, hashing ----
    bool operator==(const State& o) const {
        return labels_ == o.labels_ && runs_[0] == o.runs_[0] && runs_[1] == o.runs_[1];
    }
    bool operator!=(const State& o) const { return !(*this == o); }

    uint64_t hash() const {
        constexpr uint64_t K = 0x9E3779B97F4A7C15ull;  // ~ 2^64 / golden ratio, odd
        return (labels_ * K + runs_[0]) * K + runs_[1];
    }
};
