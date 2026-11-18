// generate_states.hpp
// All abstract states reachable from the start state (the set R of the paper).
// Every word starts with a prefix of P + 1 new letters, whose state is State<P>::start().
// A letter appended to a state is described by the set X of lags where it already occurs,
// and we use State<P>::append(X, next) to compute the successor.
#pragma once
#include <cstdint>
#include <utility>
#include <vector>

#include "state.hpp"

// --------------------------------------------------------------------
// Outside types: the family Types(A) of the paper. 
// A type X is a union of classes of joker lags, possibly empty,
// such that the symbolic word s_X is square-free and X does not contain
// both lags r_q + 1 and r_q + 1 + q for any q (maximality of the runs).
// --------------------------------------------------------------------
template <int P>
void outsideTypes(const State<P>& s, std::vector<uint32_t>& types) {
    const uint32_t jokers = s.lagsWith(State<P>::Joker);
    auto isJoker = [&](int q) { return (jokers >> q) & 1u; };
 
    // The classes of Joker lags: 
    // for every q and every i <= r_q, the letters at lags i and i + q are equal.
    // classOf[q] is the representative of the class of lag q.
    int classOf[P + 1];
    for (int q = 0; q <= P; ++q) classOf[q] = q;
    for (int q = 1; q <= P; ++q)
        for (int i = 1; i <= s.run(q) && i + q <= P; ++i) {
            if (!isJoker(i)) continue;
            const int renamed = classOf[i + q];
            for (int j = 1; j <= P; ++j)
                if (classOf[j] == renamed) classOf[j] = classOf[i];
        }

    // classAt[c] contains the class whose representative is c, as a bit mask (bit q <=> lag q), or 0.
    uint32_t classAt[P + 1] = {};      
    for (int q = 1; q <= P; ++q)
        if (isJoker(q)) classAt[classOf[q]] |= 1u << q;
 
    uint32_t classes[P];  // the classes of Joker lags
    int countClasses = 0;
    for (int c = 1; c <= P; ++c)
        if (classAt[c]) classes[countClasses++] = classAt[c];
  
    // Maximality of the runs: the run of period q stops at lag r_q + 1, so the letters at lags r_q + 1
    // and r_q + 1 + q are different, and no outside letter can occupy both of their classes.
    // conflict[j] contains the lags of the classes that cannot share a letter with classes[j].
    uint32_t conflict[P] = {};
    for (int q = 1; q <= P; ++q) {
        const int a = s.run(q) + 1, b = a + q;
        if (b > P || !isJoker(a) || !isJoker(b)) continue;
        int ja = 0, jb = 0;  // the indices of the classes of a and b
        for (int j = 0; j < countClasses; ++j) {
            if (classes[j] >> a & 1u) ja = j;
            if (classes[j] >> b & 1u) jb = j;
        }
        if (ja != jb) {
            conflict[ja] |= classes[jb];
            conflict[jb] |= classes[ja];
        }
    }

    // Two lags hold the same symbol if both are labelled Zero, both One, or both are Joker lags of the
    // same class. Bit t of equal[r] is set if lags t and t + r hold the same symbol, so the lags
    // t, ..., t + 2r - 1 hold a square if and only if the bits t, ..., t + r - 1 of equal[r] are all set.
    const uint32_t zero = s.lagsWith(State<P>::Zero), one = s.lagsWith(State<P>::One);
    uint32_t equal[P / 2 + 1];
    for (int r = 1; 2 * r <= P; ++r) {
        equal[r] = (zero & zero >> r) | (one & one >> r);
        for (int j = 0; j < countClasses; ++j) equal[r] |= classes[j] & classes[j] >> r;
    }

    // This function verifies if the symbolic word s_Y, with the new letter on the lags of Y, contains a square
    auto hasSquare = [&](uint32_t Y) {
        for (int r = 1; 2 * r <= P; ++r) {
            uint32_t E = equal[r] | (Y & Y >> r);  // with the new letter on Y: two lags of Y become equal
            // next a small trick to check if k consecutive bits of E are all 1
            // there is a faster version that double the shift but it would be interesting for values of k larger than 20...
            for (int k = 1; k < r; ++k) E &= E >> 1;
            if (E) return true;
        }
        return false;
    };
 
    // The types are constructed as unions of classes, starting from the empty union.
    // A union that has a square or a conflict is rejected with all the unions that contain it: 
    // (adding classes can only create more squares and more conflicts.)
    types.clear();
    // If s_0 already has a square, the family is empty. This never happens at a reachable state, but
    // keeping this test means that we do not need to prove it.
    if (hasSquare(0)) return;  
    // if every s_X has a square: the family is empty
    // this never happens for a reachable state and this could be proven
    // but keeping this check means we do not have to prove it
    types.push_back(0);
    static thread_local std::vector<std::pair<int, uint32_t>> stack;
    stack.assign(1, {0, 0u});
    while (!stack.empty()) {
        const auto [first, X] = stack.back();
        stack.pop_back();
        for (int j = first; j < countClasses; ++j) {
            if (X & conflict[j]) continue;  // classes[j] conflicts with a class of X
            const uint32_t Y = X | classes[j];
            if (hasSquare(Y)) continue;
            types.push_back(Y);
            stack.push_back({j + 1, Y});
        }
    }
}

// ---------------------------------------------------------------------------------------------------
// StateSet: a hash table of states, with open addressing (linear probing).
//   * Each entry holds a state, its current weight x and its next weight y,
//     so one lookup finds a state and its weights. 
//   * An entry takes 32 bytes: 24 for the state and 4 for each weight. 
//       And we need 4 more bytes for the index of each state 
//       (which might cost between 8 and 16 bytes per state because 
//        we need the array at least twice as large as the number of states).
//   * Integer weights (bounded by SCALE = 10^8) keep every computation exact: 
//     a sum of three weights stays below 2^32.  
//   * Entries are kept in a list, in order of insertion: the position of an entry is the id of its state.
//     The open-addressing table itself is index_: it stores ids and is kept at most half full.
//   * Writing our own table, instead of using eg std::unordered_set allow all these small optimizations.
// ---------------------------------------------------------------------------------------------------
template <int P>
class StateSet {
    struct Entry {
        State<P> state;
        uint32_t weight = SCALE;  // current weight; the value iteration starts from SCALE
        uint32_t nextWeight = 0;  // next weight, used in the iteration
    };
    std::vector<Entry> entries_;
    std::vector<uint32_t> index_ = std::vector<uint32_t>(1024, NONE);
    int bits_ = 10;

    // Multiplicative hashing: the top bits of hash() * (2^64 / golden ratio) select one of 2^bits_ slots.
    size_t home(const State<P>& s) const { return size_t((s.hash() * 0x9E3779B97F4A7C15ull) >> (64 - bits_)); }
    size_t next(size_t i) const { return (i + 1) & (index_.size() - 1); }

    void grow() {
        ++bits_;
        index_.assign(size_t(1) << bits_, NONE);
        for (uint32_t id = 0; id < entries_.size(); ++id) {
            size_t i = home(entries_[id].state);
            while (index_[i] != NONE) i = next(i);
            index_[i] = id;
        }
    }
public:
    static constexpr uint32_t NONE = ~0u;
    static constexpr uint32_t SCALE = 100000000;  // largest weight, 10^8

    // Adds s if it is new; returns its id.
    uint32_t insert(const State<P>& s) {
        if (2 * (entries_.size() + 1) > index_.size()) grow();
        size_t i = home(s);
        for (; index_[i] != NONE; i = next(i))
            if (entries_[index_[i]].state == s) return index_[i];
        index_[i] = uint32_t(entries_.size());
        entries_.push_back({s});
        return index_[i];
    }

    // The id of s, or NONE.
    uint32_t find(const State<P>& s) const {
        for (size_t i = home(s); index_[i] != NONE; i = next(i))
            if (entries_[index_[i]].state == s) return index_[i];
        return NONE;
    }

    size_t size() const { return entries_.size(); }
    const State<P>& operator[](size_t id) const { return entries_[id].state; }

    // Weights, by id or by state (the state must be in the set).
    uint32_t weight(size_t id) const { return entries_[id].weight; }
    void setWeight(size_t id, uint32_t x) { entries_[id].weight = x; }
    uint32_t weightOf(const State<P>& s) const { return entries_[find(s)].weight; }

    uint32_t nextWeight(size_t id) const { return entries_[id].nextWeight; }
    void setNextWeight(size_t id, uint32_t y) { entries_[id].nextWeight = y; }

    // v_D(A, X) of the paper: the weight of the successor of s by a letter of type X, or 0 if that letter
    // completes a square. The successor must be in the set.
    uint32_t extensionWeight(const State<P>& s, uint32_t X) const {
        State<P> successor;
        return s.append(X, successor) ? weightOf(successor) : 0;
    }    
};

// Breadth-first search from the start state. 
template <int P>
StateSet<P> generateStates() {
    StateSet<P> Q;
    Q.insert(State<P>::start());
    std::vector<uint32_t> outside;
    for (size_t id = 0; id < Q.size(); ++id) {
        const State<P> s = Q[id];
        State<P> next;
        if (s.append(s.lagsWith(State<P>::One), next)) Q.insert(next);
        outsideTypes(s, outside);
        for (uint32_t X : outside)
            if (s.append(X, next)) Q.insert(next);
    }
    return Q;
}
