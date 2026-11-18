# Avoiding-squares-over-lists-of-size-3
This program carries out the computer-assisted part of the proof in

> Matthieu Rosenfeld, *[title]*, preprint, [link to come].

The paper proves that every sequence of lists of size three admits an infinite square-free word, with at least (623/500)^n such words of length n. The proof reduces the problem to the existence of a *weight certificate*: a function D on a finite set of abstract states such that D(A₀) > 0 and T_D(A) ≥ α D(A) for every state A, with α large enough. This program constructs such a certificate for p = 21 and outputs the constant α it certifies.


## Compilation

A C++17 compiler with OpenMP is needed, for instance GCC:

```
g++ -O3 -march=native -fopenmp -std=c++17 -o certify certify.cpp
```

Without `-fopenmp`, the program compiles but runs on a single core. On macOS, Apple's clang does not
support `-fopenmp`; use GCC from Homebrew, or clang with `libomp`.

## Usage

```
./certify P [tolerance=1e-5] [maxIterations=1000]
```

* `P` is the largest period of the squares handled by the states (3 ≤ P ≤ 22). The paper uses P = 21.
* The iteration stops when the upper and lower bounds on the growth rate differ by at most `tolerance`,
  or after `maxIterations` iterations.
* The number of threads is set with `OMP_NUM_THREADS`; by default, all cores are used. The output does
  not depend on the number of threads, apart from the timings: every computation is exact and uses
  integers only.

For a quick test, `./certify 12` takes about ten seconds on one core. For the result of the paper:

```
OMP_NUM_THREADS=$(nproc) ./certify 21
```

This run needs about 21 GB of memory, and took 35 hours on [machine, number of cores]. On a machine
with several processors, `numactl --interleave=all ./certify 21` spreads the memory over all of them.

## Output

```
P=12: 59983 states, generated in 0.1 s (peak memory 7.0 MiB)
  iteration 10: growth rate in [1.000000, 1.666667]
  ...
value iteration: 108 iterations in 8.9 s, growth rate in [1.306550, 1.306560]
D = 91902900 at the start state, 57190 states of positive weight
peak memory 7.0 MiB, total time 8.9 s
certified growth rate: alpha = 1306550/1000000 = 1.306550
```

* The first line gives the number of states reachable from the start state (the set R of the paper).
* Every ten iterations, the program prints the current bounds on the growth rate.
* At the end, the weights D are integers between 0 and 10^8. The program computes T_D exactly, and the
  lower bound `lower` satisfies 10^6 T_D(A) ≥ `lower` · D(A) for every state A. Together with
  D(A₀) > 0, this makes D a weight certificate with constant α = `lower`/10^6.
* The program also prints the number of states of positive weight, the peak memory and the total time,
  which are the quantities reported in the paper.
* The exit code is 0 when a certificate is found, 1 when D(A₀) = 0, and 2 for invalid arguments.

To conclude, α must satisfy the criterion of the paper, (α − β)(β − 1)β^(p−1) ≥ 1 for some β > 1.
For p = 21 and β = 623/500, this holds as soon as α ≥ 1.29597. This last check is done in the paper,
not by the program.

The smaller values of P are useful for testing, but their constants are too small for the criterion;
p = 21 is the smallest value that works. Reference values, on one core:

| P | States | Iterations | α | Time |
|---|---|---|---|---|
| 10 | 7 946 | 80 | 1.317815 | 0.5 s |
| 12 | 59 983 | 108 | 1.306550 | 9 s |
| 14 | 439 401 | 100 | 1.303265 | 2 min |
| 21 | 405 488 434 | 109 | [value of the run] | 35 h on [cores] |



## Files

| File | Contents |
|---|---|
| `state.hpp` | The abstract states, the start state and the successor of a state |
| `generate_states.hpp` | The outside types, the hash table of states, and the search for the reachable states |
| `value_iteration.hpp` | The operator T_D and the iteration computing the weights |
| `certify.cpp` | The main program |

### Correspondence with the notation of the paper:

| Paper | Code |
|---|---|
| position j (distance from the end of the word) | lag j |
| labels 0, 1, ★ | `Zero`, `One`, `Joker` |
| r_q | `run(q)` |
| start state A₀ (state of the prefix ρ) | `State<P>::start()` |
| Φ(A, X), allowed types | `append(X, next)`, which returns `false` when X is not allowed |
| 𝒳(A) | `outsideTypes(s, types)` |
| v_D(A, X) | `extensionWeight(s, X)` |
| T_D(A) | `cheapestList(Q, s)` |
