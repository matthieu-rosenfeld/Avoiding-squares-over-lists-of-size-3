# Avoiding-squares-over-lists-of-size-3
This program carries out the computer-assisted part of the proof in

> Matthieu Rosenfeld, *Avoiding squares over words with lists of size three*, preprint, [link to come].

The paper proves that every sequence of lists of size three admits an infinite square-free word, with at least $c\,(5/4)^n$ such words of length $n$, for some constant $c>0$. The proof reduces the problem to the existence of a *weight certificate*: a function $D$ on a finite set of abstract states such that $D(A_0)>0$ and $T_D(A)\ge\alpha D(A)$ for every state $A$, with $\alpha$ large enough. The states handle the squares of period at most $p$, and the counting argument of the paper handles the longer ones. This program constructs such a certificate for $p = 21$ and outputs the constant $\alpha$ it certifies.


## Compilation

A C++17 compiler with OpenMP is needed, for instance GCC:

```
g++ -O3 -march=native -fopenmp -std=c++17 -o certify certify.cpp
```

Without `-fopenmp`, the program compiles but runs on a single core. 
The program was written for Linux and might not work on other operating systems.

## Usage

```
./certify P [tolerance=1e-5] [maxIterations=1000]
```

* `P` is the value of $p$, the largest period of the squares handled by the states ($3\le p\le 22$). The paper uses $p=21$.
* The iteration stops when the upper and lower bounds on the growth rate differ by at most `tolerance`,
  or after `maxIterations` iterations.
* The number of threads is set with `OMP_NUM_THREADS`; by default, all cores are used. The output does
  not depend on the number of threads, apart from the timings: every computation is exact and uses
  integers only.

For a quick test, `./certify 12` takes about ten seconds on one core. For the result of the paper:

```
OMP_NUM_THREADS=$(nproc) ./certify 21
```

This run needs about 20 GiB of memory, and took about 18 hours with 40 threads (see the section *Result for $p=21$* below).

## Output

```
P=12: 59983 states, generated in 0.1 s (peak memory 6.8 MiB)
  iteration 10: growth rate in [1.000000, 1.666667]
  [...]
  iteration 100: growth rate in [1.306543, 1.306565]
value iteration: 108 iterations in 7.0 s, growth rate in [1.306550, 1.306560]
D = 91902900 at the start state, 57190 states of positive weight
peak memory 7.0 MiB, total time 7.1 s
certified growth rate: alpha = 1306550/1000000 = 1.306550

```

* The first line gives the number of states reachable from the start state (the set $R$ of the paper).
* Every ten iterations, the program prints the current bounds on the growth rate.
* At the end, the weights $D(A)$ are integers between $0$ and $10^8$. The program computes $T_D$ exactly, and the lower bound $\ell$ satisfies $10^6 T_D(A)\ge \ell D(A)$ for every state $A$. Together with $D(A_0)>0$, this makes $D$ a weight certificate with constant $\alpha = \ell/10^6$.
* The program also prints the number of states of positive weight, the peak memory and the total time,
  which are the quantities reported in the paper.
* The exit code is `0` when a certificate is found, `1` when $D(A_0)=0$, and `2` for invalid arguments.

To conclude, $\alpha$ must satisfy the criterion of the paper, $(\alpha-\beta)(\beta-1)\beta^{p-1}\ge1$ for some $\beta>1$.
For $p = 21$ and $\beta = 5/4$, this holds as soon as $\alpha \ge 1.29612$. This last check is done in the paper,
not by the program.

The smaller values of $p$ are useful for testing. Reference values, on one core:

| $p$ | States | Iterations | $\alpha$ | Time |
|---|---|---|---|---|
| 10 | 7 946 | 80 | 1.317815 | 0.5 s |
| 12 | 59 983 | 108 | 1.306550 | 9 s |
| 14 | 439 401 | 100 | 1.303265 | 2 min |

For $p=21$, see the section *Result for $p=21$* below.

## Files

| File | Contents |
|---|---|
| `state.hpp` | The abstract states, the start state and the successor of a state |
| `generate_states.hpp` | The outside types, the hash table of states, and the search for the reachable states |
| `value_iteration.hpp` | The operator $T_D$ and the iteration computing the weights |
| `certify.cpp` | The main program |

### Correspondence with the notation of the paper

| Paper | Code |
|---|---|
| position $j$ (the $j$-th letter from the end of the word) | lag $j$ |
| labels $0$, $1$, $\star$ | `Zero`, `One`, `Joker` |
| $r_q$ | `run(q)` |
| start state $A_0$ (state of the prefix $\rho$) | `State<P>::start()` |
| $\Phi(A,X)$, allowed types | `append(X, next)`, which returns `false` when $X$ is not allowed |
| $\mathcal{X}_{\mathrm{out}}(A)$ | `outsideTypes(s, types)` |
| $v_D(A,X)$ | `extensionWeight(s, X)` |
| $T_D(A)$ | `cheapestList(Q, s)` |

## Result for $p=21$

The certificate used in the paper was computed with

```
OMP_NUM_THREADS=40 ./certify 21 1e-5 > run_21_1e-5.out
```

using 40 threads on a shared virtual machine (two AMD EPYC 9554 processors, 691 GiB of RAM, Debian 13, GCC 14.2), so the timings are only indicative. The complete output, also available in `run_21_1e-5.out`, is:

```
P=21: 405488434 states, generated in 16855.8 s (peak memory 20483.8 MiB)
  iteration 10: growth rate in [0.666666, 2.000001]
  iteration 20: growth rate in [1.173913, 1.450000]
  iteration 30: growth rate in [1.253333, 1.333334]
  iteration 40: growth rate in [1.290539, 1.311788]
  iteration 50: growth rate in [1.297129, 1.303966]
  iteration 60: growth rate in [1.299525, 1.301401]
  iteration 70: growth rate in [1.300253, 1.300844]
  iteration 80: growth rate in [1.300449, 1.300643]
  iteration 90: growth rate in [1.300513, 1.300586]
  iteration 100: growth rate in [1.300537, 1.300560]
value iteration: 109 iterations in 49366.9 s, growth rate in [1.300544, 1.300553]
D = 89182122 at the start state, 394977967 states of positive weight
peak memory 20483.8 MiB, total time 66223.1 s
certified growth rate: alpha = 1300544/1000000 = 1.300544
```

| Quantity | Value |
|---|---|
| Reachable states | 405 488 434 |
| States of positive weight | 394 977 967 |
| Weight of the start state | $D(A_0) = 89\,182\,122 > 0$ |
| Iterations | 109 |
| Certified constant | $\alpha = 1300544/10^6 = 1.300544$ |
| Peak memory | 20 484 MiB (20.0 GiB) |
| Time | 4 h 41 min to generate the states (one thread), 13 h 43 min for the iteration (40 threads), 18 h 24 min in total |

With $\beta = 5/4$, the criterion of the paper holds: $(\alpha-\beta)(\beta-1)\beta^{20} \approx 1.096 \ge 1$.
In other words, with $\beta = 5/4$, any certified constant $\alpha \ge 1.29612$ would suffice.
