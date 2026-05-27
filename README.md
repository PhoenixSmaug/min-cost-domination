# Minimum Cost Domination Problem

## Problem Description

We want to place chess pieces on an $n \times n$ board such that every cell on the board is either occupied or attacked by a piece (i.e. dominated). The aim is to find a configuration that minimizes the sum of the piece costs, so a queen counts as 9, a rook as 5, a bishop or knight as 3 and a pawn as 1.

Since the pawn attack pattern is not symmetric, it is important to mention that we fix them to always attack upwards. Additionally, as in standard chess the attacks of bishops, rooks and queens are blocked by pieces standing in their path.

## ILP Encoding

Let $q_{i,j}, r_{i,j}, b_{i,j}, k_{i,j}, p_{i,j} \in \{0, 1\}$ be binary decision variables indicating the presence of a Queen, Rook, Bishop, Knight, or Pawn at cell $(i, j)$ respectively. We enforce at most one piece per square:

$$\forall i, j \in \{1, \dots, n\}: \quad q_{i,j} + r_{i,j} + b_{i,j} + k_{i,j} + p_{i,j} \le 1$$

And for convenience we define $o_{i,j}$ as the occupancy of cell $(i,j)$:

$$o_{i,j} = q_{i,j} + r_{i,j} + b_{i,j} + k_{i,j} + p_{i,j}$$

The objective function is:

$$\min \sum_{i=1}^n \sum_{j=1}^n \left( 9 q_{i,j} + 5 r_{i,j} + 3 b_{i,j} + 3 k_{i,j} + p_{i,j} \right)$$

We implement two different approaches to model the line-of-sight propagation and blocking behavior of sliding pieces (Queens, Rooks, and Bishops).

### 1. Full Encoding (Ray Propagation)

To model sliding attacks with blocking without introducing exponential numbers of constraints, we define directional propagation variables. For a sliding piece type $x \in \{q, r, b\}$ and direction vector $d \in D_x$:

*   Let $a_{i,j,d} \in \{0, 1\}$ represent an incoming attack ray entering cell $(i, j)$ along direction $d$.
*   Let $t_{i,j,d} \in \{0, 1\}$ represent an outgoing attack ray exiting cell $(i, j)$ along direction $d$.

The incoming ray is connected to the adjacent cell's outgoing ray:

$$a_{i,j,d} = \begin{cases} t_{i-d_i, j-d_j, d} & \text{if } (i-d_i, j-d_j) \text{ is within bounds} \\ 0 & \text{otherwise} \end{cases}$$

The propagation logic is linearized using the following inequalities:

$$t_{i,j,d} \le a_{i,j,d} + x_{i,j}$$
$$t_{i,j,d} \ge x_{i,j}$$
$$t_{i,j,d} \le 1 - o_{i,j} + x_{i,j}$$
$$t_{i,j,d} \ge a_{i,j,d} - o_{i,j}$$

This ensures that:
*   An outgoing ray is generated if a source piece of type $x$ is placed on $(i, j)$ ($x_{i,j} = 1$).
*   The ray propagates freely through empty cells ($o_{i,j} = 0 \implies t_{i,j,d} = a_{i,j,d}$).
*   The ray is terminated if it hits any other piece ($o_{i,j} = 1 \land x_{i,j} = 0 \implies t_{i,j,d} = 0$).

Domination of each cell $(i, j)$ is then enforced by:

$$o_{i,j} + \text{Attacks}_{\text{sliding}}(i, j) + \text{Attacks}_{\text{knight}}(i, j) + \text{Attacks}_{\text{pawn}}(i, j) \ge 1$$

where $\text{Attacks}\_{\text{sliding}}(i, j)$ is the sum of incoming rays $a_{i,j,d}$ for all valid sliding piece types and directions.

### 2. Lazy Constraint Encoding

In practice, modeling the full physical ray propagation network requires a large number of binary variables and constraints. 

Alternatively, we can solve a relaxed model where blocking is initially ignored. Let $D_{i,j}$ represent the set of cells from which a sliding piece could attack cell $(i,j)$ if the board were otherwise empty. We enforce the initial relaxed domination constraints:

$$\forall i,j: \quad o_{i,j} + \sum_{(r, c) \in D_{i,j}} (\dots) \ge 1$$

Whenever the ILP solver finds an integer candidate solution, we verify the configuration using an early-exit board validation sweep. If any cell is found to be undefended due to blocking, we exclude this specific solution by adding the following constraint:

$$\sum_{(r, c) \in \mathcal{S}_{\text{active}}} (1 - v_{r,c}) + \sum_{(r, c) \in \mathcal{S}_{\text{inactive}}} v_{r,c} \ge 1$$

where $v_{r,c}$ represents the active and inactive piece variables of the candidate solution.

This process in then simply repeated, until even with blocking all tiles are attacked. Note that the solver is already incentivised to avoid pieces blocking each other, since the more tiles a piece attack, the less pieces we need and thus save cost. So in practice the process terminates very quickly, making the lazy approach much more efficient for larger $n$.

Both ILP approaches are implemented in `main.jl` using the commerical solver CPLEX (the open-source HiGHS solver currently has no support for lazy callbacks). Additionally, we provide the optimal solutions found by the lazy ILP solver in 24 hours of runtime in `data/sol.txt`.

## An Explicit Solution

Based on the data of the previous solver, I found the following explicit solution for all $n$ (here $c := \lfloor n/2 \rfloor$):

### Even $n \ge 8$
We place $n-2$ Bishops along column $c$:
$$\text{Bishops at } (r, c) \quad \forall r \in \{2, \dots, n-1\}$$
And place $4$ Pawns at:
$$(1, c), \quad (n, c), \quad (n + 1 - c, n), \quad (n - c, n)$$

This configuration yields a total cost of:

$$\text{Cost}_{\text{even}} = 3(n-2) + 4 = 3n - 2$$

### Odd $n \ge 7$
We place $n-2$ Bishops distributed across two columns:
$$\text{Bishops at } (r, c+1) \quad \forall r \in \{2, 4, \dots, n-1\}$$
$$\text{Bishops at } (r, c) \quad \forall r \in \{2, 4, \dots, n-3\}$$
And place $3$ Pawns at:
$$(n + 1 - c, n-1), \quad (2, c), \quad (n + 2 - c, n)$$

This configuration yields a total cost of:

$$\text{Cost}_{\text{odd}} = 3(n-2) + 3 = 3n - 3$$

These solution gives us upper bounds for the optimal value for all $n$. In fact, using the lazy ILP encoding one can verify that there are no better solution for any $n$ up to $n = 23$. I conjecture that this still holds true for all larger $n$, but it seems very hard to rigorously prove. As a reference, if we restrict ourselves to only placing queens, the problem turns into the Minimum Dominating Queen Problem which was first posed in 1862 by Carl Friedrich de Jaenisch. As explained [here](https://arxiv.org/abs/1606.02060), to this day no one could proof the seemingly obvious fact that the minimum number of queens grows monotonically with $n$, much less give an explicit formula for the objective.

## Acknowledgements

The idea to investigate the Minimum Cost Domination Problem came from [Bernardo Subercaseaux](https://bsubercaseaux.github.io/) in discussions following the SAT 2025 conference.

(c) Mia Müßig

