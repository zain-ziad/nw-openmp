# Serial & Parallel Implementation of the Needleman-Wunsch Global Sequence Alignment Algorithm in C++
In the field of bioinformatics, sequence alignment is a fundamental task for comparing and analyzing biological sequences. The Needleman-Wunsch algorithm is a widely used dynamic programming approach for finding the optimal global alignment between two sequences. However, the computational complexity of the algorithm poses challenges when dealing with long sequences. In this project, the implementation of both serial and parallel versions of the Needleman-Wunsch algorithm in C++ using OpenMP was done. I've compared their performance across various sequence lengths to demonstrate the effectiveness of parallelization.

## Needleman-Wunsch Algorithm
The Needleman-Wunsch algorithm operates on a scoring matrix to compute the optimal alignment between two sequences. The algorithm can be broken down into three main steps:

### 1. Initialization:
The first row and column of the scoring matrix are initialized with gap penalties. For a matrix $M$ of size $(m+1)*(n+1)$, where $m$ and $n$ are the lengths of the two sequences:
```math
M[i][0] = i \times \text{gap\_penalty} \quad \text{for} \quad 0 \leq i \leq m
```
```math
M[0][j] = j \times \text{gap\_penalty} \quad \text{for} \quad 0 \leq j \leq n
```
### 2. Matrix Filling:
Each cell in the matrix is filled using the recurrence relation:
```math
M[i][j] = \max \left\{
\begin{array}{l}
M[i-1][j-1] + \text{score}(seq1[i-1], seq2[j-1]),\ (match/mismatch)\\
M[i-1][j] + \text{gap\_penalty},\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ (deletion)\\
M[i][j-1] + \text{gap\_penalty}\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  \ \ \ \ \ \ \ (insertion)
\end{array}
\right.
```
This ensures that the optimal score is computed based on previous cells, capturing the best alignment up to the current position.
### 3. Traceback:
The optimal alignment is reconstructed by tracing back from the bottom-right cell of the matrix to the top-left cell, following the path of optimal scores.

## Computational Complexity
* ### Serial Complexity:
  The serial implementation of the Needleman-Wunsch algorithm has a time and space complexity of $O(m*n)$, where $m$ and $n$ are the lengths of the two sequences.
* ### Parallel Complexity:
  The parallel implementation using OpenMP divides the matrix computation into independent tasks that can be executed concurrently. By parallelizing along the anti-diagonals (wavefront parallelism), the time complexity remains $O(m*n)$, but the wall-clock time is reduced due to concurrent execution. The efficiency of the parallel implementation depends on the number of available processors $P$ and the overhead of managing parallel threads. Ideally, the execution time is reduced by a factor proportional to $P$, yielding a parallel time complexity of $O((m×n)/P)$.

## Experimental Setup
The performance of both implementations was evaluated on sequences with lengths of 50, 500, 1,000, 5,000, 10,000, 15,000, 20,000, 25,000, and 30,000. The experiments were conducted on a Ryzen 7 5700X, (8 Cores, 16 Threads @4.1GHz) with 16GB of DDR4 memory.  

## Results
<div align="center">
  
| Sequence Length | Serial Time (s) | OpenMP Time (s) |
|-----------------|-----------------|-----------------|
| 50              | 0.0001025       | 0.00009325      |
| 500             | 0.0075027       | 0.00053207      |
| 1000            | 0.0293138       | 0.0021382       |
| 5000            | 0.735076        | 0.073496        |
| 10000           | 2.91304         | 0.28411         |
| 15000           | 6.57172         | 0.65023         |
| 20000           | 11.6548         | 1.16298         |
| 25000           | 18.2623         | 1.82013         |
| 30000           | 26.3945         | 2.63911         |

![line-graph (2)](https://github.com/zain-ziad/nw-openmp/assets/28985365/3671a97e-9931-482a-8374-f27373a7ce9e)

</div>

## Discussion:
The parallel implementation of the Needleman-Wunsch algorithm using OpenMP demonstrates the potential for accelerating sequence alignment tasks in bioinformatics. By leveraging parallel computing resources, the execution time can be significantly reduced, especially for long sequences. This optimization is particularly valuable in scenarios where large-scale sequence comparisons are required, such as in genome-wide studies or comparative genomics.

The speedup achieved by the parallel implementation is influenced by various factors, including the number of available cores, the sequence lengths, and the overhead associated with parallel execution. While the parallel implementation outperforms the serial implementation for longer sequences, we should note that for very short sequences, the overhead of parallel execution may outweigh the benefits.

The parallel implementation using OpenMP is limited to shared-memory systems. For even larger-scale sequence alignment tasks, distributed computing frameworks such as MPI (Message Passing Interface) can be considered to harness the power of multiple machines or clusters.

## How to Use:
To compile the code with OpenMP support, use the following command:
```
g++ main.cpp needleman_wunsch.cpp needleman_wunsch_omp.cpp -fopenmp -o needleman_wunsch
```
Run the compiled executable:

```
./needleman_wunsch
```
Note: Sequences must be stored as seq1.txt and seq2.txt. The output would be three files, one contains the execution time for both serial and parallel implementation, the other two files are the resulting score and sequences. 
