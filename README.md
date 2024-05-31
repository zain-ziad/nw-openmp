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



![line-graph (2)](https://github.com/zain-ziad/nw-openmp/assets/28985365/3671a97e-9931-482a-8374-f27373a7ce9e)

