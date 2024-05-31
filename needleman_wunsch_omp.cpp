// needleman_wunsch_omp.cpp
#include "needleman_wunsch.h"
#include <vector>
#include <string>
#include <algorithm>
#include <omp.h>

Result needleman_wunsch_omp(const std::string &seq1, const std::string &seq2, int gap_penalty, int match_score, int mismatch_penalty) {
    int m = seq1.size();
    int n = seq2.size();
    std::vector<std::vector<int>> score_matrix(m + 1, std::vector<int>(n + 1));
    
    // Initialize the scoring matrix
    for (int i = 0; i <= m; ++i)
        score_matrix[i][0] = i * gap_penalty;
    for (int j = 0; j <= n; ++j)
        score_matrix[0][j] = j * gap_penalty;

    // Fill the scoring matrix using OpenMP with diagonal wavefront parallelization
    for (int k = 2; k <= m + n; ++k) {
        #pragma omp parallel for
        for (int i = std::max(1, k - n); i <= std::min(m, k - 1); ++i) {
            int j = k - i;
            int match = score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_penalty);
            int del = score_matrix[i - 1][j] + gap_penalty;
            int insert = score_matrix[i][j - 1] + gap_penalty;
            score_matrix[i][j] = std::max({match, del, insert});
        }
    }

    // Traceback to get the aligned sequences
    int i = m, j = n;
    std::string aligned_seq1, aligned_seq2;
    while (i > 0 && j > 0) {
        if (score_matrix[i][j] == score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_penalty)) {
            aligned_seq1 += seq1[i - 1];
            aligned_seq2 += seq2[j - 1];
            --i; --j;
        } else if (score_matrix[i][j] == score_matrix[i - 1][j] + gap_penalty) {
            aligned_seq1 += seq1[i - 1];
            aligned_seq2 += '-';
            --i;
        } else {
            aligned_seq1 += '-';
            aligned_seq2 += seq2[j - 1];
            --j;
        }
    }
    
    // Append remaining sequence if any
    while (i > 0) {
        aligned_seq1 += seq1[i - 1];
        aligned_seq2 += '-';
        --i;
    }
    while (j > 0) {
        aligned_seq1 += '-';
        aligned_seq2 += seq2[j - 1];
        --j;
    }

    // Reverse the aligned sequences to get the correct order
    std::reverse(aligned_seq1.begin(), aligned_seq1.end());
    std::reverse(aligned_seq2.begin(), aligned_seq2.end());

    return {aligned_seq1, aligned_seq2, score_matrix, score_matrix[m][n]};
}