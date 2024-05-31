// needleman_wunsch.cpp
#include "needleman_wunsch.h"
#include <vector>
#include <string>
#include <algorithm>

Result needleman_wunsch(const std::string &seq1, const std::string &seq2, int gap_penalty, int match_score, int mismatch_penalty) {
    int m = seq1.size();
    int n = seq2.size();
    std::vector<std::vector<int>> score_matrix(m + 1, std::vector<int>(n + 1));
    
    for (int i = 0; i <= m; ++i)
        score_matrix[i][0] = i * gap_penalty;
    for (int j = 0; j <= n; ++j)
        score_matrix[0][j] = j * gap_penalty;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match = score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_penalty);
            int del = score_matrix[i - 1][j] + gap_penalty;
            int insert = score_matrix[i][j - 1] + gap_penalty;
            score_matrix[i][j] = std::max({match, del, insert});
        }
    }

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

    std::reverse(aligned_seq1.begin(), aligned_seq1.end());
    std::reverse(aligned_seq2.begin(), aligned_seq2.end());

    return {aligned_seq1, aligned_seq2, score_matrix, score_matrix[m][n]};
}