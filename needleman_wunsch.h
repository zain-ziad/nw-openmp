// needleman_wunsch.h
#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#include <vector>
#include <string>

struct Result {
    std::string aligned_seq1;
    std::string aligned_seq2;
    std::vector<std::vector<int>> score_matrix;
    int score;
};

Result needleman_wunsch(const std::string &seq1, const std::string &seq2, int gap_penalty, int match_score, int mismatch_penalty);
Result needleman_wunsch_omp(const std::string &seq1, const std::string &seq2, int gap_penalty, int match_score, int mismatch_penalty);

#endif // NEEDLEMAN_WUNSCH_H