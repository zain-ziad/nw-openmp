// main.cpp
#include "needleman_wunsch.h"
#include <iostream>
#include <fstream>
#include <chrono>

std::string read_sequence(const std::string &filename) {
    std::ifstream file(filename);
    std::string sequence;
    if (file.is_open()) {
        std::getline(file, sequence);
        file.close();
    }
    return sequence;
}

void save_output(const Result &result, const std::string &filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "Aligned Sequence 1: " << result.aligned_seq1 << "\n";
        file << "Aligned Sequence 2: " << result.aligned_seq2 << "\n";
        file << "Score: " << result.score << "\n";
        file << "Score Matrix:\n";
        // for (const auto& row : result.score_matrix) {
        //     for (const auto& val : row) {
        //         file << val << " ";
        //     }
        //     file << "\n";
        // }
        // file.close();
    }
}

int main() {
    std::string seq1 = read_sequence("seq1.txt");
    std::string seq2 = read_sequence("seq2.txt");
    
    int gap_penalty = -2;
    int match_score = 1;
    int mismatch_penalty = -1;

    auto start = std::chrono::high_resolution_clock::now();
    Result result_seq = needleman_wunsch(seq1, seq2, gap_penalty, match_score, mismatch_penalty);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seq = end - start;
    std::cout << "Sequential Execution Time: " << elapsed_seq.count() << " seconds\n";

    start = std::chrono::high_resolution_clock::now();
    Result result_omp = needleman_wunsch_omp(seq1, seq2, gap_penalty, match_score, mismatch_penalty);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_omp = end - start;
    std::cout << "OpenMP Execution Time: " << elapsed_omp.count() << " seconds\n";

    save_output(result_seq, "output_seq.txt");
    save_output(result_omp, "output_omp.txt");

    std::ofstream comparison_file("output.txt");
    if (comparison_file.is_open()) {
        comparison_file << "Sequential Execution Time: " << elapsed_seq.count() << " seconds\n";
        comparison_file << "OpenMP Execution Time: " << elapsed_omp.count() << " seconds\n";
        comparison_file.close();
    }

    return 0;
}