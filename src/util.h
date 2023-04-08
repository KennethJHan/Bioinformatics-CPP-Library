#ifndef UTIL_H
#define UTIL_H

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

void reverse_string(std::string &str);
void complement_dna(std::string &str);
double calc_gc(std::string &str);
int calc_hamming_distance(std::string &str1, std::string &str2);
std::string &translate_rna_to_protein(std::string &str);
std::map<char, double> get_amino_acid_mass_table();
std::unordered_map<char, std::vector<std::string>>
get_amino_acid_to_codon_table();
void split(std::string &s, std::string delimiter, std::vector<std::string> &v);
void read_fasta(std::string filepath, std::map<std::string, std::string> &m);
void convert_dna_strings_to_profile(
    std::vector<std::string> &v, std::vector<std::vector<int>> &profile_matrix);
void print_dna_profile_matrix(std::vector<std::vector<int>> &profile_matrix);
void make_consensus_sequence_from_dna_profile_matrix(
    std::vector<std::vector<int>> &profile_matrix,
    std::string &consensus_sequence);
#endif