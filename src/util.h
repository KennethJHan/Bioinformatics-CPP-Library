#ifndef UTIL_H
#define UTIL_H

#include <map>
#include <string>

void reverse_string(std::string &str);
void complement_dna(std::string &str);
double calc_gc(std::string &str);
int calc_hamming_distance(std::string &str1, std::string &str2);
std::string &translate_rna_to_protein(std::string &str);
#endif