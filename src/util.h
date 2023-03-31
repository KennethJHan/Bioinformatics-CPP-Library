#ifndef UTIL_H
#define UTIL_H

#include <string>

void reverse_string(std::string &str);
void complement_dna(std::string &str);
double calc_gc(std::string &str);
int calc_hamming_distance(std::string &str1, std::string &str2);
#endif