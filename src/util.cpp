#include "util.h"

void reverse_string(std::string &str) {
    int n = str.length();
    for (int i = 0; i < n / 2; i++) {
        char temp = str[i];
        str[i] = str[n - i - 1];
        str[n - i - 1] = temp;
    }
}

void complement_dna(std::string &str) {
    std::string complement_str;
    for (char c : str) {
        if (c == 'A')
            complement_str += 'T';
        else if (c == 'C')
            complement_str += 'G';
        else if (c == 'G')
            complement_str += 'C';
        else if (c == 'T')
            complement_str += 'A';
    }
    str = complement_str;
}

double calc_gc(std::string &str) {
    double all_base, gc_base;
    for (char c : str) {
        if (c == 'G' || c == 'C')
            gc_base++;
        all_base++;
    }
    return gc_base / all_base * 100;
}

int calc_hamming_distance(std::string &str1, std::string &str2) {
    int hamming_distance = 0;
    for (int i = 0; i < str1.size(); i++) {
        char c1 = str1[i];
        char c2 = str2[i];
        if (c1 != c2)
            hamming_distance++;
    }

    return hamming_distance;
}