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

std::string &translate_rna_to_protein(std::string &str) {
    std::map<std::string, std::string> codon_table = {
        {"UUU", "F"}, {"UUC", "F"}, {"UUA", "L"}, {"UUG", "L"}, {"CUU", "L"},
        {"CUC", "L"}, {"CUA", "L"}, {"CUG", "L"}, {"AUU", "I"}, {"AUC", "I"},
        {"AUA", "I"}, {"AUG", "M"}, {"GUU", "V"}, {"GUC", "V"}, {"GUA", "V"},
        {"GUG", "V"}, {"UCU", "S"}, {"UCC", "S"}, {"UCA", "S"}, {"UCG", "S"},
        {"CCU", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"}, {"ACU", "T"},
        {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"}, {"GCU", "A"}, {"GCC", "A"},
        {"GCA", "A"}, {"GCG", "A"}, {"UAU", "Y"}, {"UAC", "Y"}, {"UAA", "*"},
        {"UAG", "*"}, {"CAU", "H"}, {"CAC", "H"}, {"CAA", "Q"}, {"CAG", "Q"},
        {"AAU", "N"}, {"AAC", "N"}, {"AAA", "K"}, {"AAG", "K"}, {"GAU", "D"},
        {"GAC", "D"}, {"GAA", "E"}, {"GAG", "E"}, {"UGU", "C"}, {"UGC", "C"},
        {"UGA", "*"}, {"UGG", "W"}, {"CGU", "R"}, {"CGC", "R"}, {"CGA", "R"},
        {"CGG", "R"}, {"AGU", "S"}, {"AGC", "S"}, {"AGA", "R"}, {"AGG", "R"},
        {"GGU", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"}};
    std::string ptn = "";
    std::string codon;
    for (int i = 0; i < str.size(); i += 3) {
        codon = str.substr(i, 3);
        if (codon_table.find(codon) != codon_table.end())
            ptn += codon_table[codon];
        else
            ptn += "-";
    }
    return ptn;
}