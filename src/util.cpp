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

std::map<char, double> get_amino_acid_mass_table() {
    std::map<char, double> amino_acid_mass_table;
    amino_acid_mass_table['A'] = 89.09;
    amino_acid_mass_table['R'] = 174.20;
    amino_acid_mass_table['N'] = 132.12;
    amino_acid_mass_table['D'] = 133.10;
    amino_acid_mass_table['C'] = 121.16;
    amino_acid_mass_table['Q'] = 146.15;
    amino_acid_mass_table['E'] = 147.13;
    amino_acid_mass_table['G'] = 75.07;
    amino_acid_mass_table['H'] = 155.16;
    amino_acid_mass_table['I'] = 131.17;
    amino_acid_mass_table['L'] = 131.17;
    amino_acid_mass_table['K'] = 146.19;
    amino_acid_mass_table['M'] = 149.21;
    amino_acid_mass_table['F'] = 165.19;
    amino_acid_mass_table['P'] = 115.13;
    amino_acid_mass_table['S'] = 105.09;
    amino_acid_mass_table['T'] = 119.12;
    amino_acid_mass_table['W'] = 204.23;
    amino_acid_mass_table['Y'] = 181.19;
    amino_acid_mass_table['V'] = 117.15;
    return amino_acid_mass_table;
}