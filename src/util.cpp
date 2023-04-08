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
    amino_acid_mass_table['A'] = 71.03711;
    amino_acid_mass_table['R'] = 156.10111;
    amino_acid_mass_table['N'] = 114.04293;
    amino_acid_mass_table['D'] = 115.02694;
    amino_acid_mass_table['C'] = 103.00919;
    amino_acid_mass_table['Q'] = 128.05858;
    amino_acid_mass_table['E'] = 129.04259;
    amino_acid_mass_table['G'] = 57.02146;
    amino_acid_mass_table['H'] = 137.05891;
    amino_acid_mass_table['I'] = 113.08406;
    amino_acid_mass_table['L'] = 113.08406;
    amino_acid_mass_table['K'] = 128.09496;
    amino_acid_mass_table['M'] = 131.04049;
    amino_acid_mass_table['F'] = 147.06841;
    amino_acid_mass_table['P'] = 97.05276;
    amino_acid_mass_table['S'] = 87.03203;
    amino_acid_mass_table['T'] = 101.04768;
    amino_acid_mass_table['W'] = 186.07931;
    amino_acid_mass_table['Y'] = 163.06333;
    amino_acid_mass_table['V'] = 99.06841;
    return amino_acid_mass_table;
}

std::unordered_map<char, std::vector<std::string>>
get_amino_acid_to_codon_table() {
    std::unordered_map<char, std::vector<std::string>> codon_table = {
        {'A', {"GCT", "GCC", "GCA", "GCG"}},
        {'R', {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}},
        {'N', {"AAT", "AAC"}},
        {'D', {"GAT", "GAC"}},
        {'C', {"TGT", "TGC"}},
        {'Q', {"CAA", "CAG"}},
        {'E', {"GAA", "GAG"}},
        {'G', {"GGT", "GGC", "GGA", "GGG"}},
        {'H', {"CAT", "CAC"}},
        {'I', {"ATT", "ATC", "ATA"}},
        {'L', {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}},
        {'K', {"AAA", "AAG"}},
        {'M', {"ATG"}},
        {'F', {"TTT", "TTC"}},
        {'P', {"CCT", "CCC", "CCA", "CCG"}},
        {'S', {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}},
        {'T', {"ACT", "ACC", "ACA", "ACG"}},
        {'W', {"TGG"}},
        {'Y', {"TAT", "TAC"}},
        {'V', {"GTT", "GTC", "GTA", "GTG"}},
        {'*', {"TAA", "TAG", "TGA"}}};

    return codon_table;
}

void split(std::string &s, std::string delimiter, std::vector<std::string> &v) {
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        v.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    v.push_back(s);
}

void read_fasta(std::string filepath, std::map<std::string, std::string> &m) {
    std::ifstream ifs(filepath);
    if (ifs.is_open()) {
        std::string line, header;
        bool is_content = false;
        while (std::getline(ifs, line)) {
            if (line[0] == '>') {
                header = line;
                m[header] = "";
            } else {
                m[header] += line;
            }
        }
    }
    ifs.close();
}

void convert_dna_strings_to_profile(
    std::vector<std::string> &v,
    std::vector<std::vector<int>> &profile_matrix) {
    // Initialize
    for (int i = 0; i < v[0].size(); i++) {
        profile_matrix[0].push_back(0);
        profile_matrix[1].push_back(0);
        profile_matrix[2].push_back(0);
        profile_matrix[3].push_back(0);
    }
    // Add base
    for (int i = 0; i < v.size(); i++) {
        for (int j = 0; j < v[0].size(); j++) {
            char elem = v[i][j];
            if (elem == 'A')
                profile_matrix[0][j] += 1;
            if (elem == 'C')
                profile_matrix[1][j] += 1;
            if (elem == 'G')
                profile_matrix[2][j] += 1;
            if (elem == 'T')
                profile_matrix[3][j] += 1;
        }
    }
}

void print_dna_profile_matrix(std::vector<std::vector<int>> &profile_matrix) {
    // A, C, G, T order.
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < profile_matrix[i].size(); j++) {
            std::cout << profile_matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
}

void make_consensus_sequence_from_dna_profile_matrix(
    std::vector<std::vector<int>> &profile_matrix,
    std::string &consensus_sequence) {
    for (int i = 0; i < profile_matrix[0].size(); i++) {
        char max_char = 'A';
        int max_val = profile_matrix[0][i];
        if (profile_matrix[1][i] > max_val) {
            max_char = 'C';
            max_val = profile_matrix[1][i];
        }
        if (profile_matrix[2][i] > max_val) {
            max_char = 'G';
            max_val = profile_matrix[2][i];
        }
        if (profile_matrix[3][i] > max_val) {
            max_char = 'T';
            max_val = profile_matrix[3][i];
        }
        consensus_sequence += max_char;
    }
}