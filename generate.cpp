#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <random> //c++ library for 64-bit random values

using namespace std;

#define FILENAME "inp.txt"
#define N 40 * 1000 * 1000

void generate_sequence() {
    std::ofstream file(FILENAME);
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<long long> dist(5, 150);
    
    size_t size = N;
    for (size_t i = 0; i < size; ++i) {
        file << dist(rng) << "\n";
    }
}

int main() {
    generate_sequence();
    return 1;
}
