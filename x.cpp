#include <iostream>
#include <climits>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <random>

#define ll long long
#define testsize 5000
#define NUM_TRIALS 1  // Average over multiple trials

// Matrix Chain Multiplication function
void matrixChainMultiplication(const std::vector<int>& p, int num_threads) {
    int n = p.size();
    std::vector<std::vector<int>> dp(n, std::vector<int>(n, 0));

    // Compute the cost in a bottom-up manner
    for (int len = 2; len < n; len++) { // Length of matrix chain
        omp_set_num_threads(num_threads);
        #pragma omp parallel for
        for (int i = 1; i < n - len + 1; i++) {
            int j = i + len - 1;
            dp[i][j] = INT_MAX;

            for (int k = i; k < j; k++) {
                int cost = dp[i][k] + dp[k + 1][j] + p[i - 1] * p[k] * p[j];

                if (cost < dp[i][j]) {
                    dp[i][j] = cost;
                }
            }
        }
    }

    // Output the minimum multiplication cost
    //cout << "Minimum scalar multiplications: " << dp[1][n - 1] << endl;
}

void benchmark(void(*func)(const std::vector<int>&, int), std::vector<int>& A, std::string method) {
    std::vector<int> thread_counts = {1, 2, 4, 6, 8, 10, 12, 16, 20, 32, 64};
    //std::vector<int> thread_counts = {32, 64};
    int num_tests = thread_counts.size();
    std::vector<double> times(num_tests);
    double base_time;
    
    std::cout<<method<<"\n";
    printf("Threads\tTime (s)\n");
    for (int i = 0; i < num_tests; i++) {
        int num_threads = thread_counts[i];
        double total_time = 0.0;
        
        for (int j = 0; j < NUM_TRIALS; j++) {
            double start = omp_get_wtime();
            
            //Test segment code
            func(A, num_threads);
            
            //End of Test segment
            double end = omp_get_wtime();
            total_time += (end - start);
        }
        times[i] = total_time / NUM_TRIALS;
        if (i == 0) base_time = times[i]; // Time with 1 thread
        printf("%d\t%.6f\n", num_threads, times[i]);
    }

    // Speedup calculation
    std::vector<double> parallel_fracs(num_tests);
    std::cout<<"\nSpeedup vs Processors:\n";
    std::cout<<"Threads\tSpeedup\n";
    for (int i = 0; i < num_tests; i++) {
        // Estimate Parallelization Fraction (Amdahl's Law)
        if(i>0) parallel_fracs[i] = ((times[i] / base_time) - 1.0) / ((1.0 / thread_counts[i]) - 1.0);
        printf("%d\t%.6f\t", thread_counts[i], base_time / times[i]);
        //printf("%.6f, ", base_time / times[i]);
        printf("Parallelization Fraction (Amdahl's Law): %.3f\n", parallel_fracs[i]);
    }
    
    printf("Estimated Avg Parallelization Fraction (Amdahl's Law): %.3f\n", std::accumulate(parallel_fracs.begin(), parallel_fracs.end(), 0.0)/(num_tests - 1));

}

int main() {
    std::vector<int> arr(testsize);
    for(int i=0;i<testsize;i++) {
      arr[i] = rand()% (100 + 1 - 5) + 5;
    }

    benchmark(matrixChainMultiplication, arr, "MCM");

    return 0;
}

