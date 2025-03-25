#include <iostream>
#include <climits>
#include <bits/stdc++.h>
#include <cuda_runtime.h>
#include <chrono>
#include <cstdlib>
#include <cuda_runtime.h>

#define N 30000  // Maximum number of matrices

#define FILENAME "inp.txt"

using namespace std;

__global__ void computeMCM(long long *dp, long long *dims, long long n, long long len) {
    long long i = blockIdx.x * blockDim.x + threadIdx.x + 1;  // Start from 1
    long long j = i + len - 1;

    if (j >= n) return;  // Ensure within bounds

    long long minCost = LLONG_MAX;

    for (long long k = i; k < j; k++) {
        long long cost = dp[i * N + k] + dp[(k + 1) * N + j] + 
                         (long long)dims[i - 1] * dims[k] * dims[j];
        if (cost < minCost) minCost = cost;
    }

    dp[i * N + j] = minCost;  // Store result
}

long long matrixChainCUDA(long long dims[], long long n) {
    long long *h_dp, *d_dp;
    long long *d_dims;
    
    // Allocate host memory
    h_dp = new long long[N * N];

    // Initialize DP table on host
    for (long long i = 0; i < N; i++)
        for (long long j = 0; j < N; j++)
            h_dp[i * N + j] = (i == j) ? 0 : LLONG_MAX; // 0 for single matrices

    // Allocate device memory
    cudaMalloc(&d_dp, N * N * sizeof(long long));
    cudaMalloc(&d_dims, n * sizeof(long long));

    // Copy data to device
    cudaMemcpy(d_dp, h_dp, N * N * sizeof(long long), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dims, dims, n * sizeof(long long), cudaMemcpyHostToDevice);

    // Parallel computation
    for (long long len = 2; len < n; len++) {
    	long long threadsPerBlock = 16;
        long long numBlocks = (n + 255) / threadsPerBlock;
        computeMCM<<<numBlocks, threadsPerBlock>>>(d_dp, d_dims, n, len);
        cudaDeviceSynchronize();
    }

    // Copy results back to host
    cudaMemcpy(h_dp, d_dp, N * N * sizeof(long long), cudaMemcpyDeviceToHost);

    // Prlong long the minimum cost
    //std::cout << "Minimum cost of multiplication(GPU): " << h_dp[1 * N + (n - 1)] << std::endl;
	long long result = h_dp[1 * N + (n - 1)];
    // Free memory
    delete[] h_dp;
    cudaFree(d_dp);
    cudaFree(d_dims);
    return result;
}


long long matrixChainDP(long long* dims, long long n) {
    // Allocate memory dynamically for DP table
    long long** dp = new long long*[n];
    for (long long i = 0; i < n; i++) {
        dp[i] = new long long[n];
    }

    // Initialize DP table
    for (long long i = 0; i < n; i++) {
        for (long long j = 0; j < n; j++) {
            dp[i][j] = (i == j) ? 0 : LLONG_MAX;
        }
    }

    // DP computation
    for (long long len = 2; len < n; len++) {  // Chain length
        for (long long i = 1; i < n - len + 1; i++) {
            long long j = i + len - 1;
            dp[i][j] = LLONG_MAX;

            for (long long k = i; k < j; k++) {
                long long cost = dp[i][k] + dp[k + 1][j] + 
                                 (long long)dims[i - 1] * dims[k] * dims[j];

                if (cost < dp[i][j]) {
                    dp[i][j] = cost;
                }
            }
        }
    }

    long long result = dp[1][n - 1];  // Minimum cost for multiplying matrices from 1 to n-1

    // Free allocated memory
    for (long long i = 0; i < n; i++) {
        delete[] dp[i];
    }
    delete[] dp;
    //cout << "Minimum cost of multiplication(CPU): " << result << endl;
    return result;
}

void benchmark(long long* dims, int n) {
    // Benchmark CUDA
    cudaEvent_t startCUDA, stopCUDA;
    cudaEventCreate(&startCUDA);
    cudaEventCreate(&stopCUDA);

    cudaEventRecord(startCUDA);
    long long cudaResult = matrixChainCUDA(dims, n);
    cudaEventRecord(stopCUDA);

    cudaEventSynchronize(stopCUDA);
    float cudaTime;
    cudaEventElapsedTime(&cudaTime, startCUDA, stopCUDA);

    cudaEventDestroy(startCUDA);
    cudaEventDestroy(stopCUDA);
    
    // Benchmark CPU
    auto startCPU = std::chrono::high_resolution_clock::now();
    long long cpuResult = matrixChainDP(dims, n);
    auto endCPU = std::chrono::high_resolution_clock::now();
    double cpuTime = std::chrono::duration<double, milli>(endCPU - startCPU).count();
    
    double speedup = cpuTime / cudaTime;

    // Print results
    cout << "CPU MCM Result: " << cpuResult << " | Time: " << cpuTime << " ms" << endl;
    cout << "CUDA MCM Result: " << cudaResult << " | Time: " << cudaTime << " ms" << endl;
    cout << "Speedup Ratio: " << speedup << "x" << endl;
}

long long read_data(long long arr[N]) {
    std::ifstream file(FILENAME);
    long long n = 0;
    long long x;
    while (file >> x) {
        arr[n++] = x;
    }
    return n;
}

int main() {
    long long* dims = new long long[N];
    long long n = read_data(dims);
    benchmark(dims, n);
    return 0;
}

