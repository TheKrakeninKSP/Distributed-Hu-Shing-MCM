# Hu-Shing Matrix Chain Multiplication

This repository contains implementations of **Hu-Shing's Method for Matrix Chain Multiplication (MCM)** in two parallel computing approaches:

1. **Approach 1:** OpenMP-based parallel implementation in C++.
2. **Approach 2:** pThread-based parallel implementation in C++.
3. **Approach 3:** CUDA-based parallel implementation in Cuda-C++.

## Overview
Matrix Chain Multiplication (MCM) aims to find the optimal order of multiplying a sequence of matrices to minimize computational cost. **Hu-Shing's Algorithm** provides an efficient way to approximate the optimal solution by decomposing the matrix sequence hierarchically.

This repository contains:
- **C++ OpenMP Implementation:** Uses OpenMP for parallel processing of matrix chain segments.
- **pThreads Implementation:** Uses POSIX Thread library for distributed computing of matrix chain segments.
- **CUDA Implementation:** Uses GPU acceleration to further speed up computations.

## Features
- Efficient **hierarchical decomposition** using Hu-Shing’s method.
- **Parallelized** computation for large matrix chains.
- **Optimized memory usage** to handle large-scale matrix multiplications.


### Build OpenMP Version
```bash
g++ -fopenmp hu_shing_mcm_openmp.cpp -o mcm_openmp
./mcm_openmp
```

### Build CUDA Version
```bash
nvcc hu_shing_mcm_cuda.cu -o mcm_cuda
./mcm_cuda
```

## Performance Comparison
The CUDA implementation is expected to significantly outperform the OpenMP version for large matrix chains due to GPU acceleration. The following benchmarks compare execution times:
| Approach | Execution Time (N = 1000) |
|----------|-------------------------|
| OpenMP   | X.XX seconds            |
| CUDA     | Y.YY seconds            |
