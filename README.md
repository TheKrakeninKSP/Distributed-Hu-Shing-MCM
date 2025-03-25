# Hu-Shing Matrix Chain Multiplication

This repository contains implementations of **Hu-Shing's Method for Matrix Chain Multiplication (MCM)** in three parallel computing approaches:

1. **Approach 1:** OpenMP-based parallel implementation in C++.
2. **Approach 2:** pThreads-based distributed implementation with minimized error rate.
3. **Approach 3:** CUDA-based parallel implementation of Dynamic Programming MCM algorithm.

## Overview
Matrix Chain Multiplication (MCM) aims to find the optimal order of multiplying a sequence of matrices to minimize computational cost. **Hu-Shing's Algorithm** provides an efficient way to approximate the optimal solution by decomposing the matrix sequence hierarchically.

This repository contains:
- **C++ OpenMP Implementation:** Uses OpenMP for parallel processing of matrix chain segments.
- **PThreads:** Uses Pthread library to distribute the array and compute optimal produce for each sub-chain and combine results. Has less than 0.0001% error for large matrices.
- **CUDA Implementation:** Uses GPU acceleration to speed up computations of the Dynamic-Programming based solution to Matrix Chain Multiplication.

## Features
- Efficient **hierarchical decomposition** using Hu-Shing’s method.
- **Parallelized** computation for large matrix chains.
- **Optimized memory usage** to handle large-scale matrix multiplications.

## Requirements
### For OpenMP Implementation:
- GCC or Clang with OpenMP support.
- Linux/Windows (with MinGW for Windows).

### For CUDA Implementation:
- NVIDIA GPU with CUDA Compute Capability 3.0 or higher.
- NVIDIA CUDA Toolkit.
