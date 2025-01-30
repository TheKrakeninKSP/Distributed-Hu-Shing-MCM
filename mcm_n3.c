#include <stdio.h>
#include <limits.h>
#include <time.h>

//compile command
//gcc mcm_n3.c -pg -o mcm_n3 && ./mcm_n3 && gprof mcm_n3 gmon.out > gprof_report.txt && cat gprof_report.txt

// Function declarations
void printOptimalParens(int *s, int i, int j, int n);
void matrixChainMultiplication(int *p, int n);
void runMultipleTimes(int *dimensions, int size);

// Function to print the optimal parenthesization
void printOptimalParens(int *s, int i, int j, int n) {
    if (i == j) {
        printf("A%d", i);
    } else {
        printf("(");
        printOptimalParens(s, i, *((s + i * n) + j), n);
        printOptimalParens(s, *((s + i * n) + j) + 1, j, n);
        printf(")");
    }
}

// Function to implement matrix chain multiplication
void matrixChainMultiplication(int *p, int n) {
    n = n - 1;
    int m[n][n];
    int s[n][n];
    int i, j, k, L, q;
    
    for (i = 0; i < n; i++) {
        m[i][i] = 0;
    }
    
    for (L = 2; L <= n; L++) {
        for (i = 0; i < n - L + 1; i++) {
            j = i + L - 1;
            m[i][j] = INT_MAX;
            
            for (k = i; k < j; k++) {
                q = m[i][k] + m[k + 1][j] + p[i] * p[k + 1] * p[j + 1];
                
                if (q < m[i][j]) {
                    m[i][j] = q;
                    s[i][j] = k;
                }
            }
        }
    }
    
    // Uncomment for debugging
    /*printf("Minimum scalar multiplications: %d\n", m[0][n-1]);
    printf("Optimal Parenthesization: ");
    printOptimalParens((int *)s, 0, n-1, n);
    printf("\n");
    */
}

// Function to run the algorithm multiple times for better profiling
void runMultipleTimes(int *dimensions, int size) {
    for (int i = 0; i < 1000000; i++) {  // Run many times
        matrixChainMultiplication(dimensions, size);
    }
}

int main() {
    // Example with larger matrices for better profiling
    int dimensions[] = {30, 35, 15, 5, 10, 20, 25, 40, 15, 30, 45, 20, 25, 15, 45, 25, 10, 5, 35, 45, 25, 15, 10};
    int size = sizeof(dimensions) / sizeof(dimensions[0]);
    
    clock_t start = clock();
    runMultipleTimes(dimensions, size);
    clock_t end = clock();
    
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Total execution time: %f seconds\n", time_spent);
    
    return 0;
}
