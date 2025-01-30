#include <stdio.h>
#include <limits.h>
#include <assert.h>


//compile command
//gcc mcm_n3_2.c -fprofile-arcs -ftest-coverage -o mcm_n3_2 && ./mcm_n3_2 && gcov mcm_n3_2 && cat mcm_n3_2.c.gcov

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

void matrixChainMultiplication(int *p, int n) {
    if (n < 2) {
        printf("Need at least 2 matrices to multiply\n");
        return;
    }
    
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
    
    printf("Minimum scalar multiplications: %d\n", m[0][n-1]);
    printf("Optimal Parenthesization: ");
    printOptimalParens((int *)s, 0, n-1, n);
    printf("\n");
}

// Test function to verify results
void runTests() {
    // Test Case 1: Basic case with 3 matrices
    printf("\nTest Case 1: Basic case with 3 matrices\n");
    int dimensions1[] = {10, 20, 30, 40};
    matrixChainMultiplication(dimensions1, 4);
    
    // Test Case 2: Larger case
    printf("\nTest Case 2: Larger case\n");
    int dimensions2[] = {30, 35, 15, 5, 10, 20, 25, 15, 35, 75};
    matrixChainMultiplication(dimensions2, 10);
    
    // Test Case 3: Edge case - single matrix
    printf("\nTest Case 3: Edge case - single matrix\n");
    int dimensions3[] = {10, 20};
    matrixChainMultiplication(dimensions3, 2);
    
    // Test Case 4: Edge case - empty input
    printf("\nTest Case 4: Edge case - empty input\n");
    int dimensions4[] = {10};
    matrixChainMultiplication(dimensions4, 1);
}

int main() {
    printf("Running matrix chain multiplication tests...\n");
    runTests();
    return 0;
}
