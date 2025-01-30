#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

// Structure to represent a vertex in the triangulation
typedef struct {
    int index;
    double angle;  // Angle in radians
} Vertex;

// Structure to represent an edge in the triangulation
typedef struct {
    int start;
    int end;
    double weight;
} Edge;

// Function to calculate angle between three points
double calculateAngle(double x1, double y1, double x2, double y2, double x3, double y3) {
    double angle1 = atan2(y1 - y2, x1 - x2);
    double angle2 = atan2(y3 - y2, x3 - x2);
    double angle = angle2 - angle1;
    if (angle < 0) angle += 2 * M_PI;
    return angle;
}

// Function to convert matrix dimensions to polygon vertices
Vertex* createPolygon(int* dimensions, int n, int* vertexCount) {
    *vertexCount = n - 1;
    Vertex* vertices = (Vertex*)malloc(*vertexCount * sizeof(Vertex));
    
    // Place vertices in a circle
    double radius = 1.0;
    double angleStep = 2 * M_PI / *vertexCount;
    
    for (int i = 0; i < *vertexCount; i++) {
        vertices[i].index = i;
        vertices[i].angle = i * angleStep;
    }
    
    return vertices;
}

// Function to calculate edge weight (cost of matrix multiplication)
double calculateEdgeWeight(int* dimensions, int start, int end) {
    return (double)(dimensions[start] * dimensions[start + 1] * dimensions[end + 1]);
}

// Function to check if edge intersects with another edge
int edgesIntersect(Vertex v1, Vertex v2, Vertex v3, Vertex v4) {
    double x1 = cos(v1.angle), y1 = sin(v1.angle);
    double x2 = cos(v2.angle), y2 = sin(v2.angle);
    double x3 = cos(v3.angle), y3 = sin(v3.angle);
    double x4 = cos(v4.angle), y4 = sin(v4.angle);
    
    // Check if lines intersect using cross product
    double d1 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
    double d2 = (x2 - x1) * (y4 - y1) - (y2 - y1) * (x4 - x1);
    double d3 = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3);
    double d4 = (x4 - x3) * (y2 - y3) - (y4 - y3) * (x2 - x3);
    
    return (d1 * d2 < 0) && (d3 * d4 < 0);
}

// Function to find the optimal triangulation using Hu-Shing algorithm
void huShingMCM(int* dimensions, int n) {
    int vertexCount;
    Vertex* vertices = createPolygon(dimensions, n, &vertexCount);
    
    // Create array to store edges of the triangulation
    int maxEdges = vertexCount * (vertexCount - 1) / 2;
    Edge* edges = (Edge*)malloc(maxEdges * sizeof(Edge));
    int edgeCount = 0;
    
    // Initialize with boundary edges
    for (int i = 0; i < vertexCount; i++) {
        edges[edgeCount].start = i;
        edges[edgeCount].end = (i + 1) % vertexCount;
        edges[edgeCount].weight = calculateEdgeWeight(dimensions, i, (i + 1) % vertexCount);
        edgeCount++;
    }
    
    // Find internal edges using greedy approach
    for (int len = 2; len < vertexCount - 1; len++) {
        for (int i = 0; i < vertexCount; i++) {
            int j = (i + len) % vertexCount;
            
            // Check if adding edge (i,j) creates valid triangulation
            int valid = 1;
            for (int k = 0; k < edgeCount; k++) {
                if (edges[k].start != i && edges[k].end != i &&
                    edges[k].start != j && edges[k].end != j) {
                    if (edgesIntersect(vertices[i], vertices[j],
                                     vertices[edges[k].start],
                                     vertices[edges[k].end])) {
                        valid = 0;
                        break;
                    }
                }
            }
            
            if (valid) {
                edges[edgeCount].start = i;
                edges[edgeCount].end = j;
                edges[edgeCount].weight = calculateEdgeWeight(dimensions, i, j);
                edgeCount++;
            }
        }
    }
    
    // Calculate total cost
    double totalCost = 0;
    for (int i = 0; i < edgeCount; i++) {
        totalCost += edges[i].weight;
    }
    
    // Print results
    printf("\nHu-Shing Algorithm Results:\n");
    printf("Total multiplication cost: %.0f\n", totalCost);
    printf("Triangulation edges:\n");
    for (int i = 0; i < edgeCount; i++) {
        printf("Edge %d: %d-%d (cost: %.0f)\n", 
               i + 1, edges[i].start, edges[i].end, edges[i].weight);
    }
    
    // Free allocated memory
    free(vertices);
    free(edges);
}

int main() {
    // Example usage with the same matrices as before
    int dimensions[] = {30, 35, 15, 5, 10, 20, 25};
    int n = sizeof(dimensions) / sizeof(dimensions[0]);
    
    printf("Matrix dimensions: ");
    for (int i = 0; i < n - 1; i++) {
        printf("%dx%d  ", dimensions[i], dimensions[i + 1]);
    }
    printf("\n");
    
    huShingMCM(dimensions, n);
    
    return 0;
}
