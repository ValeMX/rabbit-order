#ifndef GRAPH_H
#define GRAPH_H

#include <string>
#include <vector>
using namespace std;

class GraphBinary {
   public:
    unsigned int nNodes;  // Number of nodes
    unsigned int nEdges;  // Number of edges
    double totalWeight;   // Total weight of the graph

    // CSR representation
    vector<unsigned long> degrees;  // Cumulative degree for each node
    vector<unsigned int> edges;     // Edges for each node
    vector<double> weights;         // Weights for each edge

    GraphBinary();
    GraphBinary(char* inFile, char* weightsInFile);
    GraphBinary(int nNodes, int nEdges, double totalWeight, int* degrees, int* edges, double* weights);

    unsigned int nNeighbours(unsigned int node);  // Get the number of neighbors of a node
    double weightedDegree(unsigned int node);     // Compute the weighted degree of a node
    double selfLoops(unsigned int node);          // Compute the self-loops of a node

    pair<vector<unsigned int>::iterator, vector<double>::iterator> neighbours(unsigned int node);  // Get neighbours of a node
};

#endif  // GRAPH_H