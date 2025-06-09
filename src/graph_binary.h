#ifndef GRAPH_BINARY_H
#define GRAPH_BINARY_H

#include <algorithm>
#include <climits>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;

class GraphBinary {
   public:
    unsigned int startingNode;  // Starting node for the graph
    unsigned int nNodes;        // Number of nodes
    unsigned long nEdges;       // Number of edges
    unsigned long totalWeight;  // Total weight of the graph

    vector<pair<unsigned int, unsigned int>> edgeList;  // List of edges in the graph

    map<unsigned int, vector<unsigned int>> neighboursList;        // Adjacency list representation
    map<unsigned int, vector<unsigned int>> remoteNeighboursList;  // Remote adjacency list representation

    vector<unsigned int> localNodes;   // Local nodes in the graph
    vector<unsigned int> remoteNodes;  // Remote nodes in the graph

    GraphBinary();
    GraphBinary(char* inFile, char* weightsInFile);
    GraphBinary(unsigned int nNodes, unsigned long nEdges, double totalWeight, int* degrees, int* edges, double* weights);

    void renumber();  // Renumber the nodes in the graph
    void init();      // Initialize the graph

    bool isRemote(unsigned int node);                          // Check if a node is remote
    vector<unsigned int> neighbours(unsigned int node);        // Get neighbours of a node
    vector<unsigned int> remoteNeighbours(unsigned int node);  // Get remote neighbours of a node

    unsigned int degree(unsigned int node);           // Compute the degree of a node
    unsigned int selfLoops(unsigned int node);        // Compute the self-loops of a node
    unsigned int remoteDegree(unsigned int node);     // Compute the weighted degree of a remote node
    unsigned int remoteSelfLoops(unsigned int node);  // Compute the self-loops of a remote node

    void addEdge(unsigned int source, unsigned int destination, double weight);  // Add an edge to the graph
};

#endif  // GRAPH_BINARY_H