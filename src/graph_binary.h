#ifndef GRAPH_BINARY_H
#define GRAPH_BINARY_H

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;

class GraphBinary {
   public:
    unsigned int startingNode;  // Starting node for the graph (global index)
    unsigned int nNodes;        // Number of nodes
    unsigned long nEdges;       // Number of edges
    double totalWeight;         // Total weight of the graph

    vector<pair<unsigned int, unsigned int>> edgeList;  // List of edges in the graph (global indices)
    vector<double> weightList;                          // Weights of the edges (global indices)

    map<unsigned int, vector<pair<unsigned int, double>>> neighboursList;        // Adjacency list representation (local indices)
    map<unsigned int, vector<pair<unsigned int, double>>> remoteNeighboursList;  // Remote adjacency list representation (local indices)

    unordered_map<unsigned int, unsigned int> globalToLocal;  // Mapping from global node index to local index
    std::vector<unsigned int> localToGlobal;                  // Mapping from local node index to global index

    unordered_set<unsigned int> localNodes;   // Local nodes in the graph (local indices)
    unordered_set<unsigned int> remoteNodes;  // Remote nodes in the graph (local indices)

    GraphBinary();
    GraphBinary(char* inFile, char* weightsInFile);
    GraphBinary(unsigned int nNodes, unsigned long nEdges, double totalWeight, int* degrees, int* edges, double* weights);

    void renumber();  // Renumber the nodes in the graph
    void init();      // Initialize the graph

    bool isRemote(unsigned int node);                          // Check if a node is remote
    vector<unsigned int> neighbours(unsigned int node);        // Get neighbours of a node
    vector<unsigned int> remoteNeighbours(unsigned int node);  // Get remote neighbours of a node

    unsigned int nNeighbours(unsigned int node);  // Get the number of neighbors of a node

    double weightedDegree(unsigned int node);        // Compute the weighted degree of a node
    double selfLoops(unsigned int node);             // Compute the self-loops of a node
    double remoteWeightedDegree(unsigned int node);  // Compute the weighted degree of a remote node
    double remoteSelfLoops(unsigned int node);       // Compute the self-loops of a remote node

    void addEdge(unsigned int source, unsigned int destination, double weight);  // Add an edge to the graph (global indices)
};

#endif  // GRAPH_BINARY_H