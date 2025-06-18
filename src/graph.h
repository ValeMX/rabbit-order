#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <climits>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "partitioner.h"
using namespace std;

class Graph {
   public:
    unsigned int minNode;       // Starting node for the graph
    unsigned int maxNode;       // Ending node for the graph
    unsigned int nNodes;        // Number of nodes
    unsigned long nEdges;       // Number of edges
    unsigned long totalWeight;  // Total weight of the graph

    vector<pair<unsigned int, unsigned int>> edgeList;  // List of edges in the graph

    vector<vector<unsigned int>> neighboursList;  // Adjacency list representation
    vector<vector<int>> weights;

    set<unsigned int> localNodes;   // Local nodes in the graph
    set<unsigned int> remoteNodes;  // Remote nodes in the graph

    Graph();

    void init();                                                                                                                                           // Initialize the graph
    void coarse(Graph& g, Partitioner& p, int rank, vector<int>& n2c, map<int, vector<unsigned int>>& c2n, vector<int>& new2old, map<int, int>& old2new);  // Coarsen the graph based on communities

    bool isRemote(unsigned int node);                          // Check if a node is remote
    bool isCollected(unsigned int node);                       // Check if a node is collected
    vector<unsigned int> neighbours(unsigned int node);        // Get neighbours of a node
    vector<unsigned int> remoteNeighbours(unsigned int node);  // Get remote neighbours of a node

    unsigned int degree(unsigned int node);     // Compute the degree of a node
    unsigned int selfLoops(unsigned int node);  // Compute the self-loops of a node

    void addEdge(unsigned int source, unsigned int destination);  // Add an edge to the graph
};

#endif  // GRAPH_H