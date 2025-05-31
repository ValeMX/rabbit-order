#ifndef GRAPH_EXTENDED_H
#define GRAPH_EXTENDED_H

#include <map>
#include <string>
#include <vector>

#include "custom_vector.h"
using namespace std;

class GraphExtended {
   public:
    unsigned int nNodes;   // Number of nodes
    unsigned long nEdges;  // Number of edges
    double totalWeight;    // Total weight of the graph

    // CSR representation
    CustomVector<unsigned long>* degrees;  // Cumulative degree for each node
    CustomVector<unsigned int>* edges;     // Edges for each node
    CustomVector<double>* weights;         // Weights for each edge

    map<int, vector<unsigned int>> remoteEdges;
    map<int, vector<double>> remoteWeights;

    GraphExtended();
    GraphExtended(unsigned int nNodes, unsigned long nEdges, double totalWeight, int* degrees, int* edges, double* weights);
    virtual ~GraphExtended();

    unsigned int nNeighbours(unsigned int node);        // Get the number of neighbors of a node
    unsigned int nRemoteNeighbours(unsigned int node);  // Get the number of remote neighbors of a node

    double weightedDegree(unsigned int node);        // Compute the weighted degree of a node
    double weightedDegreeRemote(unsigned int node);  // Compute the weighted degree of a remote node
    double selfLoops(unsigned int node);             // Compute the self-loops of a node

    pair<vector<unsigned int>::iterator, vector<double>::iterator> neighbours(unsigned int node);        // Get neighbours of a node
    pair<vector<unsigned int>::iterator, vector<double>::iterator> remoteNeighbours(unsigned int node);  // Get remote neighbours of a node

    void addRemoteEdges(map<int, vector<unsigned int>> re, map<int, vector<double>> rw);  // Add remote edges to the graph
    void cleanup();

    void print();
};

#endif  // GRAPH_EXTENDED_H