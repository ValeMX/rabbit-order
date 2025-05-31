#ifndef COMMUNITY_H
#define COMMUNITY_H

#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;

#include "graph_binary.h"

class Community {
   public:
    GraphBinary g;  // Graph to compute communities on
    int size;       // Number of nodes in the graph

    vector<int> n2c;     // Node to community mapping
    vector<int> n2cNew;  // New node to community mapping
    vector<double> in;   // Internal weights for communities
    vector<double> tot;  // Incident weights to communities

    vector<double> neighbourWeights;          // Weights of neighbours
    vector<unsigned int> neighbourPositions;  // Positions of neighbours
    unsigned int neighbourLast;               // Last neighbour position

    int steps;         // Number of steps for the algorithm
    double threshold;  // Threshold for computing a step

    Community(char *inFile, char *inWeightsFile, int st, double thr);
    Community(GraphBinary gb, int st, double thr);

    void print() const;  // Print the community structure

    void insert(int node, int community, double weightNodeToCommunity);  // Insert a node into a community
    void remove(int node, int community, double weightNodeToCommunity);  // Remove a node from a community
    void neighbourCommunities(int node);                                 // Find communities of neighbours of a node

    // Compute the modularity gain for moving a node to a community (using simplified modularity formula)
    double modularityGain(int node, int community, double weightNodeToCommunity, double weight) const;
    double modularity() const;  // Compute the modularity of the current community structure

    bool step();  // Perform a single step of the community detection algorithm
};

#endif  // COMMUNITY_H
