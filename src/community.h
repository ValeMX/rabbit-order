#ifndef COMMUNITY_H
#define COMMUNITY_H

#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <vector>

#include "graph_binary.h"
#include "graph_extended.h"
using namespace std;

class Community {
   public:
    GraphBinary& g;     // Graph to compute communities on
    unsigned int size;  // Number of nodes in the graph

    map<unsigned int, unsigned int> n2c;     // Node to community mapping
    map<unsigned int, unsigned int> n2cNew;  // New node to community mapping after renumbering
    map<unsigned int, double> tot;           // Total weight of each community
    map<unsigned int, double> in;            // Internal weight of each community

    unordered_map<unsigned int, double> neighbourCommunitiesMap;  // Weights to remote communities

    int steps;         // Number of steps for the algorithm
    double threshold;  // Threshold for computing a step

    vector<unsigned int> remoteCommunities;  // Remote communities mapping
    vector<double> remoteWeights;            // Weights of remote communities

    Community(GraphBinary& gb, int st, double thr);
    void updateRemote(unsigned int node, unsigned int community, double degree);  // Update remote communities

    void print();  // Print the community structure

    void insert(int node, int community, double weightNodeToCommunity);                                   // Insert a node into a community
    void insert(int node, int community, double weightNodeToCommunity, double weight, double selfLoops);  // Insert a node into a community with weight
    void remove(int node, int community, double weightNodeToCommunity);                                   // Remove a node from a community
    void neighbourCommunities(int node);                                                                  // Find communities of neighbours of a node

    // Compute the modularity gain for moving a node to a community (using simplified modularity formula)
    double modularityGain(int community, double weightNodeToCommunity, double weight);
    double modularity();  // Compute the modularity of the current community structure

    bool step();  // Perform a single step of the community detection algorithm

    GraphBinary graph();
    GraphBinary graphNew();
};

#endif  // COMMUNITY_H
