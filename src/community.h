#ifndef COMMUNITY_H
#define COMMUNITY_H

#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>

#include "graph.h"
using namespace std;

class Community {
   public:
    Graph& g;           // Graph to compute communities on
    unsigned int size;  // Number of nodes in the graph

    vector<int> n2c;  // Node to community mapping
    vector<int> tot;  // Total weight of each community
    vector<int> in;   // Internal weight of each community

    unordered_map<unsigned int, unsigned int> neighbourCommunitiesMap;  // Weights to remote communities

    int steps;         // Number of steps for the algorithm
    double threshold;  // Threshold for computing a step

    Community(Graph& gb, int st, double thr);
    void resize();                                           // Resize the community structure based on the graph
    void updateRemote(int node, int community, int degree);  // Update remote communities

    void print();  // Print the community structure

    void insert(int node, int community, double weightNodeToCommunity);                                   // Insert a node into a community
    void insert(int node, int community, double weightNodeToCommunity, double weight, double selfLoops);  // Insert a node into a community with weight
    void remove(int node, int community, double weightNodeToCommunity);                                   // Remove a node from a community
    void neighbourCommunities(int node);                                                                  // Find communities of neighbours of a node

    int degreeN2C(int node);  // Compute the degree of a node to its community

    double modularityGain(int community, double weightNodeToCommunity, double weight);  // Compute the modularity gain for moving a node to a community
    double modularity();                                                                // Compute the modularity of the current community structure
    double modularity(int community);                                                   // Compute the modularity of a specific community

    bool step();  // Perform a single step of the community detection algorithm

    Graph graph();
    Graph graphNew();
};

#endif  // COMMUNITY_H
