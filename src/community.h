#ifndef COMMUNITY_H
#define COMMUNITY_H

#include "graph_binary.h"

class Community {
   public:
    GraphBinary g;       // Graph to compute communities on
    int size;            // Number of nodes in the graph
    vector<int> n2c;     // Node to community mapping
    vector<int> n2cNew;  // New node to community mapping
    vector<double> in;   // Internal weights for communities
    vector<double> tot;  // Incident weights to communities

    Community(string fileName = "");
};

#endif  // COMMUNITY_H
