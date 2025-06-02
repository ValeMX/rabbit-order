#ifndef PARTITIONER_H
#define PARTITIONER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

class Partitioner {
   public:
    vector<pair<unsigned int, unsigned int>> edgeList;  // List of edges in the graph
    vector<double> weights;                             // Weights of the edges

    unsigned int numberOfPartitions;
    vector<unsigned int> partitionMap;

    vector<vector<pair<unsigned int, unsigned int>>> partitionEdges;  // Edges in each partition
    vector<vector<double>> partitionWeights;                          // Weights of edges in each partition
    
    Partitioner(unsigned int numPartitions);

    void staticPartition(const char* inFile);
};

#endif  // PARTITIONER_H