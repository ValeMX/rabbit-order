#include "partitioner.h"

Partitioner::Partitioner(unsigned int numPartitions) {
    numberOfPartitions = numPartitions;
    partitionMap.resize(0);
}

void Partitioner::staticPartition(const char* inFile) {
    // Clear previous data
    edgeList.clear();
    weights.clear();
    partitionMap.clear();
    partitionEdges.clear();
    partitionWeights.clear();

    ifstream infile(inFile);
    if (!infile.is_open()) {
        cerr << "Error opening input file: " << inFile << endl;
        return;
    }

    unsigned int vid = 0;
    edgeList.resize(0);
    weights.resize(0);

    // Read edges and weights from the input file
    unsigned int src, dst, weight;
    string line;
    while (getline(infile, line)) {
        istringstream iss(line);
        if (iss >> src >> dst) {
            if (iss >> weight) {
                weights.push_back(weight);
            } else {
                weights.push_back(1.);
            }
            edgeList.emplace_back(src, dst);
            vid = max(vid, max(src, dst));  // Track the maximum vertex ID
        }
    }

    partitionMap.resize(vid + 1, -1);

    // Assign partitions to vertices
    for (unsigned int i = 0; i < partitionMap.size(); i++) {
        partitionMap[i] = i % numberOfPartitions;
    }

    partitionEdges.resize(numberOfPartitions);
    partitionWeights.resize(numberOfPartitions);

    auto weightIt = weights.begin();
    for (auto edgeIt = edgeList.begin(); edgeIt != edgeList.end(); ++edgeIt, ++weightIt) {
        unsigned int src = edgeIt->first;
        unsigned int dst = edgeIt->second;

        if ((src >= partitionMap.size()) || (dst >= partitionMap.size())) {
            cerr << " ERROR : " << dst << " : " << src << " " << endl;
            continue;
        }

        partitionEdges[partitionMap[src]].push_back(*edgeIt);
        partitionWeights[partitionMap[src]].push_back(*weightIt);
    }

    cerr << "Partitioning done" << endl;
}