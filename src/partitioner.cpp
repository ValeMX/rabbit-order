#include "partitioner.h"

#include <climits>

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
    unsigned int src, dst;
    double weight;
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

    partitionMap.resize(vid + 1, UINT_MAX);

    // Assign partitions to vertices
    for (unsigned int i = 0; i < partitionMap.size(); i++) {
        partitionMap[i] = i % numberOfPartitions;
    }

    partitionNodes.resize(numberOfPartitions);
    partitionEdges.resize(numberOfPartitions);
    partitionWeights.resize(numberOfPartitions);

    // Initialize partition nodes
    for (unsigned int i = 0; i < partitionMap.size(); i++) {
        unsigned int partitionId = partitionMap[i];
        if (partitionId < numberOfPartitions) {
            partitionNodes[partitionId].push_back(i);
        } else {
            cerr << "Partition ID out of range: " << partitionId << endl;
        }
    }

    auto weightIt = weights.begin();
    for (auto edgeIt = edgeList.begin(); edgeIt != edgeList.end(); ++edgeIt, ++weightIt) {
        unsigned int src = edgeIt->first;
        unsigned int dst = edgeIt->second;

        partitionEdges[partitionMap[src]].push_back(*edgeIt);
        partitionWeights[partitionMap[src]].push_back(*weightIt);
    }
}

void Partitioner::updatePartition(vector<int> communities) {
    vector<unsigned int> newPartitionMap(communities.size(), UINT_MAX);
    for (unsigned int i = 0; i < communities.size(); i++) {
        newPartitionMap[i] = i % numberOfPartitions;  // Reassign partitions based on communities
    }
    partitionMap = move(newPartitionMap);
}

unsigned int Partitioner::owner(unsigned int node) {
    return partitionMap[node];
}