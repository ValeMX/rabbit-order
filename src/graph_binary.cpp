#include "graph_binary.h"

#include <fstream>
#include <iostream>
using namespace std;

GraphBinary::GraphBinary() : nNodes(0), nEdges(0), totalWeight(0.0) {}

// TODO: weightsInFile is not used in this constructor, consider removing it or implementing its functionality

void GraphBinary::renumber() {
    // Renumber the nodes in the graph
    globalToLocal.clear();
    localToGlobal.clear();
    localNodes.clear();
    remoteNodes.clear();

    // Compute all nodes to remap
    std::unordered_set<unsigned int> visibleNodes;
    for (const auto& edge : edgeList) {
        unsigned int src = edge.first;
        unsigned int dst = edge.second;

        visibleNodes.insert(src);
        visibleNodes.insert(dst);
        startingNode = min(startingNode, src);
    }

    // Create mappings from global to local indices
    unsigned int localId = 0;
    for (const auto& node : visibleNodes) {
        globalToLocal[node] = localId++;
        localToGlobal.push_back(node);
    }

    // Local nodes are those that are sources in the edge list
    for (const auto& edge : edgeList) {
        localNodes.insert(globalToLocal[edge.first]);
    }

    // Remote nodes are those that are not in localNodes
    for (const auto& node : visibleNodes) {
        if (localNodes.find(globalToLocal[node]) == localNodes.end()) {
            remoteNodes.insert(globalToLocal[node]);
        }
    }
}

void GraphBinary::init() {
    renumber();  // Ensure nodes are renumbered before initialization

    nNodes = localNodes.size();
    nEdges = edgeList.size();
    totalWeight = 0.0;

    if (nEdges == 0) {
        return;  // No edges to process
    }

    // Populate neighbours list
    for (unsigned int i = 0; i < nEdges; i++) {
        unsigned int globalSrc = edgeList[i].first;
        unsigned int globalDst = edgeList[i].second;
        double weight = weightList.empty() ? 1.0 : weightList[i];

        unsigned int localSrc = globalToLocal[globalSrc];
        unsigned int localDst = globalToLocal[globalDst];

        neighboursList[localSrc].emplace_back(localDst, weight);
    }
}

bool GraphBinary::isRemote(unsigned int node) {
    return localNodes.find(node) == localNodes.end();
}

vector<unsigned int> GraphBinary::neighbours(unsigned int node) {
    if (isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    vector<unsigned int> neighbourList;
    for (const auto& neighbour : neighboursList[node]) {
        neighbourList.push_back(neighbour.first);
    }
    return neighbourList;
}

vector<unsigned int> GraphBinary::remoteNeighbours(unsigned int node) {
    if (isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    vector<unsigned int> remoteNeighbours;
    for (const auto& neighbour : neighboursList[node]) {
        if (isRemote(neighbour.first)) {
            remoteNeighbours.push_back(neighbour.first);
        }
    }
    return remoteNeighbours;
}

unsigned int GraphBinary::nNeighbours(unsigned int node) {
    if (isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    if (neighboursList.find(node) == neighboursList.end()) {
        return 0;  // No neighbours found
    }

    return neighboursList[node].size();
}

double GraphBinary::weightedDegree(unsigned int node) {
    if (isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    double degree = 0.0;
    for (const auto& neighbour : neighboursList[node]) {
        degree += neighbour.second;  // Sum the weights of the neighbours
    }
    return degree;
}

double GraphBinary::selfLoops(unsigned int node) {
    if (isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    double selfLoopWeight = 0.0;
    for (const auto& neighbour : neighboursList[node]) {
        if (neighbour.first == node) {  // Check for self-loop
            selfLoopWeight += neighbour.second;
        }
    }
    return selfLoopWeight;
}

double GraphBinary::remoteWeightedDegree(unsigned int node) {
    if (!isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    double degree = 0.0;
    for (const auto& neighbour : remoteNeighboursList[node]) {
        degree += neighbour.second;  // Sum the weights of the neighbours
    }
    return degree;
}

double GraphBinary::remoteSelfLoops(unsigned int node) {
    if (!isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    double selfLoopWeight = 0.0;
    for (const auto& neighbour : remoteNeighboursList[node]) {
        if (neighbour.first == node) {  // Check for self-loop
            selfLoopWeight += neighbour.second;
        }
    }
    return selfLoopWeight;
}

void GraphBinary::addEdge(unsigned int source, unsigned int destination, double weight) {
    // Ideally, source is already in globalToLocal, if not, we add it
    if (globalToLocal.find(source) == globalToLocal.end()) {
        unsigned int newId = localToGlobal.size();
        globalToLocal[source] = newId;
        localToGlobal.push_back(source);
    }

    if (globalToLocal.find(destination) == globalToLocal.end()) {
        unsigned int newId = localToGlobal.size();
        globalToLocal[destination] = newId;
        localToGlobal.push_back(destination);
    }

    unsigned int localSrc = globalToLocal[source];
    unsigned int localDst = globalToLocal[destination];

    // Add the edge to the remoteNeighboursList if it's a remote edge
    if (isRemote(localSrc)) {
        remoteNeighboursList[localSrc].emplace_back(localDst, weight);
    } else {
        neighboursList[localSrc].emplace_back(localDst, weight);
    }
}
