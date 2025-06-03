#include "graph_binary.h"

#include <fstream>
#include <iostream>
using namespace std;

GraphBinary::GraphBinary() : nNodes(0), nEdges(0), totalWeight(0.0) {}

void GraphBinary::init() {
    remoteNodes.clear();

    std::unordered_set<unsigned int> visibleNodes;
    for (const auto& edge : edgeList) {
        unsigned int src = edge.first;
        unsigned int dst = edge.second;

        visibleNodes.insert(src);
        visibleNodes.insert(dst);
        startingNode = min(startingNode, src);
    }

    // Remote nodes are those that are not in localNodes
    for (const auto& node : visibleNodes) {
        if (isRemote(node)) {
            remoteNodes.push_back(node);
        }
    }

    nNodes = localNodes.size();
    nEdges = edgeList.size();
    totalWeight = 0.0;

    if (nEdges == 0) {
        return;  // No edges to process
    }

    // Populate neighbours list
    for (unsigned int i = 0; i < nEdges; i++) {
        unsigned int src = edgeList[i].first;
        unsigned int dst = edgeList[i].second;
        double weight = weightList.empty() ? 1.0 : weightList[i];

        neighboursList[src].emplace_back(dst, weight);
        totalWeight += weight;
    }
}

bool GraphBinary::isRemote(unsigned int node) {
    return find(localNodes.begin(), localNodes.end(), node) == localNodes.end();
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
    // Add the edge to the remoteNeighboursList if it's a remote edge
    if (isRemote(source)) {
        remoteNeighboursList[source].emplace_back(destination, weight);
    } else {
        neighboursList[source].emplace_back(destination, weight);
    }

    totalWeight += weight;  // Update total weight
}
