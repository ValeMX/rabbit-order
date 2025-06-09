#include "graph_binary.h"

#include <fstream>
#include <iostream>
using namespace std;

GraphBinary::GraphBinary() : nNodes(0), nEdges(0), totalWeight(0) {}

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

    // Populate neighbours list
    for (unsigned int i = 0; i < nEdges; i++) {
        unsigned int src = edgeList[i].first;
        unsigned int dst = edgeList[i].second;

        neighboursList[src].push_back(dst);
        totalWeight++;
    }
}

bool GraphBinary::isRemote(unsigned int node) {
    return find(localNodes.begin(), localNodes.end(), node) == localNodes.end();
}

vector<unsigned int> GraphBinary::neighbours(unsigned int node) {
    if (isRemote(node)) {
        throw out_of_range("Neighbours: Node index out of range");
    }

    return neighboursList[node];  // Return the neighbours of the node
}

vector<unsigned int> GraphBinary::remoteNeighbours(unsigned int node) {
    if (isRemote(node)) {
        throw out_of_range("Remote neighbours: Node index out of range");
    }

    // TODO: use count
    vector<unsigned int> remoteNeighbours;
    for (const auto& neighbour : neighboursList[node]) {
        if (isRemote(neighbour)) {
            remoteNeighbours.push_back(neighbour);
        }
    }
    return remoteNeighbours;
}

unsigned int GraphBinary::degree(unsigned int node) {
    if (isRemote(node)) {
        throw out_of_range("Degree: Node index out of range");
    }

    if (neighboursList.find(node) == neighboursList.end()) {
        return 0;  // No neighbours found
    }

    return neighboursList[node].size();
}

unsigned int GraphBinary::selfLoops(unsigned int node) {
    if (isRemote(node)) {
        throw out_of_range("Self loops: Node index out of range");
    }

    // TODO: use count
    unsigned int sl = 0;
    for (const auto& neighbour : neighboursList[node]) {
        if (neighbour == node) {  // Check for self-loop
            sl++;                 // Increment self-loop weight
        }
    }
    return sl;
}

unsigned int GraphBinary::remoteDegree(unsigned int node) {
    if (!isRemote(node)) {
        throw out_of_range("Remote degree: Node index out of range");
    }

    if (remoteNeighboursList.find(node) == remoteNeighboursList.end()) {
        return 0;  // No neighbours found
    }

    return remoteNeighboursList[node].size();
}

unsigned int GraphBinary::remoteSelfLoops(unsigned int node) {
    if (!isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    // TODO: use count
    unsigned int sl = 0;
    for (const auto& neighbour : remoteNeighboursList[node]) {
        if (neighbour == node) {  // Check for self-loop
            sl++;                 // Increment self-loop weight
        }
    }
    return sl;
}

void GraphBinary::addEdge(unsigned int source, unsigned int destination, double weight) {
    // Add the edge to the remoteNeighboursList if it's a remote edge
    if (isRemote(source)) {
        remoteNeighboursList[source].push_back(destination);
    } else {
        neighboursList[source].push_back(destination);
    }

    totalWeight++;  // Update total weight
}
