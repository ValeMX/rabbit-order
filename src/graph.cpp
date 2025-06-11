#include "graph.h"

#include <fstream>
#include <iostream>
using namespace std;

Graph::Graph() : nNodes(0), nEdges(0), totalWeight(0) {}

void Graph::init() {
    neighboursList.clear();  // Clear the neighbours list
    remoteNodes.clear();

    std::unordered_set<unsigned int> visibleNodes;
    maxNode = 0;
    minNode = UINT_MAX;
    for (const auto& edge : edgeList) {
        unsigned int src = edge.first;
        unsigned int dst = edge.second;

        visibleNodes.insert(src);
        visibleNodes.insert(dst);
        minNode = min(minNode, src);
        maxNode = max(maxNode, max(src, dst));
    }

    // Remote nodes are those that are not in localNodes
    for (const auto& node : visibleNodes) {
        if (isRemote(node)) {
            remoteNodes.insert(node);  // Add remote nodes to the set
        }
    }

    minNode = min(minNode, *min_element(localNodes.begin(), localNodes.end()));
    maxNode = max(maxNode, *max_element(localNodes.begin(), localNodes.end()));

    nNodes = localNodes.size();
    nEdges = edgeList.size();
    totalWeight = 0;
    neighboursList.resize(maxNode + 1);  // Resize the neighbours list to the number of local nodes

    // Populate neighbours list
    for (unsigned int i = 0; i < nEdges; i++) {
        unsigned int src = edgeList[i].first;
        unsigned int dst = edgeList[i].second;

        neighboursList[src].push_back(dst);
        totalWeight++;
    }
}

bool Graph::isRemote(unsigned int node) {
    return find(localNodes.begin(), localNodes.end(), node) == localNodes.end();
}

vector<unsigned int> Graph::neighbours(unsigned int node) {
    if (isRemote(node)) {
        throw out_of_range("Neighbours: Node index out of range");
    }

    return neighboursList[node];  // Return the neighbours of the node
}

vector<unsigned int> Graph::remoteNeighbours(unsigned int node) {
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

unsigned int Graph::degree(unsigned int node) {
    if (isRemote(node)) {
        throw out_of_range("Degree: Node index out of range");
    }

    return neighboursList[node].size();
}

unsigned int Graph::selfLoops(unsigned int node) {
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

unsigned int Graph::remoteDegree(unsigned int node) {
    if (!isRemote(node)) {
        throw out_of_range("Remote degree: Node index out of range");
    }

    // TODO: use a single function for both local and remote degrees
    return neighboursList[node].size();
}

unsigned int Graph::remoteSelfLoops(unsigned int node) {
    if (!isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    // TODO: use count, use a single function for both local and remote self-loops
    unsigned int sl = 0;
    for (const auto& neighbour : neighboursList[node]) {
        if (neighbour == node) {  // Check for self-loop
            sl++;                 // Increment self-loop weight
        }
    }
    return sl;
}

void Graph::addEdge(unsigned int source, unsigned int destination) {
    neighboursList[source].push_back(destination);  // Add destination to the source's neighbours
    totalWeight++;                                  // Update total weight
}
