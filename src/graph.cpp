#include "graph.h"

#include <fstream>
#include <iostream>
using namespace std;

Graph::Graph() : nNodes(0), nEdges(0), totalWeight(0) {}

void Graph::init() {
    neighboursList.clear();  // Clear the neighbours list
    remoteNodes.clear();

    unordered_set<unsigned int> visibleNodes;
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

    minNode = min(minNode, *min_element(localNodes.begin(), localNodes.end()));
    maxNode = max(maxNode, *max_element(localNodes.begin(), localNodes.end()));

    nNodes = localNodes.size();
    nEdges = edgeList.size();
    neighboursList.resize(maxNode + 1);  // Resize the neighbours list to the number of local nodes

    // Populate neighbours list
    for (unsigned int i = 0; i < nEdges; i++) {
        unsigned int src = edgeList[i].first;
        unsigned int dst = edgeList[i].second;

        neighboursList[src].push_back(dst);
    }
}

void Graph::coarsen(Graph& g, Partitioner& p, int rank, vector<int>& n2c, map<int, vector<unsigned int>>& c2n, vector<int>& new2old, map<int, int>& old2new) {
    int size = c2n.size();        // Get the number of communities
    neighboursList.resize(size);  // Resize the neighbours list to the number of communities
    weights.resize(size);         // Resize the weights vector to the number of communities
    nNodes = localNodes.size();   // Update the number of local nodes
    totalWeight = g.totalWeight;  // Update the total weight of the graph
    if (g.weights.empty()) {
        for (const auto& community : localNodes) {
            vector<int> neighbourCommunities(size, 0);  // Initialize neighbour communities vector
            for (const auto& node : c2n[new2old[community]]) {
                for (const auto& neighbour : g.neighboursList[node]) {
                    int neighbourCommunity = n2c[neighbour];              // Get the community of the neighbour
                    neighbourCommunities[old2new[neighbourCommunity]]++;  // Increment the count for the neighbour community
                }
            }

            for (int i = 0; i < neighbourCommunities.size(); ++i) {
                if (neighbourCommunities[i] > 0) {
                    neighboursList[community].push_back(i);                 // Add the neighbour community to the list
                    weights[community].push_back(neighbourCommunities[i]);  // Add the weight of the neighbour community
                }
            }
        }
    } else {
        for (const auto& community : localNodes) {
            vector<int> neighbourCommunities(old2new.size(), 0);  // Initialize neighbour communities vector
            for (const auto& node : c2n[new2old[community]]) {
                auto it = g.weights[node].begin();
                for (const auto& neighbour : g.neighboursList[node]) {
                    int neighbourCommunity = n2c[neighbour];                   // Get the community of the neighbour
                    neighbourCommunities[old2new[neighbourCommunity]] += *it;  // Increment the weight for the neighbour community
                    ++it;                                                      // Move to the next weight
                }
            }

            for (int i = 0; i < neighbourCommunities.size(); ++i) {
                if (neighbourCommunities[i] > 0) {
                    neighboursList[community].push_back(i);                 // Add the neighbour community to the list
                    weights[community].push_back(neighbourCommunities[i]);  // Add the weight of the neighbour community
                }
            }
        }
    }
}

bool Graph::isRemote(unsigned int node) {
    return find(localNodes.begin(), localNodes.end(), node) == localNodes.end();
}

bool Graph::isCollected(unsigned int node) {
    return localNodes.count(node) > 0 || remoteNodes.count(node) > 0;
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
    if (!isCollected(node)) {
        throw out_of_range("Degree: Node index out of range (node: " + to_string(node) + ")");
    }

    if (weights.empty()) {
        return neighboursList[node].size();  // If weights are empty, return the size of neighbours list
    } else {
        return accumulate(weights[node].begin(), weights[node].end(), 0);  // Sum the weights for the degree
    }
}

unsigned int Graph::selfLoops(unsigned int node) {
    if (!isCollected(node)) {
        throw out_of_range("Self loops: Node index out of range");
    }

    unsigned int sl = 0;
    if (weights.empty()) {
        for (const auto& neighbour : neighboursList[node]) {
            if (neighbour == node) {  // Check for self-loop
                sl++;                 // Increment self-loop weight
            }
        }
    } else {
        if (weights[node].size() != neighboursList[node].size()) {
            throw runtime_error("Self loops: Weights and neighbours list sizes do not match for node " + to_string(node));
        }
        for (size_t i = 0; i < neighboursList[node].size(); ++i) {
            if (neighboursList[node][i] == node) {  // Check for self-loop
                sl += weights[node][i];             // Increment self-loop weight
            }
        }
    }

    return sl;  // Return the total self-loops for the node
}

void Graph::addEdge(unsigned int source, unsigned int destination) {
    neighboursList[source].push_back(destination);  // Add destination to the source's neighbours
}
