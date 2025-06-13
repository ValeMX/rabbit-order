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

void Graph::coarse(Graph& g, Partitioner& p, int rank, vector<int>& n2c, map<int, vector<unsigned int>>& c2n) {
    // Firstly renumber the communities
    map<int, int> old2new;  // Map to store new community indices
    vector<int> new2old;    // Vector to store community IDs
    for (const auto& community : c2n) {
        old2new[community.first] = old2new.size();  // Assign new indices to communities
        new2old.push_back(community.first);         // Store the old community ID
        if (p.owner(community.first) == rank)
            localNodes.insert(old2new[community.first]);  // Insert the new community index into local nodes if owned by this rank
    }

    p.updatePartition(new2old);  // Update the partition map with the new community indices

    neighboursList.resize(old2new.size());  // Resize the neighbours list to the number of communities
    weights.resize(old2new.size());         // Resize the weights vector to the number of communities
    nNodes = localNodes.size();             // Update the number of local nodes
    totalWeight = 0;                        // Reset total weight
    if (g.weights.empty()) {
        for (const auto& community : localNodes) {
            if (community >= new2old.size()) {
                throw out_of_range("Community index out of range in new2old");
            }

            vector<int> neighbourCommunities(old2new.size(), 0);  // Initialize neighbour communities vector
            for (const auto& node : c2n[new2old[community]]) {
                if (!g.isCollected(node)) {
                    std::cerr << "Rank: " << rank << " Invalid node index: " << node << " (max " << g.neighboursList.size() - 1 << ")" << std::endl;
                }

                if (node >= g.neighboursList.size()) {
                    std::cerr << "Rank: " << rank << " Invalido node index: " << node << " (max " << g.neighboursList.size() - 1 << ")" << std::endl;
                }

                for (const auto& neighbour : g.neighboursList[node]) {
                    if (neighbour >= n2c.size()) {
                        std::cerr << "Rank: " << rank << " Invalid neighbour index: " << neighbour << " (max " << n2c.size() - 1 << ")" << std::endl;
                    }
                    int neighbourCommunity = n2c[neighbour];  // Get the community of the neighbour
                    auto it = old2new.find(neighbourCommunity);
                    if (it == old2new.end()) {
                        std::cerr << "Rank: " << rank << "neighbourCommunity " << neighbourCommunity << " non trovato in old2new\n";
                    }
                    neighbourCommunities[old2new[neighbourCommunity]]++;  // Increment the count for the neighbour community
                    totalWeight++;                                        // Increment the total weight for each edge
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
                    totalWeight += *it;                                        // Increment the total weight for each edge
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
    totalWeight++;                                  // Update total weight
}
