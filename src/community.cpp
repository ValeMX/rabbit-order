#include "community.h"

Community::Community(GraphBinary& gb, int st, double thr) : g(gb) {
    size = gb.nNodes;

    neighbourWeights.resize(size, -1);
    neighbourPositions.resize(size);
    neighbourLast = 0;

    for (const auto& node : g.localNodes) {
        n2c[node] = node;                    // Initialize each node to its own community
        tot[node] = g.weightedDegree(node);  // Total weight of the community is the weighted degree of the node
        in[node] = 2 * g.selfLoops(node);    // Internal weight is twice the self-loops of the node
    }

    steps = st;
    threshold = thr;
}

void Community::updateRemote(unsigned int node, unsigned int community, double degree) {
    // Update the remote community structure
    n2c[node] = community;  // Assign the node to the specified community
    tot[community] += degree;
    in[community] += 2 * g.remoteSelfLoops(node);
}

void Community::print() const {
    cout << "Community structure: (node community)" << endl;
    for (const auto& pair : n2c) {
        cout << pair.first << " " << pair.second << endl;
    }
    cout << "Total communities: " << tot.size() << endl;
    cout << "Modularity: " << modularity() << endl;
}

void Community::insert(int node, int community, double weightNodeToCommunity) {
    if (g.isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    tot[community] += g.weightedDegree(node);
    in[community] += 2 * weightNodeToCommunity + g.selfLoops(node);
    n2c[node] = community;  // Assign the node to the specified community
}

void Community::insert(int node, int community, double weightNodeToCommunity, double weight, double selfLoops) {
    if (g.globalToLocal.find(node) == g.globalToLocal.end()) {
        unsigned int newId = g.localToGlobal.size();
        g.globalToLocal[node] = newId;    // Add the node to the global to local mapping
        g.localToGlobal.push_back(node);  // Add the node to the local to global mapping
    }

    unsigned int localNode = g.globalToLocal[node];            // Get the local index of the node
    unsigned int localCommunity = g.globalToLocal[community];  // Get the local index of the community

    tot[localCommunity] += weight;                                // Update the total weight of the community
    in[localCommunity] += 2 * weightNodeToCommunity + selfLoops;  // Update the internal weight of the community
    n2c[localNode] = localCommunity;                              // Assign the node to the specified community
}

void Community::remove(int node, int community, double weightNodeToCommunity) {
    if (g.isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    tot[community] -= g.weightedDegree(node);
    in[community] -= 2 * weightNodeToCommunity + g.selfLoops(node);
    n2c[node] = -1;  // Remove the node from its community
}

void Community::neighbourCommunities(int node) {
    if (g.isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    for (unsigned int i = 0; i < neighbourLast; i++)
        neighbourWeights[neighbourPositions[i]] = -1;
    neighbourLast = 0;

    vector<unsigned int> n = g.neighbours(node);
    unsigned int deg = g.nNeighbours(node);

    neighbourPositions[0] = n2c[node];            // Ids of the neighbour community
    neighbourWeights[neighbourPositions[0]] = 0;  // Weights of the neighbour community
    neighbourLast = 1;                            // Number of neighbour communities

    for (unsigned int i = 0; i < deg; i++) {
        unsigned int neighbour = n[i];
        unsigned int neighbourCommunity = n2c[neighbour];
        double neighbourWeight = g.weightedDegree(neighbour);

        if (neighbour != node) {
            if (neighbourWeights[neighbourCommunity] == -1) {
                neighbourWeights[neighbourCommunity] = 0.;
                neighbourPositions[neighbourLast++] = neighbourCommunity;
            }
            neighbourWeights[neighbourCommunity] += neighbourWeight;
        }
    }
}

double Community::modularityGain(int node, int community, double weightNodeToCommunity, double weight) const {
    if (g.isRemote(node)) {
        throw out_of_range("Node or community index out of range");
    }

    double totCommunity = tot.at(community);   // Total weight of the community
    double weightNode = weight;                // Weight of the node
    double totalWeight = g.totalWeight;        // Total weight of the graph
    double weightNtC = weightNodeToCommunity;  // Weight of the node to the community

    return (weightNtC - totCommunity * weightNode / totalWeight);
}

double Community::modularity() const {
    double q = .0;
    double totalWeight = g.totalWeight;

    for (const auto& node : g.localNodes) {
        if (tot.at(node) <= 0) {
            continue;  // Skip empty communities
        }

        q += in.at(node) / totalWeight - (tot.at(node) / totalWeight) * (tot.at(node) / totalWeight);
    }
}

bool Community::step() {
    bool improvement = false;
    int moves;
    int stepsDone = 0;

    vector<unsigned int> randomOrder(size);
    for (unsigned int i = 0; i < size; i++)
        randomOrder[i] = i;
    for (unsigned int i = 0; i < size - 1; i++) {
        unsigned int position = rand() % (size - i) + i;
        unsigned int tmp = randomOrder[i];
        randomOrder[i] = randomOrder[position];
        randomOrder[position] = tmp;
    }

    moves = 0;
    stepsDone++;

    // Each node is removed from its community and inserted in the best one
    for (unsigned int nodeId = 0; nodeId < size; nodeId++) {
        unsigned int node = randomOrder[nodeId];
        unsigned int community = n2c[node];
        double weight = g.weightedDegree(node);

        neighbourCommunities(node);                            // Compute neighbour communities and their weights
        remove(node, community, neighbourWeights[community]);  // Remove node from its community

        // Compute the best community to insert the node
        // The default choice is the former community
        int bestCommunity = community;
        double bestLinks = 0.;
        double bestGain = 0.;
        for (unsigned int i = 0; i < neighbourLast; i++) {
            double gain = modularityGain(node, neighbourPositions[i], neighbourWeights[neighbourPositions[i]], weight);
            if (gain > bestGain) {
                bestCommunity = neighbourPositions[i];
                bestLinks = neighbourWeights[neighbourPositions[i]];
                bestGain = gain;
            }
        }

        // Insert the node into the best community found
        insert(node, bestCommunity, bestLinks);

        if (g.isRemote(bestCommunity)) {
            remoteCommunities.push_back(bestCommunity);  // Update remote communities mapping
            remoteWeights.push_back(bestLinks);          // Update remote weights
        }

        if (bestCommunity != community)
            moves++;
    }

    double totalTot = 0;  // Total incident weights to communities
    double totalIn = 0;   // Total internal weights for communities
    for (unsigned int i = 0; i < tot.size(); i++) {
        totalTot += tot[i];
        totalIn += in[i];
    }

    if (moves > 0)
        improvement = true;
}

// GraphBinary Community::graph() {
//     // Select non-empty communities
//     vector<int> renumber(size, -1);
//     for (int node = 0; node < size; node++) {
//         renumber[n2c[node]]++;
//     }

//     // Renumber communities
//     int final = 0;
//     for (int i = 0; i < size; i++)
//         if (renumber[i] != -1)
//             renumber[i] = final++;

//     // Compute communities
//     vector<vector<int>> communityNodes(final);
//     n2cNew.clear();
//     n2cNew.resize(size);
//     for (int node = 0; node < size; node++) {
//         communityNodes[renumber[n2c[node]]].push_back(node);
//         n2cNew[node] = renumber[n2c[node]];
//     }

//     // Compute weighted graph
//     GraphBinary g2;
//     g2.nNodes = communityNodes.size();
//     g2.degrees.resize(communityNodes.size());

//     // For each community, compute the total weight of the community
//     int communityDegree = communityNodes.size();
//     for (int community = 0; community < communityDegree; community++) {
//         map<int, double> m;
//         map<int, double>::iterator it;

//         // For each node in the community, find its neighbours and their weights
//         int communitySize = communityNodes[community].size();
//         for (int node = 0; node < communitySize; node++) {
//             pair<vector<unsigned int>::iterator, vector<double>::iterator> p = g.neighbours(communityNodes[community][node]);
//             int deg = g.nNeighbours(communityNodes[community][node]);
//             for (int i = 0; i < deg; i++) {
//                 int neighbour = *(p.first + i);
//                 int neighbourCommunity = renumber[n2c[neighbour]];
//                 double neighbourWeight = (g.weights.size() == 0) ? 1. : *(p.second + i);

//                 // If the node is the first, initialize the map
//                 // Otherwise, update the weight of the neighbour community
//                 it = m.find(neighbourCommunity);
//                 if (it == m.end())
//                     m.insert(make_pair(neighbourCommunity, neighbourWeight));
//                 else
//                     it->second += neighbourWeight;
//             }
//         }

//         g2.degrees[community] = (community == 0) ? m.size() : g2.degrees[community - 1] + m.size();
//         g2.nEdges += m.size();

//         // Add the community to the graph
//         for (it = m.begin(); it != m.end(); it++) {
//             g2.totalWeight += it->second;
//             g2.edges.push_back(it->first);
//             g2.weights.push_back(it->second);
//         }
//     }

//     return g2;
// }

// GraphBinary Community::graphNew() {
//     // Select non-empty communities
//     vector<int> renumber(size, -1);
//     for (int node = 0; node < size; node++) {
//         renumber[n2c[node]]++;
//     }

//     // Renumber communities
//     int final = 0;
//     for (int i = 0; i < size; i++)
//         if (renumber[i] != -1)
//             renumber[i] = final++;

//     // Compute communities
//     vector<vector<int>> communityNodes(final);
//     n2cNew.clear();
//     n2cNew.resize(size);
//     for (int node = 0; node < size; node++) {
//         communityNodes[renumber[n2c[node]]].push_back(node);
//         n2cNew[node] = renumber[n2c[node]];
//     }

//     // Compute weighted graph
//     GraphBinary g2;
//     g2.nNodes = communityNodes.size();
//     g2.degrees.resize(communityNodes.size());

//     // For each community, compute the total weight of the community
//     int communityDegree = communityNodes.size();
//     for (int community = 0; community < communityDegree; community++) {
//         map<int, float> m;
//         map<int, float>::iterator it;

//         // For each node in the community, find its neighbours and their weights
//         int comm_size = communityNodes[community].size();
//         for (int node = 0; node < comm_size; node++) {
//             pair<vector<unsigned int>::iterator, vector<double>::iterator> p = gE.neighbours(communityNodes[community][node]);
//             int deg = gE.nNeighbours(communityNodes[community][node]);
//             for (int i = 0; i < deg; i++) {
//                 int neigh = *(p.first + i);
//                 int neigh_comm = renumber[n2c[neigh]];
//                 double neigh_weight = (gE.weights->size == 0) ? 1. : *(p.second + i);

//                 // If the node is the first, initialize the map
//                 // Otherwise, update the weight of the neighbour community
//                 it = m.find(neigh_comm);
//                 if (it == m.end())
//                     m.insert(make_pair(neigh_comm, neigh_weight));
//                 else
//                     it->second += neigh_weight;
//             }

//             p = gE.remoteNeighbours(communityNodes[community][node]);
//             deg = gE.nRemoteNeighbours(communityNodes[community][node]);

//             for (int i = 0; i < deg; i++) {
//                 int neigh = *(p.first + i);
//                 int neigh_comm = renumber[n2c[neigh]];
//                 double neigh_weight = (gE.weights->size == 0) ? 1. : *(p.second + i);

//                 it = m.find(neigh_comm);
//                 if (it == m.end())
//                     m.insert(make_pair(neigh_comm, neigh_weight));
//                 else
//                     it->second += neigh_weight;
//             }
//         }
//         g2.degrees[community] = (community == 0) ? m.size() : g2.degrees[community - 1] + m.size();
//         g2.nEdges += m.size();

//         for (it = m.begin(); it != m.end(); it++) {
//             g2.totalWeight += it->second;
//             g2.edges.push_back(it->first);
//             g2.weights.push_back(it->second);
//         }
//     }

//     return g2;
// }
