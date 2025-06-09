#include "community.h"

Community::Community(GraphBinary& gb, int st, double thr) : g(gb) {
    size = gb.nNodes;

    for (const auto& node : g.localNodes) {
        n2c[node] = node;              // Initialize each node to its own community
        tot[node] = g.degree(node);    // Total weight of the community is the weighted degree of the node
        in[node] = g.selfLoops(node);  // Internal weight of the community is the self-loops of the node
    }

    steps = st;
    threshold = thr;
}

void Community::updateRemote(unsigned int node, unsigned int community, double degree) {
    // Update the remote community structure
    n2c[node] = community;  // Assign the node to the specified community
    tot[community] += degree;
    in[community] += g.remoteSelfLoops(node);
}

void Community::print() {
    cout << "Community structure: (node community)" << endl;
    set<unsigned int> communities;  // Set to store unique communities
    for (const auto& node : g.localNodes) {
        cout << node << " " << n2c[node] << endl;  // Print each node and its community
        communities.insert(n2c[node]);             // Insert the community into the set
    }
    cout << "Total communities: " << communities.size() << endl;  // Print the total number of unique communities
}

void Community::insert(int node, int community, double weightNodeToCommunity) {
    if (g.isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    tot[community] += g.degree(node);
    in[community] += 2 * weightNodeToCommunity + g.selfLoops(node);
    n2c[node] = community;  // Assign the node to the specified community
}

void Community::insert(int node, int community, double weightNodeToCommunity, double weight, double selfLoops) {
    tot[community] += weight;                                // Update the total weight of the community
    in[community] += 2 * weightNodeToCommunity + selfLoops;  // Update the internal weight of the community
    n2c[node] = community;                                   // Assign the node to the specified community
}

void Community::remove(int node, int community, double weightNodeToCommunity) {
    if (g.isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    tot[community] -= g.weightedDegree(node);
    in[community] -= 2 * weightNodeToCommunity + g.selfLoops(node);
    n2c[node] = UINT_MAX;  // Remove the node from its community
}

void Community::neighbourCommunities(int node) {
    if (g.isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    neighbourCommunitiesMap.clear();  // Clear the neighbour weights map

    for (const auto& n : g.neighboursList[node]) {
        if (n.first == node) continue;

        unsigned int comm = n2c[n.first];           // Get the community of the neighbour
        neighbourCommunitiesMap[comm] += n.second;  // Add the weight of the edge to the community weight
    }
}

double Community::modularityGain(int community, double weightNodeToCommunity, double weight) {
    double totCommunity = 0;  // Total weight of the community

    totCommunity = tot[community];             // Total weight of the community
    double weightNode = weight;                // Weight of the node
    double totalWeight = g.totalWeight;        // Total weight of the graph
    double weightNtC = weightNodeToCommunity;  // Weight of the node to the community

    return (weightNtC - totCommunity * weightNode / totalWeight);
}

double Community::modularity() {
    double q = .0;
    double totalWeight = g.totalWeight;

    std::set<unsigned int> processed;
    for (const auto& pair : n2c) {
        unsigned int node = pair.first;
        unsigned int community = pair.second;

        if (processed.count(community)) continue;
        processed.insert(community);

        double totC = tot[community];
        double inC = in[community];

        // cerr << "Community " << community
        //      << ": tot = " << totC
        //      << ", in = " << inC << endl;

        q += inC / totalWeight - (totC / totalWeight) * (totC / totalWeight);
    }

    return q;
}

bool Community::step() {
    bool improvement = false;
    int moves;
    int stepsDone = 0;

    vector<unsigned int> randomOrder(g.localNodes.begin(), g.localNodes.end());
    random_device rd;
    mt19937 r(rd());
    shuffle(randomOrder.begin(), randomOrder.end(), r);  // Shuffle the nodes to randomize the order

    moves = 0;
    stepsDone++;

    // Each node is removed from its community and inserted in the best one
    for (unsigned int nodeId = 0; nodeId < randomOrder.size(); nodeId++) {
        unsigned int node = randomOrder[nodeId];  // Get the node in the random order
        unsigned int community = n2c[node];
        double weight = g.weightedDegree(node);

        neighbourCommunities(node);                                   // Compute neighbour communities and their weights
        remove(node, community, neighbourCommunitiesMap[community]);  // Remove node from its community

        // Compute the best community to insert the node
        // The default choice is the former community
        unsigned int bestCommunity = community;
        double bestLinks = 0.;
        double bestGain = 0.;

        // TODO: self loop se è da solo?
        for (const auto& n : neighbourCommunitiesMap) {
            double gain = modularityGain(n.first, n.second, weight);  // Compute the modularity gain for moving the node
            // cerr << "Node " << node << " to community " << n.first << " gain: " << gain << endl;
            if (gain > bestGain) {
                bestGain = gain;          // Update the best gain found
                bestCommunity = n.first;  // Update the best community
                bestLinks = n.second;     // Update the best links to the community
            }
        }

        insert(node, bestCommunity, bestLinks);  // Insert the node into the best community found

        if (g.isRemote(bestCommunity)) {
            remoteCommunities.push_back(node);   // Update remote communities mapping
            remoteWeights.push_back(bestLinks);  // Update remote weights
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

    return improvement;
}