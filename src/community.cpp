#include "community.h"

Community::Community(Graph& gb, int st, double thr) : g(gb) {
    // TODO: All'inizio ogni nodo è parte della sua community, quindi non serve condividerlo
    size = g.neighboursList.size();  // Set the size of the community structure based on the graph

    n2c.resize(size, -1);   // Initialize node to community mapping with -1
    tot.resize(size, 0.0);  // Initialize total weight of each community to 0
    in.resize(size, 0.0);   // Initialize internal weight of each community to 0

    for (const auto& node : g.localNodes) {
        n2c[node] = node;              // Initialize each node to its own community
        tot[node] = g.degree(node);    // Total weight of the community is the weighted degree of the node
        in[node] = g.selfLoops(node);  // Internal weight of the community is the self-loops of the node
    }

    steps = st;
    threshold = thr;
}

void Community::resize() {
    size = g.neighboursList.size();  // Update the size of the community structure based on the graph

    if (size > n2c.size()) n2c.resize(size, -1);   // Resize node to community mapping
    if (size > tot.size()) tot.resize(size, 0.0);  // Resize total weight of each community
    if (size > in.size()) in.resize(size, 0.0);    // Resize internal weight of each community
}

void Community::updateRemote(int node, int community, int degree) {
    if (node >= g.neighboursList.size()) {
        throw out_of_range("Node index out of range");
    }

    if (community >= tot.size()) {
        throw out_of_range("Community index out of range");
    }

    // Update the remote community structure
    n2c[node] = community;  // Assign the node to the specified community
    tot[community] += degree;
    in[community] += g.selfLoops(node) + 2 * degreeN2C(node);  // Update the internal weight of the community
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

    tot[community] -= g.degree(node);
    in[community] -= 2 * weightNodeToCommunity + g.selfLoops(node);
    n2c[node] = -1;  // Remove the node from its community
}

void Community::neighbourCommunities(int node) {
    if (g.isRemote(node)) {
        throw out_of_range("Node index out of range");
    }

    neighbourCommunitiesMap.clear();  // Clear the neighbour weights map

    // Iterate through the neighbours of the node
    for (const auto& n : g.neighboursList[node]) {
        if (n == node) continue;  // Skip self-loops

        unsigned int comm = n2c[n];       // Get the community of the neighbour
        neighbourCommunitiesMap[comm]++;  // Increment the weight of the edge to the community
    }
}

int Community::degreeN2C(int node) {
    if (!g.isCollected(node)) {
        throw out_of_range("Node index out of range");
    }

    if (n2c[node] < 0) {
        return 0;  // Node is not assigned to any community
    }

    int degree = 0;             // Initialize the degree to 0
    int community = n2c[node];  // Get the community of the node
    for (const auto& n : g.neighboursList[node]) {
        if (n != node && n2c[n] == community) {  // Check if the neighbour is not the node itself and belongs to the same community
            degree++;                            // Increment the degree if it is
        }
    }
    return degree;
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
    double totalWeight = (double)g.totalWeight;

    std::set<unsigned int> processed;
    for (const auto& node : g.localNodes) {
        unsigned int community = n2c[node];  // Get the community of the node
        if (processed.count(community) > 0 || community < 0) continue;

        processed.insert(community);  // Mark the community as processed

        double totCommunity = tot[community];
        double inCommunity = in[community];

        // Compute the modularity contribution for the community
        q += (inCommunity / totalWeight) - (totCommunity * totCommunity / (totalWeight * totalWeight));
    }

    return q;
}

double Community::modularity(int community) {
    if (community < 0 || community >= (int)tot.size()) {
        throw out_of_range("Community index out of range");
    }

    double q = 0.0;
    double totalWeight = (double)g.totalWeight;  // Total weight of the graph
    double totCommunity = tot[community];        // Total weight of the community
    double inCommunity = in[community];          // Internal weight of the community

    // Compute the modularity contribution for the community
    q += (inCommunity / totalWeight) - (totCommunity * totCommunity / (totalWeight * totalWeight));

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
        unsigned int degree = g.degree(node);

        neighbourCommunities(node);                                   // Compute neighbour communities and their weights
        remove(node, community, neighbourCommunitiesMap[community]);  // Remove node from its community

        // Compute the best community to insert the node
        // The default choice is the former community
        unsigned int bestCommunity = community;
        double bestLinks = 0.;
        double bestGain = 0.;

        // TODO: self loop se è da solo?
        for (const auto& c : neighbourCommunitiesMap) {
            double gain = modularityGain(c.first, c.second, degree);  // Compute the modularity gain for moving the node

            if (gain > bestGain) {
                bestGain = gain;          // Update the best gain found
                bestCommunity = c.first;  // Update the best community
                bestLinks = c.second;     // Update the best links to the community
            }
        }

        insert(node, bestCommunity, bestLinks);  // Insert the node into the best community found

        if (bestCommunity != community)
            moves++;
    }

    if (moves > 0)
        improvement = true;

    return improvement;
}