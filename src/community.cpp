#include "community.h"

Community::Community(char *inFile, char *inWeightsFile, int st, double thr) {
    g = GraphBinary(inFile, inWeightsFile);
    size = g.nNodes;

    neighbourWeights.resize(size, -1);
    neighbourPositions.resize(size);
    neighbourLast = 0;

    n2c.resize(size);
    in.resize(size);
    tot.resize(size);

    for (int i = 0; i < size; i++) {
        n2c[i] = i;
        tot[i] = g.weightedDegree(i);
        in[i] = g.selfLoops(i);
    }

    steps = st;
    threshold = thr;
}

Community::Community(GraphBinary gb, int st, double thr) {
    g = gb;
    size = g.nNodes;

    neighbourWeights.resize(size, -1);
    neighbourPositions.resize(size);
    neighbourLast = 0;

    n2c.resize(size);
    in.resize(size);
    tot.resize(size);

    for (int i = 0; i < size; i++) {
        n2c[i] = i;
        tot[i] = g.weightedDegree(i);
        in[i] = g.selfLoops(i);
    }

    steps = st;
    threshold = thr;
}

void Community::print() const {
    for (int i = 0; i < size; i++)
        cerr << " " << i << "/" << n2c[i] << "/" << in[i] << "/" << tot[i];
    cerr << endl;
}

void Community::insert(int node, int community, double weightNodeToCommunity) {
    if (node >= size || node < 0) {
        throw out_of_range("Node index out of range");
    }

    tot[community] += g.weightedDegree(node);
    in[community] += 2 * weightNodeToCommunity + g.selfLoops(node);
    n2c[node] = community;
}

void Community::remove(int node, int community, double weightNodeToCommunity) {
    if (node >= size || node < 0) {
        throw out_of_range("Node index out of range");
    }

    tot[community] -= g.weightedDegree(node);
    in[community] -= 2 * weightNodeToCommunity + g.selfLoops(node);
    n2c[node] = -1;
}

double Community::modularityGain(int node, int community, double weightNodeToCommunity, double weight) const {
    if (node >= size || node < 0) {
        throw out_of_range("Node or community index out of range");
    }

    double totCommunity = tot[community];      // Total weight of the community
    double weightNode = weight;                // Weight of the node
    double totalWeight = g.totalWeight;        // Total weight of the graph
    double weightNtC = weightNodeToCommunity;  // Weight of the node to the community

    return (weightNtC - totCommunity * weightNode / totalWeight);
}

void Community::neighbourCommunities(int node) {
    for (unsigned int i = 0; i < neighbourLast; i++)
        neighbourWeights[neighbourPositions[i]] = -1;
    neighbourLast = 0;

    pair<vector<unsigned int>::iterator, vector<double>::iterator> p = g.neighbours(node);
    unsigned int deg = g.nNeighbours(node);

    neighbourPositions[0] = n2c[node];            // Ids of the neighbour community
    neighbourWeights[neighbourPositions[0]] = 0;  // Weights of the neighbour community
    neighbourLast = 1;                            // Number of neighbour communities

    for (unsigned int i = 0; i < deg; i++) {
        unsigned int neighbour = *(p.first + i);
        unsigned int neighbourCommunity = n2c[neighbour];
        double neighbourWeight = (g.weights.size() == 0) ? 1. : *(p.second + i);

        if (neighbour != node) {
            if (neighbourWeights[neighbourCommunity] == -1) {
                neighbourWeights[neighbourCommunity] = 0.;
                neighbourPositions[neighbourLast++] = neighbourCommunity;
            }
            neighbourWeights[neighbourCommunity] += neighbourWeight;
        }
    }
}

double Community::modularity() const {
    double q = .0;
    double totalWeight = g.totalWeight;

    for (unsigned int i = 0; i < size; i++)
        if (tot[i] > 0) q += in[i] / totalWeight - (tot[i] * tot[i]) / (totalWeight * totalWeight);

    return q;
}

bool Community::step() {
    bool improvement = false;
    int moves;
    int stepsDone = 0;
    double newModularity = modularity();
    double currentModularity = newModularity;

    vector<int> randomOrder(size);
    for (int i = 0; i < size; i++)
        randomOrder[i] = i;
    for (int i = 0; i < size - 1; i++) {
        int position = rand() % (size - i) + i;
        int tmp = randomOrder[i];
        randomOrder[i] = randomOrder[position];
        randomOrder[position] = tmp;
    }

    // repeat while
    //   there is an improvement of modularity
    //   or there is an improvement of modularity greater than a given epsilon
    //   or a predefined number of pass have been done
    do {
        currentModularity = newModularity;
        moves = 0;
        stepsDone++;

        // Each node is removed from its community and inserted in the best one
        for (int nodeId = 0; nodeId < size; nodeId++) {
            int node = randomOrder[nodeId];
            int community = n2c[node];
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

            if (bestCommunity != community)
                moves++;
        }

        double totalTot = 0;  // Total incident weights to communities
        double totalIn = 0;   // Total internal weights for communities
        for (unsigned int i = 0; i < tot.size(); i++) {
            totalTot += tot[i];
            totalIn += in[i];
        }

        newModularity = modularity();
        if (moves > 0)
            improvement = true;
    } while (moves > 0 && ((newModularity - currentModularity) > threshold));

    return improvement;
}