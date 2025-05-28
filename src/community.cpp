#include "community.h"

Community::Community(string fileName) {
    size = g.edges.size();
    n2c.resize(size, -1);
    in.resize(size, 0.0);
    tot.resize(size, 0.0);

    // Initialize communities
    for (int i = 0; i < size; i++) {
        n2c[i] = i;                    // Each node starts in its own community
        tot[i] = g.weightedDegree(i);  // Total weight of edges incident to node i
        in[i] = g.selfLoops(i);        // Internal weight for community i
    }
}