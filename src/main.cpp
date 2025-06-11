
#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "community.h"
#include "graph.h"
#include "partitioner.h"

using namespace std;

char *fileName = NULL;

void display_time(const char *str) {
    time_t rawtime;
    time(&rawtime);
    cerr << str << ": " << ctime(&rawtime);
}

template <typename T>
void flattener(const vector<vector<T>> &vec, vector<T> &flat) {
    flat.resize(0);  // Clear the flat vector before inserting elements
    for (const auto &v : vec) {
        flat.insert(flat.end(), v.begin(), v.end());
    }
}

int main(int argc, char **argv) {
    if (argc < 2)
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
    else
        fileName = argv[1];

    time_t timeBegin, timeEnd;
    time(&timeBegin);

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ----------------------------------------------------------------
    // Distribute Graph
    // ----------------------------------------------------------------
    if (rank == 0) cerr << "Distributing graph..." << endl;

    Partitioner p(size);
    Graph localGraph;

    if (rank == 0) {
        display_time("Start");
        cerr << "Number of processes: " << size << endl;

        if (fileName == NULL) {
            cerr << "No input file specified." << endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        ifstream inputFile(fileName);
        if (!inputFile.is_open()) {
            cerr << "Could not open input file: " << fileName << endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        p.staticPartition(fileName);

        // Send partitioned data to all processes: nodes, cumulative degrees, edges, and weights
        for (unsigned int i = 1; i < size; i++) {
            int count = p.partitionNodes[i].size();
            MPI_Send(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD);                               // Send number of nodes
            MPI_Send(p.partitionNodes[i].data(), count, MPI_UNSIGNED, i, 1, MPI_COMM_WORLD);  // Send nodes

            count = p.partitionEdges[i].size();
            int bytes = count * sizeof(pair<unsigned int, unsigned int>);

            MPI_Send(&count, 1, MPI_INT, i, 2, MPI_COMM_WORLD);                                                     // Send number of edges
            MPI_Send(reinterpret_cast<void *>(p.partitionEdges[i].data()), bytes, MPI_BYTE, i, 3, MPI_COMM_WORLD);  // Send edges
        }

        // Initialize local graph with partitioned data
        localGraph.nNodes = p.partitionNodes[0].size();
        localGraph.nEdges = p.partitionEdges[0].size();
        localGraph.localNodes = set<unsigned int>(p.partitionNodes[0].begin(), p.partitionNodes[0].end());  // Set local nodes
        localGraph.edgeList = p.partitionEdges[0];
    } else {
        // Receive partitioned data from the root process
        int count;
        MPI_Status status;

        MPI_Recv(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);  // Receive number of nodes
        vector<unsigned int> localNodes(count);
        MPI_Recv(localNodes.data(), count, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, &status);  // Receive nodes
        localGraph.localNodes = set<unsigned int>(localNodes.begin(), localNodes.end());  // Set local nodes

        MPI_Recv(&count, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);  // Receive number of edges

        int bytes = count * sizeof(pair<unsigned int, unsigned int>);
        localGraph.nEdges = count;
        localGraph.edgeList.resize(count);

        MPI_Recv(localGraph.edgeList.data(), bytes, MPI_BYTE, 0, 3, MPI_COMM_WORLD, &status);  // Receive edges
    }

    // if (rank == 0) {
    //     cerr << "localGraph.edgeList:" << endl;
    //     for (const auto &edge : localGraph.edgeList) {
    //         cerr << "(" << edge.first << ", " << edge.second << ")" << endl;
    //     }
    // }

    // ----------------------------------------------------------------
    // Graph Initialization
    // ----------------------------------------------------------------
    if (rank == 0) cerr << "Initializing local graph..." << endl;
    localGraph.init();

    Community c(localGraph, -1, 0.0);  // Initialize community structure with local graph

    // ----------------------------------------------------------------
    // Collection of remote edges
    // ----------------------------------------------------------------

    vector<unsigned int> getList;         // List of remote nodes to retrieve
    vector<unsigned int> getProcessList;  // List of processes to get remote nodes from
    for (const auto &node : localGraph.localNodes) {
        vector<unsigned int> remoteNeighbours = localGraph.remoteNeighbours(node);

        for (const auto &remoteNode : remoteNeighbours) {
            if (find(getList.begin(), getList.end(), remoteNode) == getList.end()) {
                getList.push_back(remoteNode);  // Add remote node to the get list if not already present
            }
            if (find(getProcessList.begin(), getProcessList.end(), p.owner(remoteNode)) == getProcessList.end()) {
                getProcessList.push_back(p.owner(remoteNode));  // Add the owner process of the remote node to the get process list
            }
        }
    }

    vector<unsigned long> offsets;       // Offsets for each neighbour list
    vector<unsigned int> flatNodesList;  // Flattened list of neighbours
    vector<int> communities;             // Communities of the remote nodes

    MPI_Win windowOffsets;
    MPI_Win windowNodes;
    MPI_Win windowCommunities;

    offsets.push_back(0);  // Initialize the first offset to 0
    for (size_t i = 0; i < localGraph.neighboursList.size(); ++i) {
        const auto &vec = localGraph.neighboursList[i];
        flatNodesList.insert(flatNodesList.end(), vec.begin(), vec.end());  // Flatten the neighbour lists
        offsets.push_back(flatNodesList.size());                            // Update the offset for the next neighbour list
        communities.push_back(c.n2c[i]);                                    // Store the community of the local node
    }

    MPI_Win_create(flatNodesList.data(), flatNodesList.size() * sizeof(unsigned int), sizeof(unsigned int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowNodes);

    MPI_Win_create(offsets.data(), offsets.size() * sizeof(unsigned long), sizeof(unsigned long),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowOffsets);

    MPI_Win_create(communities.data(), communities.size() * sizeof(int), sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowCommunities);

    if (rank == 0) cerr << "Windows setup complete. Initiating remote edge collection..." << endl;

    // Get offsets of the neighbour lists from all processes

    vector<unsigned long> neighbourOffsets(2 * localGraph.neighboursList.size());  // Offsets for each node's neighbour list

    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowOffsets);      // Lock the window for all processes
    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowCommunities);  // Lock the communities window
    for (const auto &node : getList) {
        unsigned int owner = p.owner(node);  // Get the owner process of the remote node
        if (owner == rank) continue;         // Skip if the process is the current one

        MPI_Get(neighbourOffsets.data() + 2 * node, 2, MPI_UNSIGNED_LONG, owner,
                node, 2, MPI_UNSIGNED_LONG, windowOffsets);  // Get offsets for the remote node

        MPI_Get(&communities[node], 1, MPI_INT, owner,
                node, 1, MPI_INT, windowCommunities);  // Get the community of the remote node
    }
    MPI_Win_unlock_all(windowOffsets);      // Unlock the window after getting offsets
    MPI_Win_unlock_all(windowCommunities);  // Unlock the communities window

    if (rank == 0) cerr << "Offsets collected. Proceeding to get remote neighbours..." << endl;

    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowNodes);  // Lock the window for all processes
    for (const auto &node : getList) {
        unsigned int owner = p.owner(node);  // Get the owner process of the remote node
        if (owner == rank) continue;         // Skip if the process is the current one

        unsigned long offset = neighbourOffsets[2 * node];          // Get the starting byte for the neighbour list
        unsigned long offsetNext = neighbourOffsets[2 * node + 1];  // Get the next starting byte for the neighbour list
        unsigned long size = offsetNext - offset;                   // Calculate the size of the neighbour list

        if (size == 0) continue;  // Skip if the size is zero

        if (localGraph.neighboursList.size() <= node) {
            localGraph.neighboursList.resize(node + 1);  // Resize the neighbours list if needed
        }

        localGraph.neighboursList[node].resize(size);  // Resize the neighbour list for the remote node
        MPI_Get(localGraph.neighboursList[node].data(), size, MPI_UNSIGNED,
                owner, offset, size, MPI_UNSIGNED, windowNodes);  // Get the neighbour list for the remote node
    }
    MPI_Win_unlock_all(windowNodes);  // Unlock the window after getting remote neighbours

    offsets.resize(0);        // Clear the offsets vector as it is no longer needed
    flatNodesList.resize(0);  // Clear the flattened list of neighbours
    communities.resize(0);    // Clear the communities vector

    if (rank == 0) cerr << "Remote neighbours collected. Proceeding to update communities..." << endl;

    // Update the local graph with remote nodes and their communities
    c.resize();  // Resize the community structure based on the local graph
    for (const auto &node : getList) {
        unsigned int degree = localGraph.remoteDegree(node);  // Get the degree of the remote node
        localGraph.totalWeight += degree;                     // Update the total weight of the graph
        c.updateRemote(node, communities[node], degree);      // Update the community structure with the remote node
    }

    MPI_Win_free(&windowNodes);        // Free the window after use
    MPI_Win_free(&windowOffsets);      // Free the offsets window
    MPI_Win_free(&windowCommunities);  // Free the communities window

    // ----------------------------------------------------------------
    // Community Detection
    // ----------------------------------------------------------------
    if (rank == 0) cerr << "Starting community detection..." << endl;

    if (rank == 0) cerr << "Modularity before steps: " << c.modularity() << endl;

    c.step();

    // if (rank == 0) c.print();

    if (rank == 0) cerr << "Modularity: " << c.modularity() << endl;

    MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes before exchanging communities

    // ----------------------------------------------------------------
    // Exchange updated communities
    // ----------------------------------------------------------------
    if (rank == 0) cerr << "Exchanging updated communities..." << endl;

    vector<vector<unsigned int>> exchangingNodes(size);  // Nodes to exchange for each process
    vector<vector<int>> exchangingCommunities(size);     // Communities to exchange for each process
    vector<vector<int>> exchangingDegreeN2C(size);       // Weights to communities for each process
    vector<vector<int>> exchangingDegree(size);          // Totals for each process' communities
    vector<vector<int>> exchangingSelfLoops(size);       // Self-loops for each process' communities

    offsets.resize(0);        // Offsets for each process' community list
    flatNodesList.resize(0);  // Flattened list of nodes to get communities for

    vector<int> flatCommunitiesList;  // Flattened list of communities for the nodes
    vector<int> flatDegreeN2CList;    // Flattened list of weights to communities for the nodes
    vector<int> flatDegreeList;       // Flattened list of totals for each process' communities
    vector<int> flatSelfLoopsList;    // Flattened list of self-loops for each process' communities

    for (const auto &node : localGraph.localNodes) {
        unsigned int community = c.n2c[node];  // Get the community of the local node

        unsigned int ownerProcess = p.owner(community);  // Get the owner process of the community
        if (ownerProcess == rank) continue;              // Skip if the process is the current one

        exchangingNodes[ownerProcess].push_back(node);                   // Add the node to the remote communities list
        exchangingCommunities[ownerProcess].push_back(community);        // Add the community to the remote communities list
        exchangingDegreeN2C[ownerProcess].push_back(c.degreeN2C(node));  // Add the degree of the node to its community
        exchangingDegree[ownerProcess].push_back(c.tot[community]);      // Add the total weight of the community
        exchangingSelfLoops[ownerProcess].push_back(c.in[community]);    // Add the internal weight of the community
    }

    // Prepare the offsets and flattened list for sending communities
    offsets.push_back(0);  // Initialize the first offset to 0
    for (unsigned int i = 0; i < exchangingNodes.size(); ++i) {
        unsigned long size = exchangingNodes[i].size();
        offsets.push_back(size + offsets.back());                                                         // Update the offset for the next community list
        flatNodesList.insert(flatNodesList.end(), exchangingNodes[i].begin(), exchangingNodes[i].end());  // Flatten the community list
    }

    flattener(exchangingCommunities, flatCommunitiesList);  // Flatten the communities list
    flattener(exchangingDegreeN2C, flatDegreeN2CList);      // Flatten the weights to communities list
    flattener(exchangingDegree, flatDegreeList);            // Flatten the totals for each process' communities
    flattener(exchangingSelfLoops, flatSelfLoopsList);      // Flatten the self-loops for each process' communities

    MPI_Win windowDegreeN2C;
    MPI_Win windowDegree;
    MPI_Win windowSelfLoops;

    // Create windows for the offsets and flattened list
    MPI_Win_create(offsets.data(), offsets.size() * sizeof(unsigned long), sizeof(unsigned long),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowOffsets);

    MPI_Win_create(flatNodesList.data(), flatNodesList.size() * sizeof(unsigned int), sizeof(unsigned int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowNodes);

    MPI_Win_create(flatCommunitiesList.data(), flatCommunitiesList.size() * sizeof(int), sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowCommunities);

    MPI_Win_create(flatDegreeN2CList.data(), flatDegreeN2CList.size() * sizeof(int), sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowDegreeN2C);

    MPI_Win_create(flatDegreeList.data(), flatDegreeList.size() * sizeof(int), sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowDegree);

    MPI_Win_create(flatSelfLoopsList.data(), flatSelfLoopsList.size() * sizeof(int), sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowSelfLoops);

    if (rank == 0) cerr << "Windows for communities created. Proceeding to exchange communities..." << endl;

    // Get offsets of the community lists from all processes
    vector<unsigned long> processOffsets(size * 2);  // Offsets for each process' community list

    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowOffsets);  // Lock the window for all processes
    for (unsigned int i = 0; i < size; ++i) {
        if (i == rank) continue;  // Skip if the process is the current one

        MPI_Get(processOffsets.data() + 2 * i, 2, MPI_UNSIGNED_LONG, i,
                0, 2, MPI_UNSIGNED_LONG, windowOffsets);  // Get offsets for the remote communities
    }
    MPI_Win_unlock_all(windowOffsets);  // Unlock the window after getting offsets

    if (rank == 0) cerr << "Offsets collected. Proceeding to get remote neighbours..." << endl;

    vector<unsigned int> exchangedNodes;  // List of received nodes
    vector<int> exchangedCommunities;     // List of received communities
    vector<int> exchangedDegreeN2C;       // List of received weights to communities
    vector<int> exchangedDegree;          // List of received totals for each process' communities
    vector<int> exchangedSelfLoops;       // List of received self-loops for each process' communities

    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowNodes);        // Lock the window for all processes
    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowCommunities);  // Lock the communities window
    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowDegreeN2C);    // Lock the weights to communities window
    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowDegree);       // Lock the totals for each process' communities window
    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowSelfLoops);    // Lock the self-loops window

    for (unsigned int i = 0; i < size; ++i) {
        if (i == rank) continue;  // Skip if the process is the current one

        unsigned long offset = processOffsets[2 * i];          // Get the starting byte for the community list
        unsigned long offsetNext = processOffsets[2 * i + 1];  // Get the next starting byte for the community list
        unsigned long size = offsetNext - offset;              // Calculate the size of the community list

        if (size == 0) continue;  // Skip if the size is zero

        exchangedNodes.resize(exchangedNodes.size() + size);              // Resize the exchanged nodes list
        exchangedCommunities.resize(exchangedCommunities.size() + size);  // Resize the exchanged communities list
        exchangedDegreeN2C.resize(exchangedDegreeN2C.size() + size);      // Resize the exchanged weights to communities list
        exchangedDegree.resize(exchangedDegree.size() + size);            // Resize the exchanged totals for each process' communities list
        exchangedSelfLoops.resize(exchangedSelfLoops.size() + size);      // Resize the exchanged self-loops list

        MPI_Get(exchangedNodes.data() + exchangedNodes.size() - size, size, MPI_UNSIGNED,
                i, offset, size, MPI_UNSIGNED, windowNodes);  // Get the community list for the remote nodes

        MPI_Get(exchangedCommunities.data() + exchangedCommunities.size() - size, size, MPI_INT,
                i, offset, size, MPI_INT, windowCommunities);  // Get the communities for the remote nodes

        MPI_Get(exchangedDegreeN2C.data() + exchangedDegreeN2C.size() - size, size, MPI_INT,
                i, offset, size, MPI_INT, windowDegreeN2C);  // Get the weights to communities for the remote nodes

        MPI_Get(exchangedDegree.data() + exchangedDegree.size() - size, size, MPI_INT,
                i, offset, size, MPI_INT, windowDegree);  // Get the totals for each process' communities

        MPI_Get(exchangedSelfLoops.data() + exchangedSelfLoops.size() - size, size, MPI_INT,
                i, offset, size, MPI_INT, windowSelfLoops);  // Get the self-loops for each process' communities
    }

    MPI_Win_unlock_all(windowNodes);        // Unlock the window after getting remote communities
    MPI_Win_unlock_all(windowCommunities);  // Unlock the communities window
    MPI_Win_unlock_all(windowDegreeN2C);    // Unlock the weights to communities window
    MPI_Win_unlock_all(windowDegree);       // Unlock the totals for each process' communities window
    MPI_Win_unlock_all(windowSelfLoops);    // Unlock the self-loops window

    if (rank == 0) cerr << "Communities exchanged. Updating local communities..." << endl;

    // TODO: inserire il nodo anche in localnodes se non presente
    // Ma poi serve la neighbout list?
    for (size_t i = 0; i < exchangedNodes.size(); ++i) {
        unsigned int node = exchangedNodes[i];             // Get the global id of the node
        unsigned int community = exchangedCommunities[i];  // Get the global id of the community
        double wtc = exchangedDegreeN2C[i];                // Get the weight of the node to community connection
        double selfLoops = exchangedSelfLoops[i];          // Get the self-loops of the node
        double weightedDegree = exchangedDegree[i];        // Get the weighted degree of the node

        c.insert(node, community, wtc, weightedDegree, selfLoops);  // Insert the node into its community
        // if (localGraph.localNodes.find(node) == localGraph.localNodes.end()) {
        //     localGraph.localNodes.insert(node);        // Add the node to the local nodes if not already present
        //     localGraph.totalWeight += weightedDegree;  // Update the total weight of the graph
        // }
    }

    MPI_Win_free(&windowOffsets);      // Free the offsets window
    MPI_Win_free(&windowNodes);        // Free the nodes window
    MPI_Win_free(&windowCommunities);  // Free the communities window
    MPI_Win_free(&windowDegreeN2C);    // Free the weights to communities window
    MPI_Win_free(&windowDegree);       // Free the totals for each process' communities window
    MPI_Win_free(&windowSelfLoops);    // Free the self-loops window

    if (rank == 0) cerr << "Local communities updated. Community detection completed." << endl;

    // ----------------------------------------------------------------
    // Resolve duality conflicts
    // ----------------------------------------------------------------

    // if (rank == 0) cerr << "Resolving duality conflicts..." << endl;

    // MPI_Win_create(c.n2c.data(), c.n2c.size() * sizeof(int), sizeof(int),
    //                MPI_INFO_NULL, MPI_COMM_WORLD, &windowCommunities);  // Create a window for communities

    // MPI_Win_lock_all(MPI_MODE_NOCHECK, windowCommunities);  // Lock the communities window for all processes

    // for (auto node : localGraph.localNodes) {
    //     unsigned int comm = c.n2c[node];

    //     if (comm < node) {
    //         int comm_of_comm;
    //         unsigned int owner = p.owner(comm);

    //         if (owner == rank) {
    //             comm_of_comm = c.n2c[comm];
    //         } else {
    //             MPI_Get(&comm_of_comm, 1, MPI_INT, owner, comm, 1, MPI_INT, windowCommunities);
    //             MPI_Win_flush(owner, windowCommunities);
    //         }

    //         if (comm_of_comm == node) {
    //             int newComm = max(node, comm);

    //             if (node < comm) {
    //                 c.n2c[node] = newComm;
    //             }
    //         }
    //     }
    // }

    // MPI_Win_unlock_all(windowCommunities);
    // MPI_Win_free(&windowCommunities);

    // // Fase di normalizzazione opzionale: path compression
    // for (auto node : localGraph.localNodes) {
    //     while (c.n2c[node] != c.n2c[c.n2c[node]]) {
    //         c.n2c[node] = c.n2c[c.n2c[node]];
    //     }
    // }

    MPI_Barrier(MPI_COMM_WORLD);  // Sincronizzazione globale opzionale

    // for (const auto &node : localGraph.localNodes) {
    //     int community = c.n2c[node];  // Get the community of the local node
    //     if (community < node) {
    //         unsigned int owner = p.owner(community);  // Get the owner process of the community
    //         if (owner != rank) {
    //             // If the community is owned by another process, send the node to that process
    //             MPI_Get(&c.n2c[node], 1, MPI_INT, owner, node, 1, MPI_INT, windowCommunities);  // Get the community of the node
    //             MPI_Win_flush(owner, windowCommunities);  // Flush the window to ensure the update is visible to the owner process

    //         } else {
    //             // If the community is owned by the current process, update it locally
    //             c.n2c[node] = community;  // Update the community of the local node
    //         }
    //     }
    // }
    // MPI_Win_unlock_all(windowCommunities);  // Unlock the communities window after resolving conflicts

    // MPI_Win_free(&windowCommunities);  // Free the communities window

    // ----------------------------------------------------------------
    // Gather final community structure
    // ----------------------------------------------------------------

    MPI_Finalize();
    return 0;
}