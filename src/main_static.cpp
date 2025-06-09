
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
#include "graph_binary.h"
#include "partitioner.h"

using namespace std;

struct GraphData {
    unsigned int nNodes;
    unsigned long nEdges;
    double totalWeight;

    vector<unsigned long> degrees;
    vector<unsigned int> edges;
    vector<double> weights;
    vector<int> n2c;

    vector<int> rSource;
    vector<int> rDestination;
    vector<int> rPartition;
};

char *fileName = NULL;

void display_time(const char *str) {
    time_t rawtime;
    time(&rawtime);
    cerr << str << ": " << ctime(&rawtime);
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

    // TODO: Set number of threads for OpenMP
    // omp_set_num_threads(num_threads);

    // ----------------------------------------------------------------
    // Distribute Graph
    // ----------------------------------------------------------------
    if (rank == 0) cerr << "Distributing graph..." << endl;

    Partitioner p(size);
    GraphBinary localGraph;

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
            MPI_Send(p.partitionWeights[i].data(), count, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);                        // Send weights
        }

        // Initialize local graph with partitioned data
        localGraph.nNodes = p.partitionNodes[0].size();
        localGraph.nEdges = p.partitionEdges[0].size();
        localGraph.localNodes = p.partitionNodes[0];
        localGraph.edgeList = p.partitionEdges[0];
        localGraph.weightList = p.partitionWeights[0];
    } else {
        // Receive partitioned data from the root process
        int count;
        MPI_Status status;

        MPI_Recv(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);  // Receive number of nodes
        vector<unsigned int> localNodes(count);
        MPI_Recv(localNodes.data(), count, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, &status);  // Receive nodes
        localGraph.localNodes = localNodes;                                               // Set local nodes

        MPI_Recv(&count, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);  // Receive number of edges

        int bytes = count * sizeof(pair<unsigned int, unsigned int>);
        localGraph.nEdges = count;
        localGraph.edgeList.resize(count);
        localGraph.weightList.resize(count);

        MPI_Recv(localGraph.edgeList.data(), bytes, MPI_BYTE, 0, 3, MPI_COMM_WORLD, &status);      // Receive edges
        MPI_Recv(localGraph.weightList.data(), count, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &status);  // Receive weights
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
    cerr << "Local graph initialized with " << localGraph.nNodes << " nodes and " << localGraph.nEdges << " edges." << endl;

    Community c(localGraph, -1, 0.0);  // Initialize community structure with local graph

    // ----------------------------------------------------------------
    // Collection of remote edges
    // ----------------------------------------------------------------
    if (rank == 0) cerr << "Collecting remote edges..." << endl;

    vector<unsigned int> getList;         // List of remote nodes to retrieve
    vector<unsigned int> getProcessList;  // List of processes to get remote nodes from
    vector<unsigned int> shareList;       // List of nodes to share
    for (const auto &vertex : localGraph.neighboursList) {
        unsigned int node = vertex.first;
        vector<unsigned int> remoteNeighbours = localGraph.remoteNeighbours(node);

        for (const auto &remoteNode : remoteNeighbours) {
            if (find(getList.begin(), getList.end(), remoteNode) == getList.end()) {
                getList.push_back(remoteNode);  // Add remote node to the get list if not already present
            }
            if (find(getProcessList.begin(), getProcessList.end(), p.owner(remoteNode)) == getProcessList.end()) {
                getProcessList.push_back(p.owner(remoteNode));  // Add the owner process of the remote node to the get process list
            }
            if (find(shareList.begin(), shareList.end(), node) == shareList.end()) {
                shareList.push_back(node);  // Add the local node to the share list if not already present
            }
        }
    }

    /*
    The window will be set up as follows:
    - The first part contains the size of the share list (unsigned).
    - The second part contains a map of (id, startingByte) for each shared vertex.
    - The third part contains (community, degree, adjacency list) for each shared vertex.
    */

    // Setting up the windows
    unsigned long mapSize = (sizeof(unsigned int) + sizeof(unsigned long)) * shareList.size();
    unsigned long listSize = 0;
    for (const auto &vertex : shareList) {
        listSize += sizeof(unsigned int) + sizeof(double);                                     // community, degree
        listSize += localGraph.nNeighbours(vertex) * (sizeof(unsigned int) + sizeof(double));  // adjacency list
    }

    unsigned long windowSize = sizeof(unsigned int) +  // The number of shared verteices
                               mapSize +               // The map with (vertex, startingByte) for each shared vertex
                               listSize;               // The list with (community, degree, adjacency list) for each shared vertex

// size ||  id1: starting, id: starting.... || n1, n2, n3

    uint8_t *windowBuffer = new uint8_t[windowSize];

    unsigned int shareListSize = shareList.size();
    memcpy(windowBuffer, &shareListSize, sizeof(unsigned int));  // Store the size of the share list at the beginning of the buffer

    unsigned long lastByte = sizeof(unsigned int) + mapSize;  // Last byte written in the windowBuffer
    for (unsigned int i = 0; i < shareListSize; i++) {
        unsigned int vertex = shareList[i];
        memcpy(windowBuffer + sizeof(unsigned int) + i * (sizeof(unsigned int) + sizeof(unsigned long)),
               &vertex, sizeof(unsigned int));  // Store the vertex id in the map
        unsigned long startingByte = lastByte;  // Starting byte for the vertex in the windowBuffer
        memcpy(windowBuffer + sizeof(unsigned int) + i * (sizeof(unsigned int) + sizeof(unsigned long)) + sizeof(unsigned int),
               &startingByte, sizeof(unsigned long));  // Store the starting byte for the vertex

        unsigned int localId = vertex;  // Get the local id of the vertex

        unsigned int numberOfNeighbours = localGraph.nNeighbours(localId);  // Get the number of neighbours for the vertex
        unsigned int community = c.n2c[localId];                            // Get the community of the vertex
        double degree = localGraph.weightedDegree(localId);                 // Get the weighted degree of the vertex

        memcpy(windowBuffer + lastByte, &community, sizeof(unsigned int));  // Store the community id
        lastByte += sizeof(unsigned int);
        memcpy(windowBuffer + lastByte, &degree, sizeof(double));  // Store the weighted degree
        lastByte += sizeof(double);

        // Store the adjacency list for the vertex
        for (const auto &neighbour : localGraph.neighboursList[localId]) {
            unsigned int neighbourId = neighbour.first;  // Get the global id of the neighbour
            double weight = neighbour.second;            // Get the weight of the edge

            memcpy(windowBuffer + lastByte, &neighbourId, sizeof(unsigned int));  // Store the neighbour id
            lastByte += sizeof(unsigned int);
            memcpy(windowBuffer + lastByte, &weight, sizeof(double));  // Store the weight of the edge
            lastByte += sizeof(double);
        }
    }

    if (rank == 0) cerr << "Windows setup complete. Initiating remote edge collection..." << endl;

    // Share the window size with all processes
    vector<unsigned long> windowSizes(size);
    unsigned long localWindowSize = windowSize;

    MPI_Allgather(&localWindowSize, 1, MPI_UNSIGNED_LONG, windowSizes.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    MPI_Win window;
    MPI_Win_create(windowBuffer, windowSize, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &window);  // Create a window for shared memory

    MPI_Win_fence(0, window);  // Synchronize all processes before accessing the window

    // Gather remote edges from needed processes
    vector<uint8_t *> remoteBuffers(size);
    for (const auto &process : getProcessList) {
        if (process == rank) continue;  // Skip if the process is the current one

        unsigned long remoteWindowSize = windowSizes[process];
        remoteBuffers[process] = new uint8_t[remoteWindowSize];  // Allocate buffer for remote data

        MPI_Get(remoteBuffers[process], remoteWindowSize, MPI_BYTE, process, 0, remoteWindowSize, MPI_BYTE, window);  // Get remote data
    }

    MPI_Win_fence(0, window);  // Synchronize all processes after accessing the window

    // Update the local graph with remote edges
    for (const auto &process : getProcessList) {
        if (process == rank) continue;  // Skip if the process is the current one

        unsigned long remoteWindowSize = windowSizes[process];
        unsigned int remoteShareListSize = *reinterpret_cast<unsigned int *>(remoteBuffers[process]);
        unsigned long offset = sizeof(unsigned int);  // Offset to the start of the map

        for (unsigned int i = 0; i < remoteShareListSize; i++) {
            unsigned int vertex = *reinterpret_cast<unsigned int *>(remoteBuffers[process] + offset);
            if (find(getList.begin(), getList.end(), vertex) == getList.end()) {
                offset += sizeof(unsigned int) + sizeof(unsigned long);  // Skip if the vertex is not in the get list
                continue;
            }
            offset += sizeof(unsigned int);  // Move to the next part of the map
            unsigned long startingByte = *reinterpret_cast<unsigned long *>(remoteBuffers[process] + offset);
            offset += sizeof(unsigned long);  // Move to the next part of the map

            unsigned long nextStartingByte = (i + 1 < remoteShareListSize) ? *reinterpret_cast<unsigned long *>(remoteBuffers[process] + offset + sizeof(unsigned int)) : remoteWindowSize;

            unsigned int community = *reinterpret_cast<unsigned int *>(remoteBuffers[process] + startingByte);
            startingByte += sizeof(unsigned int);  // Move to the degree
            double degree = *reinterpret_cast<double *>(remoteBuffers[process] + startingByte);
            startingByte += sizeof(double);  // Move to the adjacency list

            unsigned long adjacencyListSize = (nextStartingByte - startingByte) / (sizeof(unsigned int) + sizeof(double));  // Calculate the size of the adjacency list

            // Add the remote edges to the local graph
            for (unsigned long j = 0; j < adjacencyListSize; j++) {
                unsigned int neighbourId = *reinterpret_cast<unsigned int *>(remoteBuffers[process] + startingByte);
                startingByte += sizeof(unsigned int);  // Move to the weight
                double weight = *reinterpret_cast<double *>(remoteBuffers[process] + startingByte);
                startingByte += sizeof(double);  // Move to the next neighbour

                localGraph.addEdge(vertex, neighbourId, weight);  // Add the edge to the local graph
            }

            c.updateRemote(vertex, community, degree);  // Update the community structure with the remote node
        }

        delete[] remoteBuffers[process];  // Deallocate remote buffer
    }

    MPI_Win_free(&window);  // Free the window after use
    delete windowBuffer;    // Deallocate the window buffer

    // ----------------------------------------------------------------
    // Community Detection
    // ----------------------------------------------------------------
    if (rank == 0) cerr << "Starting community detection..." << endl;

    // cerr << "Rank: " << rank << " Total weight of the graph: " << localGraph.totalWeight << endl;
    if (rank == 0) cerr << "Modularity before steps: " << c.modularity() << endl;

    // for (const auto &node : localGraph.localNodes) {
    //     cerr << "Node: " << node
    //          << ", Community: " << c.n2c[node]
    //          << ", Weighted Degree: " << localGraph.weightedDegree(node)
    //          << ", Self Loops: " << localGraph.selfLoops(node) << endl;
    // }

    c.step();

    // for (const auto &node : localGraph.localNodes) {
    //     cerr << "Node: " << node
    //          << ", Community: " << c.n2c[node]
    //          << ", Weighted Degree: " << localGraph.weightedDegree(node)
    //          << ", Self Loops: " << localGraph.selfLoops(node) << endl;
    // }

    // if (rank == 0) c.print();

    if (rank == 0) cerr << "Modularity: " << c.modularity() << endl;

    MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes before exchanging communities

    // ----------------------------------------------------------------
    // Exchange updated communities
    // ----------------------------------------------------------------
    // if (rank == 0) cerr << "Exchanging updated communities..." << endl;

    // MPI_Barrier(MPI_COMM_WORLD);                            // Synchronize all processes before exchanging communities
    // vector<vector<unsigned int>> sendVertexIds(size);       // List of vertex ids to send to each process
    // vector<vector<unsigned int>> sendCommunities(size);     // List of communities to send to each process
    // vector<vector<double>> sendWeightsToCommunities(size);  // List of weights to communities to send to each process
    // vector<vector<double>> sendSelfLoops(size);             // List of self-loops to send to each process
    // vector<vector<double>> sendWeights(size);               // List of weighted degrees to send to each process

    // for (unsigned int i = 0; i < c.remoteCommunities.size(); ++i) {
    //     unsigned int node = c.remoteCommunities[i];               // Get the global id of the node
    //     unsigned int community = c.n2c[node];                     // Get the global id of the community
    //     double wtc = c.remoteWeights[i];                          // Get the weight of the node to community connection
    //     double selfLoops = localGraph.selfLoops(node);            // Get the self-loops of the node
    //     double weightedDegree = localGraph.weightedDegree(node);  // Get the weighted degree of the node
    //     int ownerProcess = p.owner(community);                // Get the owner process of the community

    //     sendVertexIds[ownerProcess].push_back(node);            // Add the node to the list of vertex ids to send
    //     sendCommunities[ownerProcess].push_back(community);     // Add the community to the list of communities to send
    //     sendWeightsToCommunities[ownerProcess].push_back(wtc);  // Add the weight to the list of weights to communities to send
    //     sendSelfLoops[ownerProcess].push_back(selfLoops);       // Add the self-loops to the list of self-loops to send
    //     sendWeights[ownerProcess].push_back(weightedDegree);    // Add the weighted degree to the list of weights to send
    // }

    // if (rank == 0) cerr << "Communities collected. Sending updates to all processes..." << endl;

    // // Send the updated communities to each process
    // for (unsigned int i = 0; i < size; ++i) {
    //     if (i == rank) continue;  // Skip if the process is the current one

    //     unsigned int count = sendVertexIds[i].size();
    //     MPI_Send(&count, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD);  // Send the number of vertices to the process

    //     if (count > 0) {
    //         MPI_Send(sendVertexIds[i].data(), count, MPI_UNSIGNED, i, 1, MPI_COMM_WORLD);           // Send the vertex ids
    //         MPI_Send(sendCommunities[i].data(), count, MPI_UNSIGNED, i, 2, MPI_COMM_WORLD);         // Send the communities
    //         MPI_Send(sendWeightsToCommunities[i].data(), count, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);  // Send weights to communities
    //         MPI_Send(sendSelfLoops[i].data(), count, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);             // Send self-loops
    //         MPI_Send(sendWeights[i].data(), count, MPI_DOUBLE, i, 5, MPI_COMM_WORLD);               // Send weighted degrees
    //     }
    // }

    // if (rank == 0) cerr << "Updates sent. Receiving communities from all processes..." << endl;

    // vector<vector<unsigned int>> recvVertexIds(size);       // List of vertex ids received from each process
    // vector<vector<unsigned int>> recvCommunities(size);     // List of communities received from each process
    // vector<vector<double>> recvWeightsToCommunities(size);  // List of weights to communities received from each process
    // vector<vector<double>> recvSelfLoops(size);             // List of self-loops received from each process
    // vector<vector<double>> recvWeights(size);               // List of weighted degrees received from each process

    // vector<MPI_Request> requests;  // Requests for non-blocking communication
    // unsigned int cr = 0;           // Current request index

    // // Receive the updated communities from each process
    // for (unsigned int i = 0; i < size; ++i) {
    //     if (i == rank) continue;  // Skip if the process is the current one

    //     unsigned int count;
    //     MPI_Status status;
    //     MPI_Recv(&count, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD, &status);  // Receive the number of vertices

    //     if (count > 0) {
    //         recvVertexIds[i].resize(count);
    //         recvCommunities[i].resize(count);
    //         recvWeightsToCommunities[i].resize(count);
    //         recvSelfLoops[i].resize(count);
    //         recvWeights[i].resize(count);

    //         requests.resize(cr + 5);  // Resize requests vector to hold new requests

    //         MPI_Irecv(recvVertexIds[i].data(), count, MPI_UNSIGNED, i, 1, MPI_COMM_WORLD, requests.data() + cr++);           // Receive vertex ids
    //         MPI_Irecv(recvCommunities[i].data(), count, MPI_UNSIGNED, i, 2, MPI_COMM_WORLD, requests.data() + cr++);         // Receive communities
    //         MPI_Irecv(recvWeightsToCommunities[i].data(), count, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, requests.data() + cr++);  // Receive weights to communities
    //         MPI_Irecv(recvSelfLoops[i].data(), count, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, requests.data() + cr++);             // Receive self-loops
    //         MPI_Irecv(recvWeights[i].data(), count, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, requests.data() + cr++);               // Receive weighted degrees
    //     }
    // }

    // if (rank == 0) cerr << "Communities received. Waiting for all receive requests to complete..." << endl;

    // // Wait for all receive requests to complete
    // MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    // if (rank == 0) cerr << "All receive requests completed. Updating community structure..." << endl;

    // // Update the community structure with the received data
    // for (unsigned int i = 0; i < size; ++i) {
    //     if (i == rank) continue;  // Skip if the process is the current one

    //     for (unsigned int j = 0; j < recvVertexIds[i].size(); ++j) {
    //         unsigned int node = recvVertexIds[i][j];
    //         unsigned int community = recvCommunities[i][j];
    //         double wtc = recvWeightsToCommunities[i][j];
    //         double selfLoops = recvSelfLoops[i][j];
    //         double weightedDegree = recvWeights[i][j];

    //         c.insert(node, community, wtc, weightedDegree, selfLoops);  // Insert the node into the community structure
    //     }
    // }

    // ----------------------------------------------------------------
    // Gather final community structure
    // ----------------------------------------------------------------

    cerr << "Gathering final community structure..." << endl;

    vector<unsigned int> communitiesToSend;
    for (const auto &node : localGraph.localNodes) {
        communitiesToSend.push_back(node);
        communitiesToSend.push_back(c.n2c[node]);  // Add the node and its community to the list
    }

    // Gather sizes of communities to send
    int sizeOfCommunitiesToSend = communitiesToSend.size();
    vector<int> allSizes(size);
    MPI_Gather(&sizeOfCommunitiesToSend, 1, MPI_UNSIGNED, allSizes.data(), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    vector<int> displacements(size, 0);
    int totalSize = 0;
    if (rank == 0) {
        for (int i = 0; i < size; ++i) {
            displacements[i] = totalSize;  // Calculate displacements for each process
            totalSize += allSizes[i];      // Update total size with the size of the current process
        }
    }

    vector<unsigned int> allCommunities;
    if (rank == 0) {
        allCommunities.resize(totalSize);  // Resize the vector to hold all communities
    }

    // Gather communities from all processes
    MPI_Gatherv(communitiesToSend.data(), sizeOfCommunitiesToSend, MPI_UNSIGNED,
                allCommunities.data(), allSizes.data(), displacements.data(), MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (rank == 0) cerr << "Community structure gathered from all processes." << endl;

    MPI_Finalize();
    return 0;
}