
#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include "community.h"
#include "graph.h"
#include "partitioner.h"

using namespace std;

char *fileName = NULL;

void collectMissingNodes(Graph &localGraph, Community &c, Partitioner &p, int rank, vector<unsigned int> &getList, vector<unsigned int> &getProcessList);
void resolveDuality(vector<int> &n2c);
double computeModularity(set<unsigned int> &communities, Community &c, Partitioner &p, int rank);

int main(int argc, char **argv) {
    if (argc < 2)
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
    else
        fileName = argv[1];

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ----------------------------------------------------------------
    // Distribute Graph
    // ----------------------------------------------------------------
    double startTime = MPI_Wtime();

    Partitioner p(size);
    Graph localGraph;
    Graph coarseGraph;

    double modularity = 0.0;
    double newModularity = 0.0;
    unsigned int foundCommunities = UINT_MAX;  // Counter for found communities
    vector<unsigned int> communities;          // Vector to store communities

    if (rank == 0) {
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

            count = p.partitionMap.size();
            MPI_Send(&count, 1, MPI_INT, i, 4, MPI_COMM_WORLD);                          // Send size of partition map
            MPI_Send(p.partitionMap.data(), count, MPI_UNSIGNED, i, 5, MPI_COMM_WORLD);  // Send partition map

            int weight = p.edgeList.size();
            MPI_Send(&weight, 1, MPI_INT, i, 6, MPI_COMM_WORLD);  // Send number of edges in the graph
        }

        // Initialize local graph with partitioned data
        localGraph.nNodes = p.partitionNodes[0].size();
        localGraph.nEdges = p.partitionEdges[0].size();
        localGraph.localNodes = set<unsigned int>(p.partitionNodes[0].begin(), p.partitionNodes[0].end());  // Set local nodes
        localGraph.edgeList = p.partitionEdges[0];
        localGraph.totalWeight = p.edgeList.size();  // Set total weight based on the number of edges
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

        MPI_Recv(&count, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);  // Receive size of partition map
        p.partitionMap.resize(count);
        MPI_Recv(p.partitionMap.data(), p.partitionMap.size(), MPI_UNSIGNED, 0, 5, MPI_COMM_WORLD, &status);  // Receive partition map

        int weight;
        MPI_Recv(&weight, 1, MPI_INT, 0, 6, MPI_COMM_WORLD, &status);  // Receive number of edges in the graph
        localGraph.totalWeight = weight;                               // Set total weight based on the number of edges
    }

    double partitionTime = MPI_Wtime();  // Measure partitioning time

    localGraph.init();

    int step = 0;
    while (true) {
        // ----------------------------------------------------------------
        // Graph Initialization
        // ----------------------------------------------------------------
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

        collectMissingNodes(localGraph, c, p, rank, getList, getProcessList);  // Collect missing nodes from remote processes
        // ----------------------------------------------------------------
        // Community Detection
        // ----------------------------------------------------------------
        modularity = computeModularity(localGraph.localNodes, c, p, rank);  // Compute the initial modularity of the community structure

        c.step();

        // ----------------------------------------------------------------
        // Exchange updated communities
        // ----------------------------------------------------------------
        int nodesCount = 0;              // Count of local nodes
        vector<unsigned int> nodesList;  // Flattened list of nodes to get communities for
        vector<int> communitiesList;     // Communities of the remote nodes

        for (const auto &node : localGraph.localNodes) {
            unsigned int community = c.n2c[node];  // Get the community of the local node
            nodesList.push_back(node);             // Add the node to the flattened list
            communitiesList.push_back(community);  // Add the community to the communities list
        }
        nodesCount = nodesList.size();  // Get the count of local nodes

        // Create windows for the nodes and communities
        MPI_Win windowCount;        // Window for the count of local nodes
        MPI_Win windowNodes;        // Window for the list of local nodes
        MPI_Win windowCommunities;  // Window for the communities of the local nodes

        MPI_Win_create(&nodesCount, sizeof(int), sizeof(int),
                       MPI_INFO_NULL, MPI_COMM_WORLD, &windowCount);  // Create a window for the count of local nodes

        MPI_Win_create(nodesList.data(), nodesList.size() * sizeof(unsigned int), sizeof(unsigned int),
                       MPI_INFO_NULL, MPI_COMM_WORLD, &windowNodes);

        MPI_Win_create(communitiesList.data(), communitiesList.size() * sizeof(int), sizeof(int),
                       MPI_INFO_NULL, MPI_COMM_WORLD, &windowCommunities);

        vector<int> receivedCounts(size, 0);  // Count of received nodes
        vector<unsigned int> receivedNodes;   // List of received nodes
        vector<int> receivedCommunities;      // List of received communities

        MPI_Win_lock_all(MPI_MODE_NOCHECK, windowCount);  // Lock the window for all processes
        for (unsigned int i = 0; i < size; ++i) {
            if (i == rank) continue;                                                 // Skip if the process is the current one
            MPI_Get(&receivedCounts[i], 1, MPI_INT, i, 0, 1, MPI_INT, windowCount);  // Get the count of local nodes from the remote process
        }
        MPI_Win_unlock_all(windowCount);  // Unlock the window after getting counts

        receivedNodes.resize(accumulate(receivedCounts.begin(), receivedCounts.end(), 0));  // Resize the received nodes list
        receivedCommunities.resize(receivedNodes.size());                                   // Resize the received communities list

        MPI_Win_lock_all(MPI_MODE_NOCHECK, windowNodes);        // Lock the window for all processes
        MPI_Win_lock_all(MPI_MODE_NOCHECK, windowCommunities);  // Lock the communities window
        for (unsigned int i = 0; i < size; ++i) {
            if (i == rank) continue;  // Skip if the process is the current one

            unsigned long size = receivedCounts[i];  // Size of the community list

            if (size == 0) continue;  // Skip if the size is zero

            MPI_Get(receivedNodes.data() + accumulate(receivedCounts.begin(), receivedCounts.begin() + i, 0),
                    size, MPI_UNSIGNED, i, 0, size, MPI_UNSIGNED, windowNodes);  // Get the remote nodes

            MPI_Get(receivedCommunities.data() + accumulate(receivedCounts.begin(), receivedCounts.begin() + i, 0),
                    size, MPI_INT, i, 0, size, MPI_INT, windowCommunities);  // Get the communities of the remote nodes
        }
        MPI_Win_unlock_all(windowNodes);        // Unlock the window after getting remote communities
        MPI_Win_unlock_all(windowCommunities);  // Unlock the communities window

        for (unsigned int i = 0; i < receivedNodes.size(); ++i) {
            unsigned int node = receivedNodes[i];             // Get the global id of the node
            unsigned int community = receivedCommunities[i];  // Get the global id of the community

            if (node >= c.n2c.size()) {
                c.n2c.resize(max(node + 1, community + 1), -1);  // Resize the community mapping if needed
            }

            c.n2c[node] = community;  // Update the node's community
        }

        MPI_Win_free(&windowCount);        // Free the offsets window
        MPI_Win_free(&windowNodes);        // Free the nodes window
        MPI_Win_free(&windowCommunities);  // Free the communities window
        // ----------------------------------------------------------------
        // Resolve duality conflicts
        // ----------------------------------------------------------------
        resolveDuality(c.n2c);  // Resolve duality conflicts in the community mapping

        map<int, vector<unsigned int>> c2n;
        set<unsigned int> communitiesSet;  // Set to store unique communities
        for (size_t i = 0; i < c.n2c.size(); ++i) {
            c2n[c.n2c[i]].push_back(i);
            communitiesSet.insert(c.n2c[i]);  // Insert the community into the set
        }

        if (c2n.size() < foundCommunities) {
            foundCommunities = c2n.size();  // Update the count of found communities if it has decreased
        } else {
            break;
        }

        // ----------------------------------------------------------------
        // Collecting missing nodes for coarsening
        // ----------------------------------------------------------------
        getList.resize(0);         // List of remote nodes to retrieve
        getProcessList.resize(0);  // List of processes to get remote nodes from
        for (const auto &community : c2n) {
            for (const auto &node : community.second) {
                if (localGraph.isCollected(node)) {
                    continue;  // Skip if the node is already collected
                }

                if (find(getList.begin(), getList.end(), node) == getList.end()) {
                    getList.push_back(node);  // Add remote node to the get list if not already present
                }
                if (find(getProcessList.begin(), getProcessList.end(), p.owner(node)) == getProcessList.end()) {
                    getProcessList.push_back(p.owner(node));  // Add the owner process of the remote node to the get process list
                }
            }
        }

        collectMissingNodes(localGraph, c, p, rank, getList, getProcessList);  // Collect missing nodes from remote processes
        newModularity = computeModularity(communitiesSet, c, p, rank);         // Compute the modularity of the current community structure

        // if (rank == 0) {
        //     cerr << "Modularity before step " << step << ": " << modularity << endl;
        //     cerr << "Modularity after step " << step << ": " << newModularity << endl;
        //     cerr << "Modularity gain: " << newModularity - modularity << endl;
        //     cerr << "Communities found: " << c2n.size() << endl;
        // }

        modularity = newModularity;  // Update the modularity for the next iteration

        // ----------------------------------------------------------------
        // Communities renumbering
        // ----------------------------------------------------------------
        map<int, int> old2new;  // Map to store new community indices
        vector<int> new2old;    // Vector to store community IDs
        for (const auto &community : c2n) {
            old2new[community.first] = old2new.size();  // Assign new indices to communities
            new2old.push_back(community.first);         // Store the old community ID
        }

        // ----------------------------------------------------------------
        // Graph coarsening
        // ----------------------------------------------------------------
        Graph coarseGraph;
        getList.resize(0);         // List of remote nodes to retrieve
        getProcessList.resize(0);  // List of processes to get remote nodes from
        for (const auto &community : c2n) {
            if (old2new[community.first] % size != rank) continue;    // Skip communities not owned by the current process
            coarseGraph.localNodes.insert(old2new[community.first]);  // Add the community to the local nodes of the coarse graph

            for (const auto &node : community.second) {      // Add the node to the local nodes of the coarse graph
                if (localGraph.isCollected(node)) continue;  // Skip if the node is already collected

                if (find(getList.begin(), getList.end(), node) == getList.end()) {
                    getList.push_back(node);  // Add remote node to the get list if not already present
                }

                if (find(getProcessList.begin(), getProcessList.end(), p.owner(node)) == getProcessList.end()) {
                    getProcessList.push_back(p.owner(node));  // Add the owner process of the remote node to the get process list
                }
            }
        }

        collectMissingNodes(localGraph, c, p, rank, getList, getProcessList);   // Collect missing nodes from remote processes
        p.updatePartition(new2old);                                             // Update the partition map with the new community indices
        coarseGraph.coarse(localGraph, p, rank, c.n2c, c2n, new2old, old2new);  // Coarsen the graph
        localGraph = move(coarseGraph);                                         // Move the coarse graph to the local graph

        // if (rank == 0) cerr << "Communities found: " << c2n.size() << endl;
        // for (const auto &comm : c2n) {
        //     cerr << "Community " << comm.first << ": ";
        //     for (const auto &node : comm.second) {
        //         cerr << node << " ";  // Log the nodes in the community
        //     }
        //     cerr << endl;  // New line after each community
        // }

        step++;  // Increment the step counter for each remote node processed
    }

    double endTime = MPI_Wtime();  // Get the end time of the computation

    if (rank == 0) {
        cerr << "Finished in " << step << " steps." << endl;
        cerr << "Total time: " << endTime - startTime << " seconds." << endl;
        cerr << "Partitioning time: " << partitionTime - startTime << " seconds." << endl;
        cerr << "Community detection time: " << endTime - partitionTime << " seconds." << endl;
        cerr << "Total communities found: " << foundCommunities << endl;

        string outputFile = "out.csv";
        bool fileExists = filesystem::exists(outputFile);  // Check if the output file already exists

        ofstream outFile("out.csv", ios::app);
        if (!outFile.is_open()) {
            cerr << "Error opening output file." << endl;
            return 1;
        }

        if (!fileExists)
            outFile << "size,startTime,partitionTime,endTime,foundCommunities,fileName" << endl;  // Write header if file does not exist

        outFile << size << "," << startTime << "," << partitionTime << "," << endTime << "," << foundCommunities << "," << fileName << endl;  // Write the results to the output file
        outFile.close();                                                                                                                      // Close the output file
    }

    MPI_Finalize();
    return 0;
}

void collectMissingNodes(Graph &localGraph, Community &c, Partitioner &p, int rank, vector<unsigned int> &getList, vector<unsigned int> &getProcessList) {
    bool weighted = !localGraph.weights.empty();  // Check if the graph is weighted

    vector<unsigned long> offsets;   // Offsets for each neighbour list
    vector<unsigned int> nodesList;  // Flattened list of neighbours
    vector<int> weightsList;         // Flattened list of weights (if weighted)

    MPI_Win windowOffsets;
    MPI_Win windowNodes;
    MPI_Win windowWeights;  // Window for the weights (if weighted)

    offsets.push_back(0);  // Initialize the first offset to 0
    for (size_t i = 0; i < localGraph.neighboursList.size(); ++i) {
        const auto &vec = localGraph.neighboursList[i];
        nodesList.insert(nodesList.end(), vec.begin(), vec.end());  // Flatten the neighbour lists
        offsets.push_back(nodesList.size());                        // Update the offset for the next neighbour list
    }

    if (weighted) {
        for (size_t i = 0; i < localGraph.weights.size(); ++i) {
            const auto &vec = localGraph.weights[i];
            weightsList.insert(weightsList.end(), vec.begin(), vec.end());  // Flatten the weights lists
        }
    }

    MPI_Win_create(nodesList.data(), nodesList.size() * sizeof(unsigned int), sizeof(unsigned int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowNodes);

    MPI_Win_create(offsets.data(), offsets.size() * sizeof(unsigned long), sizeof(unsigned long),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowOffsets);

    if (weighted) {
        MPI_Win_create(weightsList.data(), weightsList.size() * sizeof(int), sizeof(int),
                       MPI_INFO_NULL, MPI_COMM_WORLD, &windowWeights);
    }

    // Get offsets of the neighbour lists from all processes
    vector<unsigned long> neighbourOffsets(2 * getList.size());  // Offsets for each node's neighbour list
    int it = 0;

    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowOffsets);  // Lock the window for all processes
    for (const auto &node : getList) {
        unsigned int owner = p.owner(node);  // Get the owner process of the remote node
        if (owner == rank) continue;         // Skip if the process is the current one

        MPI_Get(&neighbourOffsets[2 * it], 2, MPI_UNSIGNED_LONG,
                owner, node, 2, MPI_UNSIGNED_LONG, windowOffsets);

        it++;
    }
    MPI_Win_unlock_all(windowOffsets);  // Unlock the window after getting offsets

    it = 0;                                                           // Reset iterator for neighbour offsets
    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowNodes);                  // Lock the window for all processes
    if (weighted) MPI_Win_lock_all(MPI_MODE_NOCHECK, windowWeights);  // Lock the weights window if weighted
    for (const auto &node : getList) {
        unsigned int owner = p.owner(node);  // Get the owner process of the remote node
        if (owner == rank) continue;         // Skip if the process is the current one

        unsigned long offset = neighbourOffsets[2 * it];          // Get the starting byte for the neighbour list
        unsigned long offsetNext = neighbourOffsets[2 * it + 1];  // Get the next starting byte for the neighbour list
        unsigned long size = offsetNext - offset;                 // Calculate the size of the neighbour list
        it++;                                                     // Increment the iterator for the next node

        if (node >= localGraph.neighboursList.size()) {
            localGraph.neighboursList.resize(node + 1);  // Resize the neighbours list if needed
        }

        localGraph.remoteNodes.insert(node);           // Add the remote node to the set of remote nodes
        localGraph.neighboursList[node].resize(size);  // Resize the neighbour list for the remote node

        if (size == 0) continue;  // Skip if the size of the neighbour list is zero
        MPI_Get(localGraph.neighboursList[node].data(), size, MPI_UNSIGNED,
                owner, offset, size, MPI_UNSIGNED, windowNodes);  // Get the neighbour list for the remote node

        if (weighted) {
            if (node >= localGraph.weights.size()) {
                localGraph.weights.resize(node + 1);  // Resize the weights list if needed
            }

            localGraph.weights[node].resize(size);  // Resize the weights list for the remote node
            MPI_Get(localGraph.weights[node].data(), size, MPI_INT,
                    owner, offset, size, MPI_INT, windowWeights);  // Get the weights for the remote node
        }
    }
    MPI_Win_unlock_all(windowNodes);                  // Unlock the window after getting remote neighbours
    if (weighted) MPI_Win_unlock_all(windowWeights);  // Unlock the weights window if weighted

    // Update the local graph with remote nodes and their communities
    c.resize();  // Resize the community structure based on the local graph
    for (const auto &node : getList) {
        int community = c.n2c[node] < 0 ? node : c.n2c[node];  // Get the community of the remote node, or use the node itself if not assigned
        unsigned int degree = localGraph.degree(node);         // Get the degree of the remote node
        c.updateRemote(node, community, degree);               // Update the community structure with the remote node
    }

    MPI_Win_free(&windowNodes);    // Free the window after use
    MPI_Win_free(&windowOffsets);  // Free the offsets window
}

void resolveDuality(vector<int> &n2c) {
    for (int node = 0; node < static_cast<int>(n2c.size()); ++node) {
        int current = node;
        unordered_set<int> visited;

        int root = current;
        while (n2c[root] != root) {
            if (!visited.insert(root).second) {
                root = *min_element(visited.begin(), visited.end());
                break;
            }
            root = n2c[root];
        }

        current = node;
        while (n2c[current] != root) {
            int next = n2c[current];
            n2c[current] = root;
            current = next;
        }
    }
}

double computeModularity(set<unsigned int> &communities, Community &c, Partitioner &p, int rank) {
    double modularity = 0.0;
    double localModularity = 0.0;                             // Compute the local modularity
    vector<double> globalModularities(p.numberOfPartitions);  // Vector to hold global modularities from all processes

    for (const auto &community : communities) {
        if (p.owner(community) != rank) continue;    // Skip if the community is not owned by the current process
        localModularity += c.modularity(community);  // Accumulate the modularity for the local community
    }

    MPI_Win windowModularity;

    MPI_Win_create(&localModularity, sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &windowModularity);  // Create a window

    MPI_Win_lock_all(MPI_MODE_NOCHECK, windowModularity);  // Lock the window for all processes
    for (int i = 0; i < p.numberOfPartitions; ++i) {
        if (i == rank) {
            globalModularities[i] = localModularity;  // Set the local modularity for the current process
        } else {
            MPI_Get(&globalModularities[i], 1, MPI_DOUBLE, i, 0, 1, MPI_DOUBLE, windowModularity);  // Get the modularity from other processes
        }
    }
    MPI_Win_unlock_all(windowModularity);  // Unlock the window after reduction
    MPI_Win_free(&windowModularity);       // Free the window

    modularity = accumulate(globalModularities.begin(), globalModularities.end(), 0.0);  // Sum up the modularities from all processes

    return modularity;  // Return the global modularity
}