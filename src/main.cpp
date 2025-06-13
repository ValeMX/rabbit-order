
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

void collectMissingNodes(Graph &localGraph, Community &c, Partitioner &p, int rank, vector<unsigned int> &getList, vector<unsigned int> &getProcessList);
void resolveDuality(vector<int> &n2c);

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
    Graph coarseGraph;

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

            count = p.partitionMap.size();
            MPI_Send(&count, 1, MPI_INT, i, 4, MPI_COMM_WORLD);                          // Send size of partition map
            MPI_Send(p.partitionMap.data(), count, MPI_UNSIGNED, i, 5, MPI_COMM_WORLD);  // Send partition map
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

        MPI_Recv(&count, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);  // Receive size of partition map
        p.partitionMap.resize(count);
        MPI_Recv(p.partitionMap.data(), p.partitionMap.size(), MPI_UNSIGNED, 0, 5, MPI_COMM_WORLD, &status);  // Receive partition map
    }

    // if (rank == 0) {
    //     cerr << "localGraph.edgeList:" << endl;
    //     for (const auto &edge : localGraph.edgeList) {
    //         cerr << "(" << edge.first << ", " << edge.second << ")" << endl;
    //     }
    // }

    if (rank == 0) cerr << "Initializing local graph..." << endl;
    localGraph.init();

    for (int step = 0; step < 1000; step++) {
        // ----------------------------------------------------------------
        // Graph Initialization
        // ----------------------------------------------------------------
        MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes before proceeding
        if (rank == 0) cerr << "\nSTEP: " << step << endl;

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

        if (rank == 0) cerr << "Collecting remote edges..." << endl;

        // if (rank == 0 && step == 1) {
        //     cerr << "Local nodes: " << localGraph.localNodes.size() << endl;
        //     for (const auto &node : localGraph.localNodes) {
        //         cerr << "Local node: " << node << " from process " << p.owner(node) << endl;  // Log local nodes
        //     }
        // }

        // if (rank == 1 && step == 1) {
        //     cerr << "Collecting remote edges for step 1..." << endl;
        //     for (const auto &node : getList) {
        //         cerr << node << " ";  // Log remote nodes to be collected
        //     }
        //     cerr << endl;
        // }

        collectMissingNodes(localGraph, c, p, rank, getList, getProcessList);  // Collect missing nodes from remote processes

        // ----------------------------------------------------------------
        // Community Detection
        // ----------------------------------------------------------------
        if (rank == 0) cerr << "Starting community detection..." << endl;

        if (rank == 0) cerr << "Modularity before steps: " << c.modularity() << endl;

        c.step();

        // if (rank == 0) c.print();

        if (rank == 0) cerr << "Modularity: " << c.modularity() << endl;

        // ----------------------------------------------------------------
        // Exchange updated communities
        // ----------------------------------------------------------------
        if (rank == 0) cerr << "Exchanging updated communities..." << endl;

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

        if (rank == 0) cerr << "Windows for communities created. Proceeding to exchange communities..." << endl;

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

        if (rank == 0) cerr << "Communities exchanged. Updating local communities..." << endl;

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

        MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes before proceeding

        // ----------------------------------------------------------------
        // Resolve duality conflicts
        // ----------------------------------------------------------------
        if (rank == 0) cerr << "Resolving duality conflicts..." << endl;

        // Check if all communities are updated correctly
        for (const auto &comm : c.n2c) {
            if (comm < 0 || comm >= c.n2c.size()) {
                cerr << "A. Rank " << rank << " Error: Node " << &comm - &c.n2c[0] << " is not assigned to any community: " << comm << " in size " << c.n2c.size() << endl;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes before proceeding
        if (rank == 0) cerr << "All processes synchronized. Proceeding to resolve duality conflicts..." << endl;

        // Resolve duality
        resolveDuality(c.n2c);  // Resolve duality conflicts in the community mapping

        map<int, vector<unsigned int>> c2n;
        for (size_t i = 0; i < c.n2c.size(); ++i) {
            c2n[c.n2c[i]].push_back(i);
        }

        if (c2n.size() <= size) {
            if (rank == 0) {
                cerr << "Warning: Number of communities (" << c2n.size() << ") is less than or equal to the number of processes (" << size << "). This may lead to suboptimal partitioning." << endl;
            }
            break;  // Exit the loop if the number of communities is less than or equal to the number of processes
        }

        for (const auto &comm : c.n2c) {
            if (comm < 0 || comm >= c.n2c.size()) {
                cerr << "B. Rank " << rank << " Error: Node " << &comm - &c.n2c[0] << " is not assigned to any community: " << comm << " in size " << c.n2c.size() << endl;
            }
        }

        if (rank == 0) cerr << "Duality conflicts resolved." << endl;

        // if (rank == 0) {
        //     cerr << "Local communities: " << endl;
        //     for (const auto &kv : c2n) {
        //         cout << "Community " << kv.first << ": ";
        //         for (const auto &node : kv.second) {
        //             cout << node << " ";
        //         }
        //         cout << endl;
        //     }
        // }

        // ----------------------------------------------------------------
        // Collecting missing nodes for coarsening
        // ----------------------------------------------------------------
        if (rank == 0) cerr << "Collecting missing nodes for coarsening..." << endl;

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

        for (const auto &comm : c.n2c) {
            if (comm < 0 || comm >= c.n2c.size()) {
                cerr << "D. Rank " << rank << " Error: Node " << &comm - &c.n2c[0] << " is not assigned to any community: " << comm << " in size " << c.n2c.size() << endl;
            }
        }

        // ----------------------------------------------------------------
        // Graph coarsening
        // ----------------------------------------------------------------

        MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes before proceeding
        if (rank == 0) cerr << "Coarsening graph..." << endl;

        // Create a new graph for coarsening
        Graph coarseGraph;
        coarseGraph.coarse(localGraph, p, rank, c.n2c, c2n);  // Coarsen the graph based on the local graph and communities

        localGraph = move(coarseGraph);  // Move the coarse graph to the local graph
    }

    MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes before proceeding
    if (rank == 0) cerr << "Local communities updated. Community detection completed." << endl;
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

    MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes before proceeding
    if (rank == 0) cerr << "Windows setup complete. Initiating remote edge collection..." << endl;

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

    MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes before proceeding
    if (rank == 0) cerr << "Offsets collected. Proceeding to get remote neighbours..." << endl;

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

    MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes before proceeding
    if (rank == 0) cerr << "Remote neighbours collected. Proceeding to update communities..." << endl;

    // Update the local graph with remote nodes and their communities
    c.resize();  // Resize the community structure based on the local graph
    for (const auto &node : getList) {
        int community = c.n2c[node] < 0 ? node : c.n2c[node];  // Get the community of the remote node, or use the node itself if not assigned
        unsigned int degree = localGraph.degree(node);         // Get the degree of the remote node
        localGraph.totalWeight += degree;                      // Update the total weight of the graph
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
