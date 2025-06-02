
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
#include "graph_extended.h"

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

char *filename = NULL;
char *remoteEdgesFile = NULL;
char *processor_name;
// char *remoteEdgesFile = NULL;
char *filename_part = NULL;
int type = 0;
int nb_pass = 0;
double precision = 0.000001;
int display_level = -2;
int k1 = 16;

int num_threads = 1;

bool verbose = false;

vector<pair<int, int>> remoteMap(10, make_pair(-1, -1));

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

void usage(char *prog_name, const char *more) {
    cerr << more;
    cerr << "usage: " << prog_name << " input_file [-w weight_file] [-p part_file] [-q epsilon] [-l display_level] [-v] [-h]" << endl
         << endl;
    cerr << "input_file: file containing the graph to decompose in communities." << endl;
    cerr << "-w file\tread the graph as a weighted one (weights are set to 1 otherwise)." << endl;
    cerr << "-p file\tstart the computation with a given partition instead of the trivial partition." << endl;
    cerr << "\tfile must contain lines \"node community\"." << endl;
    cerr << "-q eps\ta given pass stops when the modularity is increased by less than epsilon." << endl;
    cerr << "-l k\tdisplays the graph of level k rather than the hierachical structure." << endl;
    cerr << "\tif k=-1 then displays the hierarchical structure rather than the graph at a given level." << endl;
    cerr << "-v\tverbose mode: gives computation time, information about the hierarchy and modularity." << endl;
    cerr << "-h\tshow this usage message." << endl;
    exit(0);
}

void parse_args(int argc, char **argv) {
    if (argc < 2)
        usage(argv[0], "Bad arguments number\n");

    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'r':
                    // type = WEIGHTED;
                    remoteEdgesFile = argv[i + 1];
                    i++;
                    break;
                case 'p':
                    filename_part = argv[i + 1];
                    i++;
                    break;
                case 'q':
                    precision = atof(argv[i + 1]);
                    i++;
                    break;
                case 't':
                    num_threads = atoi(argv[i + 1]);
                    i++;
                    break;
                case 'l':
                    display_level = atoi(argv[i + 1]);
                    i++;
                    break;
                case 'k':
                    k1 = atoi(argv[i + 1]);
                    i++;
                    break;
                case 'v':
                    verbose = true;
                    break;
                default:
                    usage(argv[0], "Unknown option\n");
            }
        } else {
            if (filename == NULL)
                filename = argv[i];
            else if (remoteEdgesFile == NULL)
                remoteEdgesFile = argv[i];
            else
                usage(argv[0], "More than one filename\n");
        }
    }
}

void display_time(const char *str) {
    time_t rawtime;
    time(&rawtime);
    cerr << str << ": " << ctime(&rawtime);
}

int main(int argc, char **argv) {
    time_t time_begin, time_end;
    time(&time_begin);

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    omp_set_num_threads(num_threads);

    parse_args(argc, argv);

    stringstream rankString;
    rankString << rank;

    string s(filename);
    s += "_" + rankString.str() + ".bin";

    string r(remoteEdgesFile);
    r += "_" + rankString.str() + ".remote";

    vector<pair<int, int>> remoteMap;

    ifstream remoteFileStream(r.c_str());

    // Map to store remote edges
    vector<int> rSource;       // Local vertex IDs
    vector<int> rDestination;  // Remote vertex IDs
    vector<int> rPartition;    // Partition IDs
    MPI_Status status;

    // Read remote edges from the file
    // Each line contains a local vertex and its mapping to a remote vertex and partition
    for (string line; getline(remoteFileStream, line);) {
        vector<string> parts = split(line, ' ');
        string localV = parts[0];
        string mapping = parts[1];
        unsigned int source = atoi(localV.c_str());

        vector<string> mappingParts = split(mapping, ',');

        int destination = atoi(mappingParts[0].c_str());
        int partition = atoi(mappingParts[1].c_str());

        if (remoteMap.size() <= source) {
            remoteMap.resize(source + 1, make_pair(-1, -1));
        }

        remoteMap[source] = make_pair(destination, partition);

        rSource.push_back(source);
        rDestination.push_back(destination);
        rPartition.push_back(partition);
    }

    char *tmp = new char[s.length() + 1];
    strcpy(tmp, s.c_str());

    Community c(tmp, NULL, -1, precision);

    GraphBinary g;
    bool improvement = true;
    double modularity = c.modularity(), newModularity;
    int level = 0;

    improvement = c.step();
    newModularity = c.modularity();

    // if (++level == display_level)
    //     g.display();
    // if (display_level == -1)
    //     c.display_partition();

    g = c.graph();

#pragma omp parallel for
    for (unsigned int i = 0; i < rSource.size(); i++) {
        rSource[i] = c.n2cNew[rSource[i]];  // Re-number the source vertices
    }

    int *levelOneNodes = new int[size];
    unsigned long *levelOneDegree = new unsigned long[size];

    levelOneNodes[rank] = g.nNodes;
    levelOneDegree[rank] = g.degrees[g.degrees.size() - 1];

    MPI_Request request[2 * size - 2];  // Requests for non-blocking communication
    MPI_Status st[2 * size - 2];        // Status for the requests

    // Gather information about the number of nodes and final degree from all processes
    int id = 0;
    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Irecv(&levelOneNodes[i], 1, MPI_INT, i, 1, MPI_COMM_WORLD, &request[id++]);
            MPI_Irecv(&levelOneDegree[i], 1, MPI_UNSIGNED_LONG, i, 2, MPI_COMM_WORLD, &request[id++]);
        }
    }

    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Send(&g.nNodes, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&g.degrees[g.degrees.size() - 1], 1, MPI_UNSIGNED_LONG, i, 2, MPI_COMM_WORLD);
        }
    }

    MPI_Waitall(2 * size - 2, request, st);

    // Create gaps for re-numbering nodes across processes
    // For each process, calculate the gap based on the number of nodes
    int gap = 0;
    int *gaps = new int[size];
    for (int i = 0; i < size; i++) {
        gaps[i] = gap;
        if (i != rank)
            gap += levelOneNodes[i];
        else
            gap += g.nNodes;
    }

#pragma omp parallel for
    for (unsigned int i = 0; i < g.edges.size(); i++) {
        g.edges[i] += gaps[rank];
    }

#pragma omp parallel for
    for (unsigned int i = 0; i < rSource.size(); i++) {
        rSource[i] += gaps[rank];
    }

#pragma omp parallel for
    for (unsigned int i = 0; i < c.n2cNew.size(); i++) {
        c.n2cNew[i] += gaps[rank];
    }

    // Update the degrees of the graph based on the gaps calculated
    unsigned long degreeGap = 0;
    for (int i = 0; i < rank; i++) {
        degreeGap += levelOneDegree[i];
    }

    if (rank != 0)
        for (unsigned int i = 0; i < g.degrees.size(); i++) {
            g.degrees[i] += degreeGap;
        }

    modularity = newModularity;
    // Each process sends its graph data to the root process (rank 0)
    if (rank != 0) {
        MPI_Send(&g.nEdges, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD);          // Send number of edges
        MPI_Send(&g.nNodes, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);           // Send number of nodes
        MPI_Ssend(&g.totalWeight, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);  // Send total weight of graph

        int nEdges = g.edges.size();
        MPI_Ssend(&nEdges, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);                             // Send number of edges
        MPI_Ssend(&g.edges.front(), g.edges.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);       // Send edges
        MPI_Ssend(&g.degrees.front(), g.degrees.size(), MPI_LONG, 0, 1, MPI_COMM_WORLD);  // Send degrees

        int nWeights = g.weightList.size();
        MPI_Ssend(&nWeights, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);                            // Send number of weights
        MPI_Ssend(&g.weightList.front(), g.weightList.size(), MPI_FLOAT, 0, 1, MPI_COMM_WORLD);  // Send weights

        int rSize = rSource.size();
        MPI_Ssend(&rSize, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);                                   // Send number of remote edges
        MPI_Ssend(&rSource.front(), rSource.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);            // Send remote source vertices
        MPI_Ssend(&rDestination.front(), rDestination.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);  // Send remote destination vertices
        MPI_Ssend(&rPartition.front(), rPartition.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);      // Send remote partition IDs

        int n2cSize = c.n2cNew.size();
        MPI_Ssend(&n2cSize, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);                         // Send size of n2cNew
        MPI_Ssend(&c.n2cNew.front(), c.n2cNew.size(), MPI_INT, 0, 1, MPI_COMM_WORLD);  // Send n2cNew
    }
    if (rank == 0) {
        // Rank 0 collects data from all processes
        GraphData data[size - 1];

        int total_nodes = 0;
        for (int i = 1; i < size; i++) {
            MPI_Recv(&data[i - 1].nEdges, 1, MPI_LONG, i, 1, MPI_COMM_WORLD, &status);         // Receive number of edges
            MPI_Recv(&data[i - 1].nNodes, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);          // Receive number of nodes
            MPI_Recv(&data[i - 1].totalWeight, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);  // Receive total weight of graph

            int nEdges;
            MPI_Recv(&nEdges, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);  // Receive number of edges

            data[i - 1].edges.resize(nEdges);
            data[i - 1].degrees.resize(data[i - 1].nNodes);
            MPI_Recv(&data[i - 1].edges.front(), nEdges, MPI_INT, i, 1, MPI_COMM_WORLD, &status);                 // Receive edges
            MPI_Recv(&data[i - 1].degrees.front(), data[i - 1].nNodes, MPI_LONG, i, 1, MPI_COMM_WORLD, &status);  // Receive degrees

            int nWeights;
            data[i - 1].weights.resize(nWeights);
            MPI_Recv(&nWeights, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);                              // Receive number of weights
            MPI_Recv(&data[i - 1].weights.front(), nWeights, MPI_FLOAT, i, 1, MPI_COMM_WORLD, &status);  // Receive weights

            int rSize;
            data[i - 1].rSource.resize(rSize);
            data[i - 1].rDestination.resize(rSize);
            data[i - 1].rPartition.resize(rSize);
            MPI_Recv(&rSize, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);                                 // Receive number of remote edges
            MPI_Recv(&data[i - 1].rSource.front(), rSize, MPI_INT, i, 1, MPI_COMM_WORLD, &status);       // Receive remote source vertices
            MPI_Recv(&data[i - 1].rDestination.front(), rSize, MPI_INT, i, 1, MPI_COMM_WORLD, &status);  // Receive remote destination vertices
            MPI_Recv(&data[i - 1].rPartition.front(), rSize, MPI_INT, i, 1, MPI_COMM_WORLD, &status);    // Receive remote partition IDs

            int n2cSize;
            data[i - 1].n2c.resize(n2cSize);
            MPI_Recv(&n2cSize, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);                        // Receive size of n2cNew
            MPI_Recv(&data[i - 1].n2c.front(), n2cSize, MPI_INT, i, 1, MPI_COMM_WORLD, &status);  // Receive n2cNew

            total_nodes += data[i - 1].nNodes;
        }

        // Construct new graph for next iteration
        GraphExtended newG = GraphExtended();

        newG.nNodes = g.nNodes;
        newG.nEdges = 0;
        newG.totalWeight = 0;

        // Merge received graphs into newG
        for (int i = 0; i < size; i++) {
            if (i == 0) {  // Local graph
                newG.degrees->extend(g.degrees);
                newG.edges->extend(g.edges);
                newG.weights->extend(g.weightList);
                newG.nEdges += g.nEdges;
                newG.totalWeight += g.totalWeight;
            } else {
                newG.degrees->extend(data[i - 1].degrees);
                newG.edges->extend(data[i - 1].edges);
                newG.weights->extend(data[i - 1].weights);
                newG.nEdges += data[i - 1].nEdges;
                newG.totalWeight += data[i - 1].totalWeight;
                newG.nNodes += data[i - 1].nNodes;
            }
        }

        map<int, vector<unsigned int>> re;
        map<int, vector<double>> rw;
        for (int i = 0; i < size; i++) {
            // For root process, merge remote edges
            if (i == 0) {
                // Create a map to store remote edges and their weights
                // Key is a pair of (node, community), value is the weight
                map<pair<int, unsigned int>, double> m;
                map<pair<int, unsigned int>, double>::iterator it;
                for (unsigned int j = 0; j < rSource.size(); j++) {
                    int destination = rDestination[j];
                    int destinationPartition = rPartition[j];
                    int target = data[destinationPartition - 1].n2c[destination];

                    it = m.find(make_pair(rSource[j], target));

                    if (it == m.end()) {
                        m.insert(make_pair(make_pair(rSource[j], target), 1.0f));
                    } else {
                        it->second += 1.0f;
                    }
                }

                newG.nEdges += m.size();

                map<int, vector<unsigned int>>::iterator remoteEdgeIt;
                map<int, vector<double>>::iterator remoteWeightsIt;

                for (it = m.begin(); it != m.end(); it++) {
                    pair<int, unsigned int> e = it->first;
                    double w = it->second;

                    remoteEdgeIt = re.find(e.first);
                    remoteWeightsIt = rw.find(e.first);

                    if (remoteEdgeIt == re.end()) {
                        vector<unsigned int> list;
                        list.push_back(e.second);
                        re.insert(make_pair(e.first, list));

                        vector<double> wList;
                        wList.push_back(w);
                        rw.insert(make_pair(e.first, wList));
                    } else {
                        remoteEdgeIt->second.push_back(e.second);
                        if (remoteWeightsIt != rw.end()) {
                            remoteWeightsIt->second.push_back(w);
                        }
                    }
                }
            } else {
                map<pair<int, unsigned int>, float> m;
                map<pair<int, unsigned int>, float>::iterator it;

                for (unsigned int j = 0; j < data[i - 1].rSource.size(); j++) {
                    int destination = data[i - 1].rDestination[j];
                    int destinationPartition = data[i - 1].rPartition[j];
                    int target;

                    if (destinationPartition == 0) {
                        target = c.n2cNew[destination];
                    } else {
                        target = data[destinationPartition - 1].n2c[destination];
                    }

                    it = m.find(make_pair(data[i - 1].rSource[j], target));
                    if (it == m.end()) {
                        m.insert(make_pair(make_pair(data[i - 1].rSource[j], target), 1.0f));
                    } else {
                        it->second += 1.0f;
                    }
                }

                newG.nEdges += m.size();

                map<int, vector<unsigned int>>::iterator remoteEdgeIt;
                map<int, vector<double>>::iterator remoteWIt;

                for (it = m.begin(); it != m.end(); it++) {
                    pair<int, int> e = it->first;
                    double w = it->second;

                    remoteEdgeIt = re.find(e.first);
                    remoteWIt = rw.find(e.first);

                    if (remoteEdgeIt == re.end()) {
                        vector<unsigned int> list;
                        list.push_back(e.second);
                        re.insert(make_pair(e.first, list));

                        vector<double> wList;
                        wList.push_back(w);
                        rw.insert(make_pair(e.first, wList));

                    } else {
                        remoteEdgeIt->second.push_back(e.second);
                        if (remoteWIt != rw.end()) {
                            remoteWIt->second.push_back(w);
                        }
                    }
                }
            }
        }

        // Update the graph with remote edges
        newG.addRemoteEdges(re, rw);

        // Print results
        newG.print();
    }

    MPI_Finalize();

    time(&time_end);

    display_time("End");
    cerr << "Total duration: " << (time_end - time_begin) << " sec." << endl;

    return 0;
}
