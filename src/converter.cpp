// Credits: https://github.com/usc-cloud/parallel-louvain-modularity/blob/master/src/parallel/converter.cpp

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "graph.h"

using namespace std;

char *infile = NULL;
char *outfile = NULL;
char *partition_file = NULL;
int numberOfPartitions = 1;

void usage(char *prog_name, const char *more) {
    cerr << more;
    cerr << "usage: " << prog_name << " -i input_file -o outfile -p partition file -n number of partitions" << endl
         << endl;
    exit(0);
}

void parse_args(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'i':
                    if (i == argc - 1)
                        usage(argv[0], "Infile missing\n");
                    infile = argv[i + 1];
                    i++;
                    break;
                case 'o':
                    if (i == argc - 1)
                        usage(argv[0], "Outfile missing\n");
                    outfile = argv[i + 1];
                    i++;
                    break;
                case 'p':
                    if (i == argc - 1)
                        usage(argv[0], "partition missing\n");
                    partition_file = argv[i + 1];
                    i++;
                    break;
                case 'n':
                    if (i == argc - 1)
                        usage(argv[0], "number of partitions missing\n");
                    numberOfPartitions = atoi(argv[i + 1]);
                    i++;
                    break;
                default:
                    usage(argv[0], "Unknown option\n");
            }
        } else {
            usage(argv[0], "More than one fileName\n");
        }
    }
    if (infile == NULL || outfile == NULL)
        usage(argv[0], "In or outfile missing\n");
}

/*
 * This will accept original graph file, partition file and number of partitions.
 * This will partition the graph in to multiple bin files based on number of partitions.
 * Each file will will have its local numbering
 * There is a remote file which contain the remote edges.
 * This will have the remote edges and associated partition number
 */
int main(int argc, char **argv) {
    parse_args(argc, argv);

    //  infile = "~/metis_test/test_louvan.graph";
    //  outfile = "~/metis_test/test_part";
    //  partition_file = "~/metis_test/test.graph.part.2";
    //  numberOfPartitions = 2;

    cout << "Intput file: " << infile << " Output " << outfile << " Partitions : " << partition_file << " number of Part : " << numberOfPartitions << endl;
    ifstream finput;
    finput.open(infile, fstream::in);

    vector<pair<int, int>> edgeList;
    int nb_links = 0;

    // Read the input file and create the edge list
    unsigned int psrc = 0, pdest = 0;
    while (!finput.eof()) {
        unsigned int src, dest;
        finput >> src >> dest;

        if (finput.eof())
            break;

        // Check if duplicate or blank line
        if ((src != psrc) || (dest != pdest)) {
            // cout << src << " " << dest << endl;
            edgeList.push_back(make_pair(src, dest));
            nb_links++;
            psrc = src;
            pdest = dest;
        } else {
            continue;
        }

        if (finput.eof())
            break;
    }
    finput.close();

    cout << "Loading graph done" << edgeList.size() << endl;

    // Load the partition file: each line contains the partition number for each vertex
    vector<int> partitionMap;
    ifstream fpartition;
    fpartition.open(partition_file, fstream::in);

    // vid will be used to track the number of vertices
    int vid = 0;
    while (!fpartition.eof()) {
        int partition;
        fpartition >> partition;
        if (fpartition.eof())
            break;
        partitionMap.push_back(partition);
        vid++;
    }

    cout << "Loading Partition done size : " << vid << endl;

    // Load the edges into partitions based on the partition map:
    // if both source and destination vertices belong to the same partition,
    // the edge is added to the local partition; otherwise, it is added to the remote edges.
    vector<vector<pair<int, int>>> partitions;
    partitions.resize(numberOfPartitions);
    vector<vector<pair<int, int>>> remoteEdges;
    remoteEdges.resize(numberOfPartitions);

    for (unsigned int i = 0; i < edgeList.size(); i++) {
        unsigned int src = edgeList[i].first;
        unsigned int dst = edgeList[i].second;

        if ((src > partitionMap.size()) || (dst > partitionMap.size())) {
            cout << " ERROR : " << dst << " : " << src << " " << endl;
        }

        if (partitionMap[src - 1] == partitionMap[dst - 1]) {
            // in same partition
            partitions[partitionMap[src - 1]].push_back(edgeList[i]);
        } else {
            // in different partitions
            remoteEdges[partitionMap[src - 1]].push_back(edgeList[i]);
        }
    }

    cout << "Partitioning done" << endl;
    vector<pair<int, int>> oldToNewMap;
    oldToNewMap.resize(vid + 1);

    // Renumber the vertices in each partition and create a mapping from old to new vertex IDs
    // Each partition will have its own numbering starting from 1
    // The mapping will be stored in oldToNewMap, where the index is the old vertex ID
    // and the value is a pair containing the new vertex ID and the partition number
    for (int i = 0; i < numberOfPartitions; i++) {
        vector<int> map;
        map.resize(vid + 1, 0);
        unsigned int v = 1;
        // Renumber the vertices in each partition
        for (unsigned int j = 0; j < partitions[i].size(); j++) {
            int src = partitions[i][j].first;
            int dst = partitions[i][j].second;

            // Skip self-loops
            if (src == dst) {
                dst = -1;
            }

            // Check if the source vertex is already mapped
            // If not, assign a new vertex ID
            if (map[src] != 0) {
                src = map[src];
            } else {
                map[src] = v;
                oldToNewMap[src - 1] = make_pair(v, i);
                src = v++;
            }

            if (dst == -1) {
                partitions[i][j].first = src;
                partitions[i][j].second = src;
                continue;
            }

            // Check if the destination vertex is already mapped
            // If not, assign a new vertex ID
            if (map[dst] != 0) {
                dst = map[dst];
            } else {
                map[dst] = v;
                oldToNewMap[dst - 1] = make_pair(v, i);
                dst = v++;
            }
            partitions[i][j].first = src;
            partitions[i][j].second = dst;
        }
    }

    cout << "Mapping done" << endl;

    // TODO: Capire la questione dei weights
    for (int i = 0; i < numberOfPartitions; i++) {
        // Write the remote edges to separate files for each partition
        // Each file will contain the edges that are remote to that partition
        ofstream fremoteList;
        string s(outfile);
        stringstream ss;
        ss << i;
        string remoteName = s + "_" + ss.str() + ".remote";
        char *cstr = new char[remoteName.length() + 1];
        strcpy(cstr, remoteName.c_str());
        fremoteList.open(cstr, fstream::out);
        delete[] cstr;
        for (unsigned int j = 0; j < remoteEdges[i].size(); j++) {
            int src = remoteEdges[i][j].first;
            int dst = remoteEdges[i][j].second;
            fremoteList << (oldToNewMap[src - 1].first - 1) << " " << (oldToNewMap[dst - 1].first - 1) << "," << (oldToNewMap[dst - 1].second) << endl;
        }

        fremoteList.close();

        cout << "Partition " << i << " remote done" << endl;

        Graph g(partitions[i]);
        g.clean();
        string st(outfile);
        stringstream sst;
        sst << i;
        string sout = st + "_" + sst.str() + ".bin";
        string soutweight = st + "_" + sst.str() + ".weight";
        cstr = new char[sout.length() + 1];
        char *cstrweight = new char[soutweight.length() + 1];
        strcpy(cstrweight, soutweight.c_str());
        strcpy(cstr, sout.c_str());

        g.saveBinary(cstr, cstrweight);
        delete[] cstrweight;
        delete[] cstr;

        cout << "Partition " << i << " local done" << endl;
    }

    return 0;
}