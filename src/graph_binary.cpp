#include "graph_binary.h"

#include <fstream>
#include <iostream>
#include <map>
using namespace std;

GraphBinary::GraphBinary() : nNodes(0), nEdges(0), totalWeight(0.0) {}

GraphBinary::GraphBinary(char* inFile, char* weightsInFile) {
    ifstream finput(inFile, fstream::in | fstream::binary);
    if (!finput.is_open()) {
        cerr << "Error opening file: " << inFile << endl;
        return;
    }

    // Read number of nodes
    finput.read((char*)(&nNodes), 4);
    cout << "Number of Nodes: " << nNodes << endl;

    // Read cumulative degree for each node
    degrees.resize(nNodes);
    finput.read((char*)(&degrees[0]), nNodes * 8);
    cout << "Degrees read successfully." << endl;

    // Read edges and weights
    nEdges = degrees[nNodes - 1];
    edges.resize(nEdges);
    finput.read((char*)(&edges[0]), nEdges * 4);

    weights.resize(nEdges);
    finput.read((char*)(&weights[0]), nEdges * 8);

    finput.close();

    // Compute total weight
    for (unsigned int i = 0; i < nEdges; i++) {
        totalWeight += weightedDegree(i);
    }
}

unsigned int GraphBinary::nNeighbours(unsigned int node) {
    if (node >= nNodes || node < 0) {
        throw out_of_range("Node index out of range");
    }

    if (node == 0)
        return degrees[0];
    else
        return degrees[node] - degrees[node - 1];
}

double GraphBinary::weightedDegree(unsigned int node) {
    if (node >= nNodes || node < 0) {
        throw out_of_range("Node index out of range");
    }

    pair<vector<unsigned int>::iterator, vector<double>::iterator> p = neighbours(node);
    double res = 0;
    for (unsigned int i = 0; i < nNeighbours(node); i++) {
        res += (double)*(p.second + i);
    }
    return res;
}

double GraphBinary::selfLoops(unsigned int node) {
    if (node >= nNodes || node < 0) {
        throw out_of_range("Node index out of range");
    }

    pair<vector<unsigned int>::iterator, vector<double>::iterator> p = neighbours(node);
    double res = 0;
    for (unsigned int i = 0; i < nNeighbours(node); i++) {
        if (*(p.first + i) == node) {
            res += *(p.second + i);
        }
    }
    return res;
}

pair<vector<unsigned int>::iterator, vector<double>::iterator> GraphBinary::neighbours(unsigned int node) {
    if (node >= nNodes || node < 0) {
        throw out_of_range("Node index out of range");
    }

    if (node == 0)
        return make_pair(edges.begin(), weights.begin());
    else if (weights.size() != 0)
        return make_pair(edges.begin() + degrees[node - 1], weights.begin() + degrees[node - 1]);
}