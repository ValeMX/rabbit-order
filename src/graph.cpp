#include "graph.h"

#include <fstream>
#include <iostream>
#include <map>
using namespace std;

Graph::Graph(string fileName = "") {
    // Load graph data from a file
    ifstream file(fileName);

    if (!file.is_open()) {
        cerr << "Error opening file: " << fileName << endl;
        return;
    }

    int src, dst;
    double weight = 1.;
    while (file >> src >> dst >> weight) {
        addEdge(src, dst, weight);
    }

    file.close();
}

Graph::Graph(vector<pair<int, int>> edgeList) {
    double w = 1.;
    for (unsigned int i = 0; i < edgeList.size(); i++) {
        addEdge(edgeList[i].first - 1, edgeList[i].second - 1, w);
    }
}

void Graph::addEdge(int src, int dst, double weight) {
    // Ensure the list is large enough to hold the nodes
    if (edges.size() <= max(src, dst) + 1) {
        edges.resize(max(src, dst) + 1);
    }

    edges[src].emplace_back(dst, weight);
    // TODO: aggiungere anche il contrario?
}

// Remove duplicate edges and self-loops
void Graph::clean() {
    for (unsigned int i = 0; i < edges.size(); i++) {
        map<int, double> m;
        map<int, double>::iterator it;

        // For each edge in the list, add it to the map
        // If the edge already exists, update its weight
        for (unsigned int j = 0; j < edges[i].size(); j++) {
            it = m.find(edges[i][j].first);
            if (it == m.end())
                m.insert(make_pair(edges[i][j].first, edges[i][j].second));
            else
                it->second += edges[i][j].second;
        }

        // Clear the list and copy unique edges back from the map
        // Convert the map back to a vector of pairs
        vector<pair<int, double>> v;
        for (it = m.begin(); it != m.end(); it++)
            v.push_back(*it);
        edges[i].clear();
        edges[i] = v;
    }
}

void Graph::renumber() {
    vector<int> linked(edges.size(), -1);
    vector<int> renum(edges.size(), -1);
    int n = 0;

    // For each node, check if it is linked to another node
    for (unsigned int i = 0; i < edges.size(); i++) {
        for (unsigned int j = 0; j < edges[i].size(); j++) {
            linked[i] = 1;
            linked[edges[i][j].first] = 1;
        }
    }

    // Renumber each node that is linked
    for (unsigned int i = 0; i < edges.size(); i++) {
        if (linked[i] == 1) {
            renum[i] = n++;
        }
    }

    // Create a new list with renumbered nodes
    for (unsigned int i = 0; i < edges.size(); i++) {
        if (linked[i] == 1) {
            for (unsigned int j = 0; j < edges[i].size(); j++) {
                edges[i][j].first = renum[edges[i][j].first];
            }
            edges[renum[i]] = edges[i];
        }
    }
    edges.resize(n);
}

void Graph::saveBinary(char *outFile, char *weightsOutFile) {
    ofstream foutput;
    foutput.open(outFile, fstream::out | fstream::binary);

    unsigned int s = edges.size();
    cout << "Number of Nodes " << s << endl;

    // Number of nodes in 4 bytes
    foutput.write((char *)(&s), 4);

    // Cumulative degree in 8 bytes per node
    long tot = 0;
    for (unsigned int i = 0; i < s; i++) {
        tot += (long)edges[i].size();
        foutput.write((char *)(&tot), 8);
    }

    // Edges in 4 bytes per edge
    // The starting node of each edge is determined by the cumulative degree
    // The destination node is written in 4 bytes
    for (unsigned int i = 0; i < s; i++) {
        for (unsigned int j = 0; j < edges[i].size(); j++) {
            int dest = edges[i][j].first;
            foutput.write((char *)(&dest), 4);
        }
    }
    foutput.close();

    // Weights in a different file in 4 bytes per edge
    ofstream foutput_w;
    foutput_w.open(weightsOutFile, fstream::out | fstream::binary);
    for (unsigned int i = 0; i < s; i++) {
        for (unsigned int j = 0; j < edges[i].size(); j++) {
            double weight = edges[i][j].second;
            foutput_w.write((char *)(&weight), 8);
        }
    }
    foutput_w.close();
}
