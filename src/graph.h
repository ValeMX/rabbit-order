#ifndef GRAPH_H
#define GRAPH_H

#include <string>
#include <vector>
using namespace std;

class Graph {
   public:
    vector<vector<pair<int, double>>> edges;

    Graph(string fileName = "");
    Graph(vector<pair<int, int>> edgeList);
    Graph(vector<pair<int, int>> edgeList, vector<double> weights);

    void addEdge(int src, int dst, double weight);

    void clean();  // Clean the graph by removing duplicate edges
    void renumber();
    void saveBinary(char* outFile, char* weightsOutFile);  // Save the graph in binary format
};

#endif  // GRAPH_H