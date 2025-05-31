#include "graph_extended.h"

#include <fstream>
#include <iostream>
using namespace std;

GraphExtended::GraphExtended() : nNodes(0), nEdges(0), totalWeight(0) {
    degrees = new CustomVector<unsigned long>();
    edges = new CustomVector<unsigned int>();
    weights = new CustomVector<double>();
}

GraphExtended::~GraphExtended() {
    delete degrees;
    delete edges;
    delete weights;
}

unsigned int GraphExtended::nNeighbours(unsigned int node) {
    if (node >= nNodes || node < 0) {
        throw out_of_range("Node index out of range");
    }

    if (node == 0)
        return degrees->get(0);
    else
        return degrees->get(node) - degrees->get(node - 1);
}

unsigned int GraphExtended::nRemoteNeighbours(unsigned int node) {
    if (node >= nNodes || node < 0) {
        throw out_of_range("Node index out of range");
    }

    return remoteEdges[node].size();
}

double GraphExtended::weightedDegree(unsigned int node) {
    if (node >= nNodes || node < 0) {
        throw out_of_range("Node index out of range");
    }

    if (weights->size == 0) {
        return (double)nNeighbours(node);
    } else {
        pair<vector<unsigned int>::iterator, vector<double>::iterator> p = neighbours(node);
        double res = 0;
        for (unsigned int i = 0; i < nNeighbours(node); i++) {
            res += (double)*(p.second + i);
        }
        return res;
    }
}

double GraphExtended::weightedDegreeRemote(unsigned int node) {
    if (node >= nNodes || node < 0) {
        throw out_of_range("Node index out of range");
    }

    if (remoteWeights.size() == 0) {
        return (double)nRemoteNeighbours(node);
    } else {
        pair<vector<unsigned int>::iterator, vector<double>::iterator> p = remoteNeighbours(node);
        int nRemote = nRemoteNeighbours(node);
        double res = 0;
        for (unsigned int i = 0; i < nRemote; i++) {
            res += (double)*(p.second + i);
        }
        return res;
    }
}

double GraphExtended::selfLoops(unsigned int node) {
    if (node >= nNodes || node < 0) {
        throw out_of_range("Node index out of range");
    }

    pair<vector<unsigned int>::iterator, vector<double>::iterator> p = neighbours(node);
    double res = 0;
    for (unsigned int i = 0; i < nNeighbours(node); i++) {
        if (*(p.first + i) == node) {
            if (weights->size == 0)
                res += 1.0;
            else
                res += *(p.second + i);
        }
    }
    return res;
}

pair<vector<unsigned int>::iterator, vector<double>::iterator> GraphExtended::neighbours(unsigned int node) {
    if (node >= nNodes || node < 0) {
        throw out_of_range("Node index out of range");
    }

    if (node == 0)
        return make_pair(edges->getPointer(0), weights->getPointer(0));
    else if (weights->size != 0) {
        if (degrees->get(node - 1) == degrees->get(node))  // No neighbours
            return make_pair(edges->getPointer(edges->size - 1), weights->getPointer(weights->size - 1));
        return make_pair(edges->getPointer(degrees->get(node - 1)), weights->getPointer(degrees->get(node - 1)));
    } else {
        if (degrees->get(node - 1) == degrees->get(node))  // No neighbours
            return make_pair(edges->getPointer(edges->size - 1), weights->getPointer(0));
        return make_pair(edges->getPointer(degrees->get(node - 1)), weights->getPointer(0));
    }
}

void GraphExtended::addRemoteEdges(map<int, vector<unsigned int>> re, map<int, vector<double>> rw) {
    remoteEdges = re;
    remoteWeights = rw;

    double res = 0;
    map<int, vector<double>>::iterator itw = rw.begin();
    for (; itw != rw.end(); itw++) {
        int node = itw->first;
        vector<double> ws = itw->second;

        if (ws.size() == 0) {
            int rDeg = nRemoteNeighbours(node);
            res += rDeg;
        } else
            for (unsigned int i = 0; i < ws.size(); i++) {
                res += ws[i];
            }
    }
    totalWeight += res;

    map<int, vector<unsigned int>>::iterator ite = re.begin();
    int nRemoteEdges = 0;
    for (; ite != re.end(); ite++) {
        vector<unsigned int> res = ite->second;
        nRemoteEdges += res.size();
    }
    nEdges += nRemoteEdges;
}

void GraphExtended::cleanup() {
    delete degrees;
    delete edges;
    delete weights;
}

void GraphExtended::print() {
    for (unsigned int node = 0; node < nNodes; node++) {
        pair<vector<unsigned int>::iterator, vector<double>::iterator> p = neighbours(node);
        cout << node << ":";
        for (unsigned int i = 0; i < nNeighbours(node); i++) {
            if (true) {
                if (weights->size != 0)
                    cout << " (" << *(p.first + i) << " " << *(p.second + i) << ")";
                else
                    cout << " " << *(p.first + i);
            }
        }

        p = remoteNeighbours(node);
        for (unsigned int i = 0; i < nRemoteNeighbours(node); i++) {
            if (true) {
                if (remoteWeights.find(node)->second.size() != 0)
                    cout << " (" << *(p.first + i) << " " << *(p.second + i) << ")";
                else
                    cout << " " << *(p.first + i);
            }
        }

        cout << endl;
    }
}
