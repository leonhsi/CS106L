/*
 * TODO: complete this file comment.
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include "SimpleGraph.h"

using namespace std;

const double kPi = 3.14159265358979323;   
const double kAttract = 0.001;
const double kRepel = 0.001;

void Welcome();
void GetFileName(ifstream &ifs);
void GetThreshold(chrono::milliseconds &threshold);
void CreateGraphByFile(ifstream &ifs, SimpleGraph &graph);
void IterateSimpleGraph(chrono::milliseconds threshold, SimpleGraph &graph);
bool CheckRunAgain();

// Main method
int main() {
    Welcome();

    SimpleGraph graph;
    ifstream ifs;
    chrono::milliseconds threshold;

    while (true) {
        GetFileName(ifs);
        CreateGraphByFile(ifs, graph);

        GetThreshold(threshold);
        IterateSimpleGraph(threshold, graph);

        if (CheckRunAgain() == false) break;
    }

    cout << "Simple Graph Visualization Program Finished" << endl;
    return 0;
}

/* Prints a message to the console welcoming the user and
 * describing the program. */
void Welcome() {
    cout << "Welcome to CS106L GraphViz!" << endl;
    cout << "This program uses a force-directed graph layout algorithm" << endl;
    cout << "to render sleek, snazzy pictures of various graphs." << endl;
    cout << endl;
}

void GetFileName(ifstream &ifs) {
    string filename;
    cout << "Enter filename : ";
    cin.clear();
    while (getline(cin, filename)) {
        ifs.open(filename);

        if (ifs.fail()) {
            cout << "Wrong file name. Enter again: ";
            ifs.clear(); // reset state bits
            continue;
        }
        else {
            return;
        }
    }
    return;
}

void GetThreshold(chrono::milliseconds &threshold) {
    cout << "Enter maximum duration : ";
    string line;
    int val;
    while (true) {
        if (!getline(cin, line)) throw domain_error("cin error");

        istringstream iss(line);
        char remain;
        if ((iss >> val) && !(iss >> remain)) {
            break;
        }
        cout << "Wrong duration. Enter again : ";
    }
    
    threshold = chrono::milliseconds(val);
    return;
}

void CreateGraphByFile(ifstream &ifs, SimpleGraph &graph) {
    int numNodes;
    ifs >> numNodes;

    graph.nodes.resize(numNodes);

    // initialize position of each node
    for (int k = 0; k < numNodes; k++) {
        graph.nodes[k].x = cos(2 * kPi * k / numNodes);
        graph.nodes[k].y = sin(2 * kPi * k / numNodes);
    }
    
    size_t n1, n2;
    while (ifs >> n1 >> n2) {
        graph.edges.push_back({n1, n2});
    }

    DrawGraph(graph);
    ifs.close();
    cout << "Number of nodes : " << numNodes << ", edges : " << graph.edges.size() << endl;
}

void CalculateNextGraphPosition(SimpleGraph &graph) {
    int numNodes = graph.nodes.size();
    vector<Node> a_delta(numNodes, {0, 0}), r_delta(numNodes, {0, 0}); 

    // compute attractive force between edges
    for (auto [n0Idx, n1Idx] : graph.edges) {
        auto node0 = graph.nodes[n0Idx], node1 = graph.nodes[n1Idx];
        double x0 = node0.x, y0 = node0.y, x1 = node1.x, y1 = node1.y; 

        double fAttract = kAttract * ((y1 - y0) * (y1 - y0) + (x1 - x0) * (x1 - x0));
        double theta = atan2(y1 - y0, x1 - x0);

        a_delta[n0Idx].x += (fAttract * cos(theta));    
        a_delta[n0Idx].y += (fAttract * sin(theta));    
        a_delta[n1Idx].x -= (fAttract * cos(theta));    
        a_delta[n1Idx].y -= (fAttract * sin(theta));    
    }

    // compute repulsive force between every pair of nodes
    for (int n0Idx = 0; n0Idx < numNodes; n0Idx++) {
        for (int n1Idx = n0Idx + 1; n1Idx < numNodes; n1Idx++) {
            auto node0 = graph.nodes[n0Idx], node1 = graph.nodes[n1Idx];
            double x0 = node0.x, y0 = node0.y, x1 = node1.x, y1 = node1.y; 

            double fRepel = kRepel / (double)sqrt((y1 - y0) * (y1 - y0) + (x1 - x0) * (x1 - x0));
            double theta = atan2(y1 - y0, x1 - x0);

            r_delta[n0Idx].x -= (fRepel * cos(theta));
            r_delta[n0Idx].y -= (fRepel * sin(theta));
            r_delta[n1Idx].x += (fRepel * cos(theta));
            r_delta[n1Idx].y += (fRepel * sin(theta));
        }
    }

    for (int i = 0; i < numNodes; i++) {
        graph.nodes[i].x += a_delta[i].x + r_delta[i].x;
        graph.nodes[i].y += a_delta[i].y + r_delta[i].y;
    }
}

void IterateSimpleGraph(chrono::milliseconds threshold, SimpleGraph &graph) {
    cout << "Start iterating graph, given threshold : " << threshold.count() << endl;
    auto start = chrono::high_resolution_clock::now();
    auto cur = chrono::high_resolution_clock::now();

    auto curMs = chrono::duration_cast<chrono::milliseconds>(cur - start);

    while (curMs < threshold) {
        CalculateNextGraphPosition(graph); 
        DrawGraph(graph);

        cur = chrono::high_resolution_clock::now();
        curMs = chrono::duration_cast<chrono::milliseconds>(cur - start);
    }
    
    cur = chrono::high_resolution_clock::now();
    int ms = chrono::duration_cast<chrono::milliseconds>(cur - start).count();

    cout << "Finish iterating graph, took " << ms << " ms " << endl;
}

bool CheckRunAgain() {
    cout << "\nPress \'y\' or \'Y\' to run again, others to finish : ";
    string line;

    while (true) {
        if (!getline(cin, line)) throw domain_error("cin error");

        istringstream iss(line);
        char isEnd, remain;
        if ((iss >> isEnd) && !(iss >> remain)) {
            return (isEnd == 'y') || (isEnd == 'Y');
        }
        cout << "Invalid command. Enter only \'y\' or \'Y\' : ";
    }
    return false;
}
