#include <iostream>
#include <fstream>

#include <ars/definitions.h>
#include <ars/ars2d.h>
#include <ars/ArsGraph.h>
#include <ars/ArsGraphInterval.h>

using namespace std;

int csvToCoeffs(const std::string&, std::vector<double>&);

const vector<int> clouds {0, 400, 800, 1100};
const vector<vector<int>> edges {{0,1}, {1,2}, {2,3}, {1,3}, {3,0}};

int main(){
    ars::ArsGraph::Ptr graph(new ars::ArsGraph);
    graph->setFourierOrder(30);
    ars::ArsGraphIntervalFull::Ptr interval(new ars::ArsGraphIntervalFull);

    int nodeCount = 0;
    for(auto& cloud : clouds){
        vector<double> coeffs;

        stringstream ss;
        ss << setw(5) << setfill('0') << cloud;
        string name = "/home/frencio/uni/ARS/dataset/csv_2d/coeff_" + ss.str() + ".csv";

        if(csvToCoeffs(name, coeffs) == -1){
            cout << "Invalid cloud!" << endl;
            return -1;
        }
        graph->addNode(coeffs);
    }
    cout << "Added " << nodeCount << " nodes" << endl;

    int edgeCount = 0;
    for(auto& edge : edges){
        graph->addEdge(edge[0], edge[1]);
        edgeCount++; 
    }
    cout << "Added " << edgeCount << " edges" << endl;

    interval->init(graph);

    cout << "State before split: " << endl;
    for(int i = 0; i < clouds.size(); i++){
        cout << "node " << i << ": lower theta " << interval->stateLower(i) << 
            " upper theta " << interval->stateUpper(i) << endl; 
    }
    for(int i = 0; i < edges.size(); i++){
        cout << "edge " << edges[i][0]  << "-" << edges[i][1] << ": lower fourier " << 
            interval->edgeLower(i) << " upper fourier " << interval->edgeUpper(i) << endl; 
    }

    cout << endl << endl << "State after splits: " << endl;
    for(int i = 1; i < clouds.size(); i++){
        cout << "split on node " << i << endl;
        /*ars::ArsGraphIntervalFull::Ptr intervLower(new ars::ArsGraphIntervalFull(graph));
        ars::ArsGraphIntervalFull::Ptr intervUpper(new ars::ArsGraphIntervalFull(graph));*/
        ars::ArsGraphIntervalFull::Ptr intervLower;
        ars::ArsGraphIntervalFull::Ptr intervUpper;
        interval->split(i, intervLower, intervUpper);
        cout << "nodes: " << std::endl;
        for(int j = 0; j < clouds.size(); j++){
            cout << "node " << j <<": lower: " << "lower theta " << intervLower->stateLower(j) << 
                " upper theta " << intervLower->stateUpper(j) << endl; 

            cout << "node " << j <<": upper: " << "lower theta " << intervUpper->stateLower(j) << 
                " upper theta " << intervUpper->stateUpper(j) << endl; 
        }
        cout << "edges: " << endl;
        for(int j = 0; j < edges.size(); j++){
            cout << "edge " << edges[j][0]  << "-" << edges[j][1] << " lower: lower fourier " << 
                intervLower->edgeLower(j) << " upper fourier " << intervLower->edgeUpper(j) << endl;

            cout << "edge " << edges[j][0]  << "-" << edges[j][1] << " upper: lower fourier " << 
                intervUpper->edgeLower(j) << " upper fourier " << intervUpper->edgeUpper(j) << endl; 
        }
        cout << endl << endl;
    }

    return 0;
}

int csvToCoeffs(const std::string& filename, std::vector<double>& coeffs){
    std::ifstream file(filename);
    std::string line, comment;
    coeffs.clear();
    while (!file.eof()) {
        double tmp = 0;
        std::getline(file, line, ',');
        std::stringstream ssline(line);
        if (ssline >> tmp) {
            coeffs.push_back(tmp);
        }
    }
    file.close();
    if(coeffs.size() > 0)
        return 0;
    else
        return -1;
}