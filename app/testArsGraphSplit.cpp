#include <iostream>
#include <fstream>

#include <ars/definitions.h>
#include <ars/ars2d.h>
#include <ars/ArsGraph.h>
#include <ars/ArsGraphInterval.h>

using namespace std;

int csvToCoeffs(const string&, vector<double>&);
int readConfig(const string& filename, vector<string>& cloudFiles, vector<vector<int>>& edges);

/*const vector<int> clouds {0, 400, 800, 1100};
const vector<vector<int>> edges {{0,1}, {1,2}, {2,3}, {1,3}, {3,0}};*/

int main(int argc, char** argv){
    vector<vector<int>> edges;
    vector<string> clouds;
    ars::ArsGraph::Ptr graph(new ars::ArsGraph);
    graph->setFourierOrder(30);
    ars::ArsGraphIntervalFull::Ptr interval(new ars::ArsGraphIntervalFull);

    if(argc == 0 || readConfig(argv[1], clouds, edges) == -1){
        std::cout << "Couldn't read config file!" << std::endl;
        return -1;
    }

    int nodeCount = 0;
    for(auto& cloud : clouds){
        vector<double> coeffs;

        if(csvToCoeffs(cloud, coeffs) == -1){
            std::cout << "Invalid cloud!" << endl;
            return -1;
        }
        graph->addNode(coeffs);
        nodeCount++;
    }
    std::cout << "Added " << nodeCount << " nodes" << endl;

    int edgeCount = 0;
    for(auto& edge : edges){
        graph->addEdge(edge[0], edge[1]);
        edgeCount++; 
    }
    std::cout << "Added " << edgeCount << " edges" << endl;

    interval->init(graph);

    std::cout << "State before split: " << endl;
    for(int i = 0; i < clouds.size(); i++){
        std::cout << "node " << i << ": lower theta " << interval->nodeLower(i) << 
            " upper theta " << interval->nodeUpper(i) << endl; 
    }
    for(int i = 0; i < edges.size(); i++){
        std::cout << "edge " << edges[i][0]  << "-" << edges[i][1] << ": lower fourier " << 
            interval->edgeLower(i) << " upper fourier " << interval->edgeUpper(i) << endl; 
    }

    std::cout << endl << endl << "State after splits: " << endl;
    for(int i = 1; i < clouds.size(); i++){
        std::cout << "split on node " << i << endl;

        ars::ArsGraphIntervalFull::Ptr intervLower(new ars::ArsGraphIntervalFull);
        ars::ArsGraphIntervalFull::Ptr intervUpper(new ars::ArsGraphIntervalFull);
        interval->split(i, intervLower, intervUpper);
        std::cout << "nodes: " << std::endl;
        for(int j = 0; j < clouds.size(); j++){
            std::cout << "node " << j <<": lower: " << "lower theta " << intervLower->nodeLower(j) << 
                " upper theta " << intervLower->nodeUpper(j) << endl; 

            std::cout << "node " << j <<": upper: " << "lower theta " << intervUpper->nodeLower(j) << 
                " upper theta " << intervUpper->nodeUpper(j) << endl; 
        }
        std::cout << "edges: " << endl;
        for(int j = 0; j < edges.size(); j++){
            std::cout << "edge " << edges[j][0]  << "-" << edges[j][1] << " lower: lower fourier " << 
                intervLower->edgeLower(j) << " upper fourier " << intervLower->edgeUpper(j) << endl;

            std::cout << "edge " << edges[j][0]  << "-" << edges[j][1] << " upper: lower fourier " << 
                intervUpper->edgeLower(j) << " upper fourier " << intervUpper->edgeUpper(j) << endl; 
        }
        std::cout << endl << endl;
    }

    return 0;
}

int readConfig(const string& filename, vector<string>& cloudFiles, vector<vector<int>>& edges){
    std::ifstream file(filename);
    std::string line;
    cloudFiles.clear();
    edges.clear();

    std::getline(file, line);
    string tmp;
    {
        std::stringstream ssline(line);
        while (std::getline(ssline, tmp, ' ')) {
            cloudFiles.push_back(tmp);
        }
    }

    std::getline(file, line);
    {
    std::stringstream ssline(line);
        while (std::getline(ssline, tmp, ' ')) {
            vector<int> tmpEdge(2);
            int pos = tmp.find(',');
            if (pos != std::string::npos) {
                tmp.replace(pos, 1, " ");
            }
            std::stringstream ssEdge(tmp);
            if (ssEdge >> tmpEdge[0] >> tmpEdge[1] ) {
                edges.push_back(tmpEdge);
            }
        }
    }
    file.close();
    if(cloudFiles.size() == 0 || edges.size() == 0)
        return -1;
    else
        return 0;
}

int csvToCoeffs(const std::string& filename, std::vector<double>& coeffs){
    std::ifstream file(filename);
    std::string line;
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