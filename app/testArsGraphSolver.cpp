#include <iostream>
#include <fstream>

#include <ars/definitions.h>
#include <ars/ars2d.h>
#include <ars/ArsGraph.h>
#include <ars/ArsGraphSolver.h>

using namespace std;

int csvToCoeffs(const string&, vector<double>&);
int readConfig(const string& filename, vector<string>& cloudFiles, vector<vector<int>>& edges);

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

    std::cout << "Initial state: " << endl;
    for(int i = 0; i < clouds.size(); i++){
        std::cout << "node " << i << ": lower theta " << interval->nodeLower(i) << 
            " upper theta " << interval->nodeUpper(i) << endl; 
    }
    for(int i = 0; i < edges.size(); i++){
        std::cout << "edge " << edges[i][0]  << "-" << edges[i][1] << ": lower fourier " << 
            interval->edgeLower(i) << " upper fourier " << interval->edgeUpper(i) << endl; 
    }

    ars::ArsGraphSolver solver(graph);

    std::vector<double> solution;
    double cost;
    solver.solve(solution, cost);

    std::cout << "Solution: " << std::endl;
    for(int i = 0; i < solution.size(); i++){
        std::cout << "Node " << i << ": " 
            << solution[i]*(180.0/M_PI) << std::endl; 
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