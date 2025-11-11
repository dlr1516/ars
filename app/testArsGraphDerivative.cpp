#include <iostream>
#include <fstream>

#include <ars/definitions.h>
#include <ars/ars2d.h>
#include <ars/ArsGraph.h>
#include <ars/thirdparty/gnuplot-iostream.h>

using namespace std;

#define RAD2DEG(X) (180.0/M_PI*(X))

int csvToCoeffs(const string&, vector<double>&);
int readConfig(const string& filename, vector<string>& cloudFiles, vector<vector<int>>& edges);

int main(int argc, char** argv){
    vector<vector<int>> edges;
    vector<string> clouds;
    std::chrono::system_clock::time_point timeStart, timeStop;
    double time, timeD;
    ars::ArsGraph::Ptr graph(new ars::ArsGraph);
    ars::ArsGraph::Ptr graphD(new ars::ArsGraph);
    graph->setFourierOrder(30);
    graphD->setFourierOrder(30);

    if(argc == 1 || readConfig(argv[1], clouds, edges) == -1){
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
        graphD->addNode(coeffs);
        nodeCount++;
    }
    std::cout << "Added " << nodeCount << " nodes" << endl;

    int edgeCount = 0;
    timeStart = std::chrono::system_clock::now();
    for(auto& edge : edges){
        graph->addEdge(edge[0], edge[1]);
        edgeCount++; 
    }
    timeStop = std::chrono::system_clock::now();
    time = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    std::cout << "Added " << edgeCount << " edges without derivative in time: " 
        << time << "ms" << std::endl;

    edgeCount = 0;
    timeStart = std::chrono::system_clock::now();
    for(auto& edge : edges){
        graphD->addEdgeWithDerivative(edge[0], edge[1]);
        edgeCount++; 
    }
    timeStop = std::chrono::system_clock::now();
    timeD = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    std::cout << "Added " << edgeCount << " edges with derivative in time: " 
        << timeD << "ms" << std::endl;

    std::cout << "method with derivate takes " << timeD-time << "ms longer" << std::endl;

    int id = 100;
    int thnum = 360;
    for(auto& edge : graphD->edges()){
        Gnuplot gp("gnuplot -persist");
        gp << "set term wxt " << ++id << "\n";
        
        gp << "plot '-' title \"fourier\" w l, " 
            "'-' title \"stationary Points\" w p pt 7 ps 0.5\n";

        for (int i = 0; i < thnum; ++i) {
            double theta = (M_PI / thnum) * i;
            double fourier = ars::evaluateFourier(edge.coeffs, 2.0 * theta);
            gp << RAD2DEG(theta) << " " << fourier << "\n";
        }
        gp << "e" << std::endl;

        for(size_t i = 0; i < edge.sPoints.size(); i++){
            double theta = edge.sPoints[i].theta;
            double fourier = edge.sPoints[i].val;
            std::cout << "id " << id << ": " << RAD2DEG(theta) << "," << fourier << "; ";
            gp << RAD2DEG(theta) << " " << fourier << "\n";
        }
        gp << "e" << std::endl;
        std::cout << std::endl;
    }

    return 0;
}

int readConfig(const string& filename, vector<string>& cloudFiles, vector<vector<int>>& edges){
    std::ifstream file(filename);
    if(file.fail()){
       return -1; 
    }
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
    if(file.fail()){
       return -1; 
    }
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