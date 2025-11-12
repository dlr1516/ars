#include <iostream>
#include <fstream>

#include <ars/definitions.h>
#include <ars/ars2d.h>
#include <ars/ArsGraph.h>
#include <ars/ArsGraphSolver.h>
#include <ars/thirdparty/gnuplot-iostream.h>

using namespace std;

#define RAD2DEG(X) (180.0/M_PI*(X))

int csvToCoeffs(const string&, vector<double>&);
int readConfig(const string& filename, vector<string>& cloudFiles, vector<vector<int>>& edges);


int main(int argc, char** argv){
    vector<vector<int>> edges;
    vector<string> clouds;
    bool useDiff = false;
    ars::ArsGraph::Ptr graph(new ars::ArsGraph);
    ars::ArsGraph::Ptr graphD(new ars::ArsGraph);
    std::chrono::system_clock::time_point timeStart, timeStop;
    double time, timeS;
    graph->setFourierOrder(30);
    graphD->setFourierOrder(30);


    if(argc == 1 || readConfig(argv[1], clouds, edges) == -1){
        std::cout << "Couldn't read config file!" << std::endl;
        return -1;
    }

    if (argc > 2 && std::string(argv[2]) == "true") {
        useDiff = true;
        std::cout << "Enabled differential Intervals" << std::endl;
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
    
    std::cout << std::endl << "With stationary points:" << std::endl;
    timeStart = std::chrono::system_clock::now();
    {
        ars::ArsGraphIntervalFull::Ptr interval(new ars::ArsGraphIntervalFull);        
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
       //solver.initialSolutionFromTree();

        std::vector<double> solution;
        double cost = 0;
        ars::ArsGraphSolver::Statistics stats;
        solver.solve(solution, cost, stats, useDiff);

        std::cout << "Solution: " << std::endl;
        for(int i = 0; i < solution.size(); i++){
            std::cout << "Node " << i << ": " 
                << RAD2DEG(solution[i]) << std::endl; 
        }

        std::cout << "Created nodes: " << stats.createdNodes << std::endl;
        std::cout << "Min interval size: " << stats.minIntervalSize << std::endl;
        std::cout << "Max interval size: " << stats.maxIntervalSize << std::endl;
        std::cout << "Average interval size: " << stats.avgIntervalSize << std::endl;
    }
    timeStop = std::chrono::system_clock::now();
    time = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    std::cout << "Solver without stationary points ended in time: " 
        << time << "ms" << std::endl;

    std::cout << std::endl << "With stationary points:" << std::endl;
    timeStart = std::chrono::system_clock::now();
    {
        ars::ArsGraphIntervalFull::Ptr interval(new ars::ArsGraphIntervalFull); 
        int edgeCount = 0;
        for(auto& edge : edges){
            graphD->addEdgeWithDerivative(edge[0], edge[1]);
            edgeCount++; 
        }
        std::cout << "Added " << edgeCount << " edges" << endl;

        interval->initWithStationary(graphD);

        std::cout << "Initial state: " << endl;
        for(int i = 0; i < clouds.size(); i++){
            std::cout << "node " << i << ": lower theta " << interval->nodeLower(i) << 
                " upper theta " << interval->nodeUpper(i) << endl; 
        }
        for(int i = 0; i < edges.size(); i++){
            std::cout << "edge " << edges[i][0]  << "-" << edges[i][1] << ": lower fourier " << 
                interval->edgeLower(i) << " upper fourier " << interval->edgeUpper(i) << endl; 
        }

        ars::ArsGraphSolver solver(graphD);
        //solver.initialSolutionFromTree();

        std::vector<double> solution;
        double cost = 0;
        ars::ArsGraphSolver::Statistics stats;
        solver.solveWithStationary(solution, cost, stats, useDiff);

        std::cout << "Solution: " << std::endl;
        for(int i = 0; i < solution.size(); i++){
            std::cout << "Node " << i << ": " 
                << RAD2DEG(solution[i]) << std::endl; 
        }

        std::cout << "Created nodes: " << stats.createdNodes << std::endl;
        std::cout << "Min interval size: " << stats.minIntervalSize << std::endl;
        std::cout << "Max interval size: " << stats.maxIntervalSize << std::endl;
        std::cout << "Average interval size: " << stats.avgIntervalSize << std::endl;
    }
    timeStop = std::chrono::system_clock::now();
    timeS = (double) std::chrono::duration_cast<std::chrono::milliseconds>(timeStop - timeStart).count();
    std::cout << "Solver with stationary points ended in time: " 
        << timeS << "ms" << std::endl;

    std::cout << "Method with stationaty points takes " << time-timeS << "ms less" << std::endl;

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
            gp << RAD2DEG(theta) << " " << fourier << "\n";
        }
        gp << "e" << std::endl;
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