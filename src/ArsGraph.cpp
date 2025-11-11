/**
 * ARS - Angular Radon Spectrum
 * Copyright (C) 2025 Dario Lodi Rizzini - Ernesto Fontana.
 *
 * ARS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ARS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ARS.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <ars/ArsGraph.h>
#include <ars/ars2d.h>

namespace ars {

// --------------------------------------------------------
// ARS GRAPH
// --------------------------------------------------------

ArsGraph::ArsGraph() : nodes_(), edges_() {}

ArsGraph::~ArsGraph() {}

size_t ArsGraph::getNodeNum() const {
    return nodes_.size();
}

size_t ArsGraph::getEdgeNum() const {
    return edges_.size();
}

size_t ArsGraph::getFourierOrder() const {
    return fourierOrder_;
}

void ArsGraph::setFourierOrder(size_t order) {
    fourierOrder_ = order;
}

int ArsGraph::addNode(const std::vector<double>& coeffs) {
    if (coeffs.size() == 2 * (fourierOrder_ + 1)) {
        Node node;
        node.coeffs = coeffs;
        nodes_.push_back(node);
        return ((int)nodes_.size() - 1);
    } else {
        ARS_ERROR("invalid Fourier order: there are "
                  << coeffs.size() << " coefficients instead of "
                  << 2 * (fourierOrder_ + 1));
        return -1;
    }
}

int ArsGraph::addEdge(int isrc, int idst, double weight) {
    int n = getNodeNum();
    if (isrc < 0 || isrc >= n || idst < 0 || idst >= n) {
        ARS_ERROR("invalid node indices: isrc "
                  << "isrc" << isrc << " idst " << idst
                  << " must be in interval [0, " << n << "[");
        return -1;
    }

    Edge edge;
    edge.isrc = isrc;
    edge.idst = idst;
    edge.weight = weight;
    computeFourierCorr(nodes_[isrc].coeffs, nodes_[idst].coeffs, edge.coeffs);
    edges_.push_back(edge);

    nodes_[isrc].incidents.push_back((int)edges_.size() - 1);
    nodes_[idst].incidents.push_back((int)edges_.size() - 1);

    return ((int)edges_.size() - 1);
}

int ArsGraph::addEdgeWithDerivative(int isrc, int idst, double weight) {
    int n = getNodeNum();
    if (isrc < 0 || isrc >= n || idst < 0 || idst >= n) {
        ARS_ERROR("invalid node indices: isrc "
                  << "isrc" << isrc << " idst " << idst
                  << " must be in interval [0, " << n << "[");
        return -1;
    }

    Edge edge;
    edge.isrc = isrc;
    edge.idst = idst;
    edge.weight = weight;
    computeFourierCorr(nodes_[isrc].coeffs, nodes_[idst].coeffs, edge.coeffs);
    //Finding roots of the derivative for each edge and its value in the original function
    std::vector<double> dCoeffs, rootsTheta;
    std::vector<StationaryPoint> sPoints;
    fourierDerivative(edge.coeffs, dCoeffs);
    fourierRootsCCM(dCoeffs, rootsTheta);
    for(int i = 0; i < rootsTheta.size(); i++){
        //since the ars function is periodic on PI and not on 2*PI, the angle found by
        //the CCM method is double the actual angle
        double theta = rootsTheta[i]*0.5;
        sPoints.push_back(StationaryPoint(
            theta, evaluateFourier(edge.coeffs, 2*theta)));
    }
    //Stationary points sorted in descending order considering their value
    std::sort(sPoints.begin(), sPoints.end(), statPointSorter);
    edge.sPoints = sPoints;
    
    edges_.push_back(edge);

    nodes_[isrc].incidents.push_back((int)edges_.size() - 1);
    nodes_[idst].incidents.push_back((int)edges_.size() - 1);

    return ((int)edges_.size() - 1);
}

const std::vector<ArsGraph::Node>& ArsGraph::nodes() const {
    return nodes_;
}

const std::vector<ArsGraph::Edge>& ArsGraph::edges() const {
    return edges_;
}

size_t ArsGraph::size() const{
    size_t sizeNodes = 0;
    for(auto& elem : nodes_){
        sizeNodes += elem.coeffs.size()*sizeof(elem.coeffs.front())
            + elem.incidents.size()*sizeof(elem.incidents.front());
    }
    size_t sizeEdges = 0;
    for(auto& elem : edges_){
        sizeEdges += sizeof(elem.idst) + sizeof(elem.isrc) + sizeof(elem.weight) 
            + elem.coeffs.size()*sizeof(elem.coeffs.front());
    }
    size_t sizeOrder = sizeof(fourierOrder_);
    return sizeNodes + sizeEdges + sizeOrder;
}

void ArsGraph::fourierDerivative(const std::vector<double>& coeffs, std::vector<double>& dCoeffs){
    int fourierOrder = (coeffs.size() - 2) * 0.5;
    dCoeffs.resize(coeffs.size());
    //The average value of the function doesn't influence it's derivative
    dCoeffs[0] = 0;
    dCoeffs[1] = 0;
    //in order to compute the derivative of the fourier series for each level k we have
    //a*cos(k*theta) + b*sin(k*theta) => k*b*cos(k*theta) - k*a*sin(k*theta)
    for(int i = 1; i < fourierOrder+1; i++){
        int idx = 2*i;
        dCoeffs[idx] = i*coeffs[idx+1];
        dCoeffs[idx+1] = -(i*coeffs[idx]);
    }
}
//Method CCM described in the paper https://www.sciencedirect.com/science/article/pii/S0021999113002003#s0070
void ArsGraph::fourierRootsCCM(const std::vector<double>& coeffs, std::vector<double>& roots){
    int fOrder = (coeffs.size() - 2) * 0.5;
    Eigen::MatrixXcd m(2*fOrder, 2*fOrder);
    std::vector<compT> h(2*fOrder + 1);
    roots.clear();

    for (int k = 0; k < h.size(); k++){
        if (k < fOrder)
            h[k] = compT(coeffs[2*(fOrder-k)], coeffs[2*(fOrder-k) + 1]);
        else if (k == fOrder)
            h[k] = 2*coeffs[0];
        else
            h[k] = compT(coeffs[2*(k-fOrder)], -coeffs[2*(k-fOrder) + 1]);
    }

    for (int i = 0; i<2*fOrder-1; i++){
        for (int j = 0; j<2*fOrder; j++){
            if(i == j-1)    m(i,j) = 1.0;
            else            m(i,j) = .0;
        }
    }

    int i = 2*fOrder-1;
    compT den = 1.0/compT(coeffs[2*fOrder], -coeffs[2*fOrder+1]);
    for(int j = 0; j<2*fOrder; j++){
        m(i,j) = -(h[j]*den);
    }

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eigensolver(m);
    std::complex<double> eig;
    double norm;
    for(int i = 0; i < m.rows(); ++i){
        eig = eigensolver.eigenvalues().col(0)[i];
        norm = std::norm(eig);
        if(1.0 - e <= norm && norm <= 1.0 + e){
            double arg = std::arg(eig);
            double root = arg < 0 ? arg + 2*M_PI : arg;
            roots.push_back(root);
        }
    }
}
} // namespace ars
