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
    if (coeffs.size() == 2 * fourierOrder_) {
        Node node;
        node.coeffs = coeffs;
        nodes_.push_back(node);
        return ((int)nodes_.size() - 1);
    } else {
        ARS_ERROR("invalid Fourier order: there are "
                  << coeffs.size() << " coefficients instead of "
                  << (2 * fourierOrder_));
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

    return ((int)edges_.size() - 1);
}

}  // namespace ars
