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
#ifndef ARSGRAPH_H
#define ARSGRAPH_H

#include <memory>
#include <vector>

#include <ars/definitions.h>
#include <ars/functions.h>

namespace ars {

// --------------------------------------------------------
// ARS GRAPH
// --------------------------------------------------------

class ArsGraph {
   public:
    using Self = ArsGraph;
    using Ptr = std::shared_ptr<Self>;

    struct Node {
        std::vector<double> coeffs;
        std::vector<int> incidents; 
    };

    struct Edge {
        int isrc;
        int idst;
        double weight;
        std::vector<double> coeffs;
    };

    /**
     * @brief Constructor of ArsGraph.
     * @param nodeNum number of angle state variables.
     */
    ArsGraph();

    /**
     * @brief Destructor.
     */
    virtual ~ArsGraph();

    /**
     * @brief Returns the number of nodes aka of the state variables.
     */
    size_t getNodeNum() const;

    /**
     * @brief Returns the number of edges aka of the constrains between state
     * variables.
     */
    size_t getEdgeNum() const;

    size_t getFourierOrder() const;

    void setFourierOrder(size_t order);

    int addNode(const std::vector<double>& coeffs);

    int addEdge(int isrc, int idst, double weight = 1.0);

    const std::vector<Node>& nodes() const;

    const std::vector<Edge>& edges() const;

    size_t size() const;

   private:
    std::vector<Node> nodes_;
    std::vector<Edge> edges_;
    size_t fourierOrder_;
};

}  // namespace ars

#endif
