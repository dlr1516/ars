#ifndef VECTORTREE_H
#define VECTORTREE_H

#include <iostream>
#include <vector>

namespace ars {

    /**
     * Class VectorTree is a self-contained simple implementation of a data structure
     * for storing a D-ary tree (i.e. a tree whose nodes have D children) into 
     * a vector or array. 
     * The tree is organized so that:
     * 1) every node with id N has D children nodes in N*D, N*D+1, ..., N*D + D - 1;
     * 2) the parent of node with id N is node with id floor(N / D). 
     * 
     */
    template <unsigned int D, typename Item, typename Container = std::vector<Item> >
    class VectorTree {
    public:
        using NodeId = int;
        
        /**
         * No node id. 
         */
        static const NodeId INVALID = -1;

        VectorTree() : data_(), nodeNum_(0) {
        }

        ~VectorTree() {
        }
        
        /**
         * Returns the current size of internal data structure. 
         * @return 
         */
        size_t size() const {
            return data_.size();
        }
        
        /**
         * Returns the number of nodes stored inside the tree. 
         * @return 
         */
        size_t nodeNum() const {
            return nodeNum_;
        }
        
        
        NodeId add(NodeId par, unsigned int nchild, const Item& item) {
            NodeId nid;
            
            if (data_.empty()) {
                data_.insert(data_.end(), item);
                return 0;
            }
            else if (par < 0 || par >= data_.size()) {
                std::cerr << __FILE__ << "," << __LINE__ << ": invalid parent node " << par 
                        << ": shoould be in interval [0, " << (data_.size()-1) << "]" << std::endl;
                return INVALID;
            }
            else if (nchild < 0 || nchild >= D) {
                std::cerr << __FILE__ << "," << __LINE__ << ": invalid child index " << nchild
                        << ": shoould be in interval [0, " << (D-1) << "]" << std::endl;
                return INVALID;
            }
            
            nid = D * par + nchild;
            if (nid >= data_.size()) {
                data_.size()
            }
            
            return ((NodeId)data_.size() - 1);
        }
        
        /**
         * Returns the node id of the parent node of the given node. 
         * The returned index is -1 if 
         * @param nid the input node id
         * @return 
         */
        NodeId parent(NodeId nid) const {
            if (nid == )
        }

        /**
         * Returns the nchild-th child node
         * @param nid the id of the node 
         * @param nchild the index of the child 0, 1, ..., D-1
         * @return 
         */
        NodeId child(NodeId nid, unsigned int nchild) const;

        
    private:
        Container data_;
        size_t nodeNum_;
    };
}

#endif /* VECTORTREE_H */

