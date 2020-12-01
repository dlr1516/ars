#include <iostream>
#include <ars/MortonOrderedPoints.h>
#include <ars/definitions.h>

int main(int argc, char** argv) {
    using Morton = ars::MortonOrderedPoints2d4h;
    using SimpleIndex = Morton::SimpleIndex;
    using MultiIndex = Morton::MultiIndex;
    using ArrayIndex = Morton::ArraySimpleIndex;
    Morton mop;
    ArrayIndex ai;
    MultiIndex mi;
    Morton::ConstIterator beg, end;
    ars::VectorVector2 points;
    ars::Vector2 p, octMin, octMax;

    ai[0] = SimpleIndex("0110");
    ai[1] = SimpleIndex("1011");
    mi = Morton::encode(ai);

    ARS_VARIABLE2(ai[0], ai[1]);
    //ARS_VARIABLE2(Morton::enlarge(ai[0]), Morton::enlarge(ai[1]))
    std::cout << "Morton encoded: mi " << mi << std::endl;

    ai[0].reset();
    ai[1].reset();
    std::cout << "Cleared ai[]:\n";
    ARS_VARIABLE2(ai[0], ai[1]);
    std::cout << "After ai = Morton::decode(mi):\n";
    ai = Morton::decode(mi);
    ARS_VARIABLE2(ai[0], ai[1]);
    ARS_VARIABLE(mi);
    
    ARS_VARIABLE(ars::lessBitsets<Morton::MULTI_INDEX_BITNUM>(mi, mi));

    // Insert test points
    p << -2.0, -2.0;
    points.push_back(p);
    p << -1.5, -1.8;
    points.push_back(p);
    p << -1.5, -1.5;
    points.push_back(p);
    p << 0.5, 0.3;
    points.push_back(p);
    p << 0.1, 0.1;
    points.push_back(p);
    p << 0.1, 0.7;
    points.push_back(p);
    p << 1.6, -0.4;
    points.push_back(p);
    p << 2.0, -0.8;
    points.push_back(p);
    p << -2.0, 2.0;
    points.push_back(p);
    p << 1.5, 1.8;
    points.push_back(p);
    p << 1.8, 1.8;
    points.push_back(p);
    
    mop.insert(points);
    
    std::cout << "\nTest point set:\n";
//    for (auto& pi : points) {
//        std::cout << "  [" << pi.transpose() << "] -> " << mop.pointToMorton(pi) << "\n";
//    }

    std::cout << "\nVisit levels with levelMax " << mop.getLevelMax() << std::endl;
    for (int level = 0; level < mop.getLevelMax() && level < 3; ++level) {
        std::cout << "level " << level << ", octant num " << mop.getOctantNum(level) << ":\n";
        for (int octant = 0; octant < mop.getOctantNum(level); ++octant) {
            mop.getOctantBounds(level, octant, octMin, octMax);
            std::cout << "  octant " << octant << " / " << mop.getOctantNum(level) 
                    << ", bounds [" << octMin.transpose() << "][" << octMax.transpose() << "]\n";
            mop.getOctantPoints(level, octant, beg, end);
            for (auto it = beg; it != end; ++it) {
                std::cout << "  [" << it->second.transpose() << "] code "  
                        << it->first << " (check " << mop.pointToMorton(it->second) << ")"<< std::endl;
            }
        }
    }
    
    
    

    return 0;
}
