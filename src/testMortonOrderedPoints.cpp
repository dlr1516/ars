#include <iostream>
#include <ars/MortonOrderedPoints.h>
#include <ars/MortonTree.h>
#include <ars/definitions.h>
#include <deque>

int main(int argc, char** argv) {
    ars::VectorVector2 points;
    ars::Vector2 p, octMin, octMax;

    {
        using Morton = ars::MortonOrderedPoints2d4h;
        using SimpleIndex = Morton::SimpleIndex;
        using MultiIndex = Morton::MultiIndex;
        using ArrayIndex = Morton::ArraySimpleIndex;
        Morton mop;
        ArrayIndex ai;
        MultiIndex mi;
        Morton::ConstIterator beg, end;


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
                            << it->first << " (check " << mop.pointToMorton(it->second) << ")" << std::endl;
                }
            }
        }
    }


    {
        std::cout << "\n\n-----\nTEST MORTON_TREE (different implementation):\n";
        using MyMortonTree = ars::MortonTree<2, ars::Vector2, uint16_t>;
        using MyIntPoint = MyMortonTree::IntPoint;
        using MyLessIntPoint = MyMortonTree::LessIntPoint;
        using MyOctant = MyMortonTree::Octant;
        using MyOctantChildrenArray = MyMortonTree::OctanChildrenArray;

        MyMortonTree mt;
        MyIntPoint ip, ip1, ip2;
        MyLessIntPoint lip;
        std::deque<MyOctant> queue;
        MyMortonTree::ConstIterator octBeg, octEnd;
        MyOctantChildrenArray children;
        float value = 0.0f;

        ip1 << 1, 2;
        ip2 << 3, 0;
        ARS_VARIABLE4(ip1.transpose(), ip2.transpose(), lip(ip1, ip2), lip(ip2, ip1));
        ARS_VARIABLE3(MyMortonTree::INTEGER_BITNUM, MyMortonTree::INTEGER_MAX, MyMortonTree::CHILDREN_NUM);

        for (auto& p : points) {
            ip(0) = floor((p(0) + 2.0) * MyMortonTree::INTEGER_MAX / 4.0);
            if (ip(0) > MyMortonTree::INTEGER_MAX) ip(0) = MyMortonTree::INTEGER_MAX;
            ip(1) = floor((p(1) + 2.0) * MyMortonTree::INTEGER_MAX / 4.0);
            if (ip(1) > MyMortonTree::INTEGER_MAX) ip(1) = MyMortonTree::INTEGER_MAX;

            std::cout << "  inserting [" << p.transpose() << "] -> [" << ip.transpose() << "] p " << p.transpose()
                    << ": mt.size() " << mt.size() << "\n";
            mt.insert(ip, p);
            value += 1.0f;
        }

        std::cout << "Items in MortonTree:\n";
        for (auto it = mt.begin(); it != mt.end(); ++it) {
            std::cout << "  [" << it->first.transpose() << "]: " << it->second.transpose() << "\n";
        }

        std::cout << "\nVisit levels with levelMax " << mt.getLevelMax() << std::endl;
        queue.push_back(mt.root());
        while (!queue.empty()) {
            MyOctant octCur = queue.front();
            queue.pop_front();
            //std::cout << "Octant level " << octCur.level << " mask " << octCur.mask.transpose() << "\n";
            std::cout << "Octant " << octCur.getString() << "\n";
            mt.getOctant(octCur, octBeg, octEnd);
            for (auto it = octBeg; it != octEnd; ++it) {
                std::cout << "   [" << it->first.transpose() << "], [" << it->second.transpose() << "]\n";
            }
            if (octCur.level < 2) {
                children = mt.children(octCur);
                queue.insert(queue.end(), children.begin(), children.end());
            }
        }
    }



    return 0;
}
