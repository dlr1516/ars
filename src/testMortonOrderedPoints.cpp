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
    
    ai[0] = SimpleIndex("0110");
    ai[1] = SimpleIndex("1011");
    mi = Morton::encode(ai);
    
    ARS_VARIABLE2(ai[0], ai[1]);
    ARS_VARIABLE2(Morton::enlarge(ai[0]), Morton::enlarge(ai[1]))
    std::cout << "Morton encoded: mi " << mi << std::endl;
    
    ai[0].reset();
    ai[1].reset();
    std::cout << "Cleared ai[]:\n";
    ARS_VARIABLE2(ai[0], ai[1]);
    std::cout << "After ai = Morton::decode(mi):\n";
    ai = Morton::decode(mi);
    ARS_VARIABLE2(ai[0], ai[1]);
    ARS_VARIABLE(Morton::truncate(mi));
            
    return 0;
}
