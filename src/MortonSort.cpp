#include <ars/MortonSort.h>

namespace ars {
    
    // ---------------------------------------------------------------
    // NUMBER OF LEADING ZEROS
    // ---------------------------------------------------------------
    
        // ---------------------------------------------------------------
    // NUMBER OF LEADING ZEROS
    // ---------------------------------------------------------------
    
    int8_t nlz8(uint8_t x) {
        uint8_t y;
        int n;

        n = 8;
        y = x >> 4;
        if (y != 0) {
            n = n - 4;
            x = y;
        }
        y = x >> 2;
        if (y != 0) {
            n = n - 2;
            x = y;
        }
        y = x >> 1;
        if (y != 0) return n - 2;
        return n - x;
    }

    int16_t nlz16(uint16_t x) {
        uint16_t y;
        int n;

        n = 16;
        y = x >> 8;
        if (y != 0) {
            n = n - 8;
            x = y;
        }
        y = x >> 4;
        if (y != 0) {
            n = n - 4;
            x = y;
        }
        y = x >> 2;
        if (y != 0) {
            n = n - 2;
            x = y;
        }
        y = x >> 1;
        if (y != 0) return n - 2;
        return n - x;
    }

    int32_t nlz32(uint32_t x) {
        uint32_t y;
        int n;

        n = 32;
        y = x >> 16;
        if (y != 0) {
            n = n - 16;
            x = y;
        }
        y = x >> 8;
        if (y != 0) {
            n = n - 8;
            x = y;
        }
        y = x >> 4;
        if (y != 0) {
            n = n - 4;
            x = y;
        }
        y = x >> 2;
        if (y != 0) {
            n = n - 2;
            x = y;
        }
        y = x >> 1;
        if (y != 0) return n - 2;
        return n - x;
    }

    int64_t nlz64(uint64_t x) {
        uint64_t y;
        int n;

        n = 64;
        y = x >> 32;
        if (y != 0) {
            n = n - 32;
            x = y;
        }
        y = x >> 16;
        if (y != 0) {
            n = n - 16;
            x = y;
        }
        y = x >> 8;
        if (y != 0) {
            n = n - 8;
            x = y;
        }
        y = x >> 4;
        if (y != 0) {
            n = n - 4;
            x = y;
        }
        y = x >> 2;
        if (y != 0) {
            n = n - 2;
            x = y;
        }
        y = x >> 1;
        if (y != 0) return n - 2;
        return n - x;
    }
    
}

