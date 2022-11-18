#pragma once

#include "Settings.hpp"
#include "Vector2d.hpp"
#include "misc.hpp"

/*
 * The abundance region class is a wrapper around a vector of integers
 * denoting the cell numbers that make up one abundance region. 
 */
class abundanceRegion {
    static Settings* _sim;
    int _id {-1};
    std::vector<int> _cells;
    friend Settings;
public:
    int _births{0};
    int _deaths{0};
    int _bycatch{0};

    abundanceRegion() {};
    void addCell(const int cell);
    bool empty();
    int randomCell();
    Vector2df randomPoint();
    void reset_stats();
    int id();
    void setId(int id);
};

