#ifndef __BLOCK__
#define __BLOCK__

#include "Settings.hpp"
#include "GridCell.hpp"
#include "Vector2d.hpp"
#include "misc.hpp"

class GridCell;
class Settings;

class Block {
private:
    static int m_xmax, m_ymax, m_ncell;
    static Settings* sim;
    int m_id {}; // id = index in vector of blocks
    Vector2df m_center; // block center coordinates, calculated once all traversable cells are known
    std::vector<int> m_cells; // traversable cells in block
    std::vector<int> m_patches; // food patches in block
    int m_cellcount {}; // number of traversable cells in block
    int m_patchcount {}; // number of food patches in block
    float m_total_food[4] {}; // total food per season
    float m_density[4] {}; // total food per season / number of cells
    float m_value[4] {}; // density / number of porps in block
    //float m_relativeValue[4] {}; // value / best value among all blocks
    float m_block_xmin, m_block_xmax, m_block_ymin, m_block_ymax;
    void setId(int id);
    friend Settings;
    
public:
    //static float bestBlockValue[4];
    Block() = default;
    ~Block() = default;
    int m_N { 1 }; // number of porpoises starting out in this block, used to calculate perceived block value
    void addPatch(int cellnum);
    void addCell(int cellnum);
    void calcCenter();
    void calcDensity();
    void calcValue();
    void setN(int N);
    int randomCell() const;
    Vector2df center() const;
    int id() const;
    float density(int season) const;
    float value(int season) const;
};


#endif // __BLOCK__
