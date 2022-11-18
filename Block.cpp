#include <Rcpp.h>
#include <algorithm>
#include "Settings.hpp"
#include "Block.hpp"
#include "GridCell.hpp"
#include "misc.hpp"

int Block::m_xmax;
int Block::m_ymax;
int Block::m_ncell;
Settings* Block::sim {};

void Block::addPatch(int cellnum) {
    m_patches.push_back(cellnum);
    ++m_patchcount;
}
void Block::addCell(int cellnum) {

    GridCell &cell { sim->Grid[cellnum] };
    m_cells.push_back(cellnum);
    ++m_cellcount;
    
    for (int season = 0; season <= 3; ++season) {
        // we don't need to multiply by maxU or divide by mean maxent value, because all values are treated equally in each season
        m_total_food[season] += cell.maxent(season);
        //if (sm > m_max_food[season]) m_max_food[season] = sm;
    }
    
    if (cell.Id <= m_ncell) {
        float row = floor(cell.Id / m_xmax);
        float col = cell.Id - row * m_xmax;
        float x = col + 0.5;
        float y = m_ymax - (row + 0.5);
        
        if (m_cellcount == 1) {
            m_block_xmin = x;
            m_block_xmax = x;
            m_block_ymin = y;
            m_block_ymax = y;
        } else {
            m_block_xmin = std::min(x, m_block_xmin);
            m_block_xmax = std::max(x, m_block_xmax);
            m_block_ymin = std::min(y, m_block_ymin);
            m_block_ymax = std::max(y, m_block_ymax);
        }
    }
    
}
void Block::setId(int id) { 
    m_id = id; 
}

void Block::calcCenter() {
    if (m_cellcount > 0) {
        m_center.x = m_block_xmin + (m_block_xmax - m_block_xmin) / 2;
        m_center.y = m_block_ymin + (m_block_ymax - m_block_ymin) / 2;
    }
}
void Block::calcDensity() {
    for (int season = 0; season <= 3; ++season) {
        m_density[season] = (m_total_food[season] / m_cellcount );
        // update 01.10.2022: change from using abundances, only use food value
        m_value[season] = m_density[season] / GridCell::meanSeasonalMaxent[season];
    }
}
void Block::calcValue() {
    int divisor = m_N;
    if (m_N == 0) divisor = 1;
    for (int season = 0; season <= 3; ++season) {
        m_value[season] = m_density[season] / divisor;
    }
}

void Block::setN(int N) {
    m_N = N;
}

int Block::randomCell() const {
    if (m_cellcount > 0) {
        return (*select_randomly(m_cells.begin(), m_cells.end()));
    } else {
        return -1;
    }
}

int Block::id() const {
    return m_id;
}

Vector2df Block::center() const {
    return m_center;
}

float Block::density(int season) const {
    return m_density[season];
}
float Block::value(int season) const {
    return m_value[season];
}
/*
float Block::relativeValue(int season) const {
    return m_relativeValue[season];
}
 */