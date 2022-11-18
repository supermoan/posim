#include <Rcpp.h>
#include "Gillnet.h"
#include "Vector2d.hpp"
#include "Settings.hpp"
#include "misc.hpp"
#include "Logger.h"

// todo: handle tries >= 10 in posOK routine better

std::vector<float> Gillnet::interaction_probability;
float Gillnet::m_catchability_by_type[3];

Settings* Gillnet::sim;

Gillnet::Gillnet(int block, int type, int soaktime, float length, bool pingered) {
    m_block = block;
    m_type = type;
    m_max_soaktime = soaktime;
    m_length = length;
    m_pingered = pingered;
    
    m_catchability = m_catchability_by_type[type];
    if (pingered == true) {
        m_catchability *= sim->pinger_effect;
    }

    Vector2df start;
    Vector2df end;

    bool posOK = false;
    int tries = 0;
    
    while (posOK == false && tries < 10) {
        // placement of gillnets: pick a random cell within the fishery block that was specified,
        // randomly pick xy coordinates for one end of the gillnet string within that cell, then 
        // pick a random angle, and find the xy coordinates of the point that lies a distance m_length
        // grid away in the direction determined by the angle. Once we have the coordinates of the 
        // start and end of our line segment, check if any part of it overlaps with an intraversable
        // cell (i.e. a cell on land). Repeat until we find a line segment that does not have such overlap.
        int cellnum = *select_randomly(sim->FisheryBlocks[m_block].begin(), sim->FisheryBlocks[m_block].end());
        float theta = getRandomFloat(0, 1) * 2*PI; // random angle
        start = sim->pointFromCell(cellnum); // coordinates of cell center
        start.x += getRandomFloat(-0.49, 0.49); // offset x from center by random distance
        start.y += getRandomFloat(-0.49, 0.49); // offset y from center by random distance
        // calculate xy coordinates for the end of our segment
        end = Vector2df(start.x + m_length * cos(theta), start.y + m_length * sin(theta));
        //Logger::debug(0, "gillnet: theta = %f, m_length = %f, x1 = %f, y1 = %f", theta, m_length, end.x, end.y);
        
        // list of cells intersected by gillnet
        m_cellnums = sim->GetCellsIntersected(start.x, start.y, end.x, end.y);
        
        // if both coordinates are within the map extent,
        // then start ut assuming all cells are traversable (i.e. in water)
        posOK = sim->isCoordValid(start) && sim->isCoordValid(end);
        
        if (posOK) {
            // then check the depth of each cell, and if we find at least one that is not deep
            // enough, generate a new position for the gillnet
            for (const auto& cell : m_cellnums) {
                if (cell == -1 || sim->Grid[cell].Bathymetry > sim->MinimumWaterDepth) {
                    posOK = false;
                }
                break;
            }
        }
        ++tries;
    }
    
    m_coords.first = start;
    m_coords.second = end;
    sort(m_cellnums.begin(), m_cellnums.end());
    
    for (const auto& cell : m_cellnums) {
        sim->Grid[cell].gillnets.push_front(this);
        m_cells.push_back(sim->Grid[cell].gillnets.begin());
    }
}

Gillnet::~Gillnet() {
    for (int i = 0; i < m_cellnums.size(); ++i) {
        sim->Grid[i].gillnets.erase(m_cells[i]);
    }
}

// check: when this gets called, that means a porpoise has entered a cell
// with a gillnet. This function calculates the squared distance of the 
// porpoise's current position to the nearest point of the gillnet. 
// if the square distance is less than 0.015625, which is 50 meters squared
// in map units, then we calculate an interaction probability. 
// Finally, we draw a real number from unif(0,1) and check that against the
// entanglement probability, which is catchability of gillnet type 
// multiplied with interaction probability
bool Gillnet::check4(Vector2df xy, int cell) {
    float dist = distanceSquared(m_coords, xy);
    // if distance is closer than 50 meters 
    // (in map units, (50m)^2 = (1/400*50)^2)
    if (dist > 0.015625f) return false;
    // then proceed to do further checks
    dist = sqrt(dist); // calculate distance
    int index = round(dist * 800) * 0.5; // find index in lookup table of probabilities
    float prob_interaction = interaction_probability[index];
    float prob_entangled = prob_interaction * m_catchability;
    
    // finally check if porp was entangled
    if (getRandomFloat(0, 1) < prob_entangled) {
        m_catch++;
        return true;
    }
    return false;
}
bool Gillnet::check3(Vector2df xy, int cell) {
    
    // if porp's current cell is smaller or larger than the smallest and largest
    // cell numbers intersected by this gillnet, then there's no chance of interaction
    if (cell < m_cellnums[0] || cell > m_cellnums.back()) {
        return false;
    }
    // otherwise, there could be an interaction. check to see whether the porp's
    // current cell intersects with any of the gillnet cells
    if (std::find(m_cellnums.begin(), m_cellnums.end(), cell) != m_cellnums.end()) {
        // if it does, calculate the minimum squared distance from porp's pos to the gillnet string
        float dist = distanceSquared(m_coords, xy);
        
        // if distance > 50 meters, stop here.
        if (dist >= 0.015625) return false;
        
        // Distance is closer than 50 meters (in map units, (50m)^2 = (1/400*50)^2)
        // Proceed to do further checks
        dist = sqrt(dist); // calculate distance
        int index = round(dist * 800); // find index in lookup table of probabilities
        const float& prob_interaction = interaction_probability[index];
        float prob_entangled = prob_interaction * m_catchability;
        
        // finally check if porp was entangled
        if (getRandomFloat(0, 1) < prob_entangled) {
            m_catch++;
            return true;
        }
    }
    return false;
}
bool Gillnet::check2(std::pair<Vector2df, Vector2df> path) {
    /*
    //float dist2 = distanceSquared(this->m_coords, path); // distance in map units

    //if (dist2 > 0.015625) { // if distance is > 50 meters
        return false;
    }
    
    int index = round(dist2 * 800) * 0.5;
    
    //float maxp { 0.0531923f };
    //float prob_interaction = R::dnorm4(dist, 0.0, 7.5, false) / maxp;
    
    float prob_interaction = interaction_probability[index];
    float prob_entangled = prob_interaction * m_catchability;
     */
    
    // check if porp's latest position corresponds to a cell
    
    if (intersects(this->m_coords, path)) {
            
        float prob_entangled = m_catchability;
        
        if (getRandomFloat(0, 1) < prob_entangled) {
            m_catch++;
            return true;
        }
    }
    return false;
    
}
bool Gillnet::check(std::pair<Vector2df, Vector2df> path) {
    float dist = distance(this->m_coords, path); // distance in map units
    
    // the probability of getting entangled depends on the catchability of the
    // net type and the distance between the gillnet and porpoise's path of travel
    // the catchability of the gillnet is scaled by a distance factor e^-(400*0.75*distance)
    float p = m_catchability;
    if (dist > 0) p *= pow(exp(1), (400*-0.075 * dist));
    
    if (getRandomFloat(0, 1) < p) {
        m_catch++;
        return true;
    } else {
        return false;
    }
}


