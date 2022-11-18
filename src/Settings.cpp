#include <Rcpp.h>
#include <random>
#include "Settings.hpp"
#include "misc.hpp"
#include "IO.hpp"
#include "Logger.h"
#include "Vector2d.hpp"
#include "Gillnet.h"
#include "Block.hpp"
#include "Porpoise.hpp"

Settings::Settings(Rcpp::List& conf) {
    using namespace Rcpp;
    
    Porpoise::sim = { this };
    Block::sim = { this };
    Gillnet::sim = { this };
    abundanceRegion::_sim = { this };
    
    nFisheryDataSampleSize= as<int>(conf["nFisheryDataSampleSize"]);
    nGillnetTypes = as<int>(conf["nGillnetTypes"]);
    MinimumWaterDepth = as<float>(conf["minTraversableWaterDepth"]);
    FoodGrowthRate = as<float>(conf["foodGrowthRate"]);
    offGridCellsTraversable = as<bool>(conf["offGridCellsTraversable"]);
    maxU = as<float>(conf["maxU"]);
    pinger_effect = as<float>(conf["pingerEffect"]);
    Porpoise::nextId = 0;
    // set porpoise static members
    Porpoise::AgeOfMaturity = as<float>(conf["ageOfMaturity"]);
    
    NumericVector tmp = conf["monthlyEnergyMultiplier"];
    for (int i = 0; i < 12; i++) {
        Porpoise::monthlyEnergyMultiplier[i] = tmp[i];
    }
    Porpoise::withCalfEnergyMultiplier = as<float>(conf["withCalfEnergyMultiplier"]);
    Porpoise::distEnergyMultiplier = as<float>(conf["distEnergyMultiplier"]);
    Porpoise::stepEnergyMultiplier = as<float>(conf["stepEnergyMultiplier"]);
    Porpoise::m_mort_prob = as<float>(conf["mMortProb"]);
    Porpoise::x_surv_prob = as<float>(conf["xSurvProb"]);
    Porpoise::dispersalInertia = as<int>(conf["dispersalInertia"]) - 1;
    Porpoise::meanDispersalDistance = as<float>(conf["meanDispersalDistance"]);
    Porpoise::minDispersalDistance = as<float>(conf["minDispersalDistance"]);
    Porpoise::maxDispersalDistance = as<float>(conf["maxDispersalDistance"]);
    Porpoise::minDispersalDistanceSquared = Porpoise::minDispersalDistance * Porpoise::minDispersalDistance;
    Porpoise::maxDispersalDistanceSquared = Porpoise::maxDispersalDistance * Porpoise::maxDispersalDistance;
        
    Porpoise::minDispersalDepth = as<float>(conf["minDispersalDepth"]);
    Porpoise::minDispersalDistanceToLand = as<float>(conf["minDispersalDistanceToLand"]);
    Porpoise::pregnancy_prob = as<float>(conf["pregnancyProb"]);
    Porpoise::max_age = as<float>(conf["maxAge"]);
    Porpoise::maxMemory = as<int>(conf["memoryMax"]);
    Porpoise::inertia_const = as<float>(conf["inertiaConst"]);
    Porpoise::corrLogmov = as<float>(conf["corrLogmov"]);
    Porpoise::corrAngle = as<float>(conf["corrAngle"]);
    Porpoise::m = as<float>(conf["m"]);
    Porpoise::maxLogmov = as<float>(conf["maxLogmov"]);
    Porpoise::age_dist = as<std::vector<int>>(conf["ageDist"]);
    Porpoise::ref_mem_strength = as<std::vector<float>>(conf["refMemStrength"]);
    Porpoise::work_mem_strength = as<std::vector<float>>(conf["workMemStrength"]);
    GridCell::foodGrowthRate = as<float>(conf["foodGrowthRate"]);
    Gillnet::interaction_probability = as<std::vector<float>>(conf["interaction_probability"]);
    Gillnet::m_catchability_by_type[0] = as<float>(conf["catchabilitySmall"]);
    Gillnet::m_catchability_by_type[1] = as<float>(conf["catchabilityMedium"]);
    Gillnet::m_catchability_by_type[2] = as<float>(conf["catchabilityLarge"]);
        
    int yday = as<int>(conf["start"]);
    float lastDayofSeason[4] = { 90, 181, 273, 365 };
    for (int i = 0; i < 4; ++i) {
        if (yday <= lastDayofSeason[i]) {
            GridCell::setSeason(i);
            break;
        }
    }
    
    // read data
    Logger::debug(0, "Reading spatial data from %s", as<std::string>(conf["sasc"]).c_str());
    read_gis(as<std::string>(conf["sasc"]), this);
    
    fishery.initialize(nFisheryBlocks, nGillnetTypes, nFisheryDataSampleSize);
    
    if (conf["fish"] == R_NilValue || conf["effort"] == R_NilValue) {
        Logger::debug(0, "Fishery or effort data not specified, gillnet agents disabled");
    } else {
        Logger::debug(0, "Reading fishery data from %s", as<std::string>(conf["fish"]).c_str());
        stats_size = read_fishery_data(as<std::string>(conf["fish"]), fishery);
        
        Logger::debug(0, "Reading effort data from %s", as<std::string>(conf["effort"]).c_str());
        effort_size = read_fishing_effort(as<std::string>(conf["effort"]), fishery);
    }
    
    postProcessData();
}

void Settings::initData(int x, int y, int nsurv, int nfish, int nblock, int ntrav, int bsize, int npatch, float meanMaxent1, float meanMaxent2, float meanMaxent3, float meanMaxent4) {
    xmx = x;
    ymx = y;
    ncell = xmx * ymx;
    block_size = bsize;
    nSurveyBlocks = nsurv;
    nFisheryBlocks = nfish;
    nBlocks = nblock;
    
    Block::m_xmax = xmx;
    Block::m_ymax = ymx;
    Block::m_ncell = ncell;
    
    Blocks.reserve(nBlocks);
    for (int i = 0; i < nBlocks; ++i) {
        Blocks.emplace_back();
    }
    Grid.reserve(ncell);
    Patches.reserve(npatch);
    TraversableCells.reserve(ntrav);
    
    abundanceRegions.resize(nSurveyBlocks);
    FisheryBlocks.resize(nFisheryBlocks);
    
    GridCell::meanSeasonalMaxent[0] = meanMaxent1;
    GridCell::meanSeasonalMaxent[1] = meanMaxent2;
    GridCell::meanSeasonalMaxent[2] = meanMaxent3;
    GridCell::meanSeasonalMaxent[3] = meanMaxent4;
}

void Settings::postProcessData() {

    // Go over blocks:
    // 1) delete blocks that contain no traversable cells
    // 2) calculate block seasonal values. These are assumed to be known by porpoises.
    //    block seasonal values are the standardized block seasonal averages. Block-specific
    //    seasonal averages are standardized by dividing by the max seasonal average
    //   i.e. val[b,s] = mean_maxent[b,s] / max_maxent[s], where b = block and s = season
    
    for (int i = 0; i < Blocks.size();) {
        if (Blocks[i].m_cellcount == 0) {
            Blocks.erase(Blocks.begin()+i);
        } else {
            for (auto& cell : Blocks[i].m_cells) {
                Grid[cell].Block = i;
            }
            Blocks[i].calcCenter();
            Blocks[i].calcDensity();
            Blocks[i].setId(i);
            ++i;
        }
    }
    
    for (int i = 0; i < abundanceRegions.size();) {
        if (abundanceRegions[i].empty()) { 
            abundanceRegions.erase(abundanceRegions.begin()+i);
        } else {
            for (auto& cell : abundanceRegions[i]._cells) {
                Grid[cell].abundanceBlock = i;
            }
            abundanceRegions[i].setId(i);
            ++i;
        }
    }
    
    // delete fishing areas which contain no traversable cells,
    // or where the fishing effort is zero for all seasons
    for (int i = 0; i < FisheryBlocks.size();) {
        if (FisheryBlocks[i].size() == 0) { 
            FisheryBlocks.erase(FisheryBlocks.begin() + i);
            fishery.removeBlock(i);
        } else {
            ++i;
        }
    }

    Blocks.shrink_to_fit();
    FisheryBlocks.shrink_to_fit();
    abundanceRegions.shrink_to_fit();
    
    nBlocks = Blocks.size();
    nSurveyBlocks = abundanceRegions.size();
    nFisheryBlocks = FisheryBlocks.size();
    
}

int Settings::cellFromXY(int xmax, int ymax, float x, float y) {
    
    if (!isCoordValid(Vector2df(x, y))) {
        return -1;
    }

    int row = floor(ymax-y);
    // points on the edge between rows resolve to the top row,
    // except for the first row, where there is no top row
    if (row == ymax) --row;
    float col = floor(x);
    if (col == xmax) --col;
    
    int cell = row * xmx + col;
    if (cell < 0 || cell >= ncell) cell = -1;
    return cell;
}

int Settings::cellFromPoint(Vector2df pt) {
    return cellFromXY(xmx, ymx, pt.x, pt.y);
}
                                
int Settings::cellFromXY(float x, float y) {
    return cellFromXY(xmx, ymx, x, y);
}
Vector2df Settings::pointFromCell(int cell) {
    assert(cell >= 0 && cell < ncell);
    float row = floor(cell / xmx);
    float col = cell - row * xmx;
    float x = col + 0.5;
    float y = ymx - (row + 0.5);
    return Vector2df(x, y);
}

Vector2df Settings::randomPoint(int a) {
    int cellnum;
    
    if (a == -1 || a >= nSurveyBlocks) {
        cellnum = *select_randomly(TraversableCells.begin(), TraversableCells.end());
    } else {
        cellnum = abundanceRegions[a].randomCell();
    }
    
    Vector2df pos{pointFromCell(cellnum)};
    pos.x += getRandomFloat(-0.49, 0.49);
    pos.y += getRandomFloat(-0.49, 0.49);
    return pos;
}

// replaced by xyFromCell - remove this
Vector2df Settings::GetXYFromCell(int cell) {
    assert(cell >= 0 && cell < ncell);
    float row = floor(cell / xmx);
    float col = cell - row * xmx;
    float x = col + 0.5;
    float y = ymx - (row + 0.5);
    return Vector2df(x, y);
}

bool Settings::isCoordValid(Vector2df coords) {
    return (coords.x >= 0.0f && coords.x <= xmx && coords.y >= 0.0f && coords.y <= ymx);
}

float Settings::subtract_headings(const float origin, const float destination) {
    float left = origin - destination; // left-hand turn
    float right = destination - origin; // right-hand turn
    if (left < 0.0f) left += 360.0f;
    if (right < 0.0f) right += 360.0f;
    return left < right ? left : right;
}

Vector2df Settings::findPathDeepest(Vector2df currentPos, Vector2df mov, float offset, float step, float lookAhead) {
    float bestDepth = 0;
    Vector2df newMov(mov); // use current move by default (no change)
    offset *= PI/180; // degrees to radians
    step *= PI/180;
    
    // start at turning angle 0 (not turning at all) to avoid biasing porpoise turns in one direction
    for (float angle = 0; angle <= offset; angle += step) {
        // loop to check both left and right turns
        for (int i = -1; i < 2; i = i + 2) {
            float theta = angle * i;
            Vector2df candidateMov = mov.rotate(theta);
            Vector2df distantPt = candidateMov * lookAhead;
            
            if (IsPathTraversable(currentPos.x, currentPos.y, currentPos.x + distantPt.x, currentPos.y + distantPt.y)) {
                int cell = cellFromPoint(currentPos + candidateMov);
                if (cell != -1) {
                    float bathy = Grid[cell].Bathymetry;
                    if (bathy < bestDepth) {
                        bestDepth = bathy;
                        newMov = candidateMov;
                    }
                }
            }
            if (theta == 0) break; // avoid checking angle zero twice
        }
    }
    
    return newMov;
}

Vector2df Settings::findPathFarthestFromShore(Vector2df currentPos, Vector2df mov, float offset, float step) {
    float bestDist = 0;
    Vector2df newMov(mov);
    offset *= PI/180;
    step *= PI/180;
    
    // start at turning angle 0 (not turning at all) to avoid biasing porpoise turns in one direction
    for (float angle = 0; angle <= offset; angle += step) {
        // loop to check both directions
        for (int i = -1; i < 2; i = i+ 2) {
            float theta = offset * i;
            Vector2df candidateMov = mov.rotate(theta);
            
            if (IsPathTraversable(currentPos.x, currentPos.y, currentPos.x + candidateMov.x, currentPos.y + candidateMov.y)) {
                int cell = cellFromPoint(candidateMov);
                if (cell != -1) {
                    float dist = Grid[cell].DistanceToCoast;
                    if (dist < bestDist) {
                        bestDist = dist;
                        newMov = candidateMov;
                    }
                }
            }
            if (theta == 0) break; // avoid checking angle zero twice
        }
    }
    
    return newMov;
}

Vector2df Settings::findPathParallellToCoast(Vector2df currentPos, Vector2df mov, float offset, float step, float min, float max) {
    Vector2df newMov(mov); // use current move by default (no change)
    int cell = cellFromPoint(currentPos);
    if (cell == -1) return newMov; // for now.
    
    float currentDist = Grid[cell].DistanceToCoast;
    float bestDist = -1; // impossible value. We'll use this to check if this var is "initalized" or not
    offset *= PI/180; // degrees to radians
    step *= PI/180;
    
    // start at turning angle 0 (not turning at all) to avoid biasing porpoise turns in one direction
    for (float angle = 0; angle <= offset; angle += step) {
        // loop to check both left and right turns
        for (int i = -1; i < 2; i = i + 2) {
            float theta = angle * i;
            Vector2df candidateMov = mov.rotate(theta);
            
            if (IsPathTraversable(currentPos.x, currentPos.y, currentPos.x + candidateMov.x, currentPos.y + candidateMov.y)) {
                int cellno = cellFromPoint(currentPos + candidateMov);
                if (cellno == -1) continue;
                
                float dist = Grid[cellno].DistanceToCoast;
                if (bestDist == -1) {
                    bestDist = dist;
                    newMov = candidateMov;
                    continue;
                }
                
                //if ((currentDist > max && dist > min && dist < bestDist) || (currentDist < min && dist < max && dist > bestDist)) {
                if ((currentDist > max && dist < bestDist) || (currentDist < min && dist > bestDist)) {
                    // optimizing for returning to 1-4 km wide corridor parallel to coast, from position outside corridor
                    bestDist = dist;
                    newMov = candidateMov;
                } else if (currentDist >= min && currentDist <= max && abs(dist - currentDist) < bestDist) {
                    // optimizing for staying at current depth in corridor parallel to coast
                    bestDist = abs(dist - currentDist);
                    newMov = candidateMov;
                }
            }
            if (theta == 0) break; // avoid checking angle zero twice
        }
    }
    
    return newMov;
}

void Settings::adjustMoveToAvoidShallowWater(Porpoise* porp, Vector2df& newPos, float& TurningAngle, float& Distance) {
    
    int currentCell = porp->currentCell;
    int destinationCell = cellFromPoint(newPos);
    
    // if the porpoise does not leave the current cell, we don't have a shallow water problem.
    if (destinationCell == currentCell) return;
    
    // if the porpoise is in a cell from which the distance to the nearest 
    // cell with intraversable (i.e. shallow) water AND the boundaries of
    // the landscape are greater than its desired move distance, then it's
    // not possible for the porpoise to reach land or the map edge with that move,
    // so we can skip any shallow water checks
    /*
    if (currentCell != -1) {
        auto& cell = Grid[currentCell];
        int dist = 1.25 * Distance;
        
        if (cell.DistanceToCoast > dist && cell.DistanceToEdge > dist) {
            return;
        }
    }
    */
    // Otherwise, we need to check the water depth ahead. 
    const Vector2df& currentPos{porp->currentPos};
    const float& heading{porp->Heading};
    
    bool pathTraversable = IsPathTraversable(currentPos.x, currentPos.y, newPos.x, newPos.y);
    
    if (pathTraversable == true) return;
    
    // if we get to this point, that means at least one of the cells in the
    // porpoise's current trajectory would be intraversable/too shallow. 
    // So we'll reiteratively alter the turning angle until we find a traversable path.
    float adjX, adjY; // holds XY coordinates for the end point of adjusted paths

    std::vector<int> tryAngles = { 10, 25, 50, 75, 100, 125, 150 };
    int step = 0;
    int newAngle;
    
    while (pathTraversable == false && step < tryAngles.size()) {
        std::vector<int> angleTypes;
        angleTypes.push_back(getRandomInt(0, 1));
        angleTypes.push_back(1 - angleTypes[0]);
        
        for (int i = 0; i < 2; i++) {
            // pick the next turning angle candidate
            newAngle = tryAngles[step];
            if (angleTypes[i] == 0) newAngle = newAngle * -1;
            // if we've not exhausted our list of candidates < 180 degrees,
            // add a random component to the turning angle
            if (step < tryAngles.size()-1) newAngle += getRandomNormal(0, 2.5);
            // calculate a new trajectory and recheck traversability 
            adjX = currentPos.x + sin((heading + TurningAngle + newAngle) * PI/180) * Distance;
            adjY = currentPos.y + cos((heading + TurningAngle + newAngle) * PI/180) * Distance;
            pathTraversable = IsPathTraversable(currentPos.x, currentPos.y, adjX, adjY);
            if (pathTraversable) {
                newPos.x = adjX;
                newPos.y = adjY;
                TurningAngle += newAngle;
                return;
            }
        }
        step++;
    }
    
    // if a path couldn't be found, have the porpoise go back to the last visited location
    if (porp->track.size() > 1) {
        newPos = porp->lastPos;
        Distance = newPos.distanceFrom(currentPos);
        TurningAngle = 180.0f;
    } else {
        // if for some reason that didn't work out either, move the porp around a bit in its current cell
        newPos = pointFromCell(currentCell);
        newPos.x +=  getRandomFloat(-0.49, 0.49);
        newPos.y += getRandomFloat(-0.49, 0.49);
        Distance = newPos.distanceFrom(currentPos);
    }
}


// IsPathTraversable: iterates over all the individual cells found by 
// running the GetCellsIntersected function (i.e. those cells that
// intersect with a line segment as defined by the point XY and the 
// point newXnewY) and checks the water depth in those cells. If at
// least one cell is too shallow, the function return false.
// Cells that fall outside the grid are considered traversable if
// offGridCellsTraversable is true.
bool Settings::IsPathTraversable(float x, float y, float newX, float newY) {
    
    if (!isCoordValid(Vector2df(newX, newY))) return false;

    std::vector<int> cellsInPath = GetCellsIntersected(x, y, newX, newY);
    if (cellsInPath.size() == 0) return false;
    
    for (const int& cell : cellsInPath) {
        if (cell == -1 || Grid[cell].Bathymetry > MinimumWaterDepth) {
            return false;
        }
    }
    return true;
}
// Fast Voxel Traversal Algorithm
// https://www.researchgate.net/publication/2611491_A_Fast_Voxel_Traversal_Algorithm_for_Ray_Tracing
// https://github.com/cgyurgyik/fast-voxel-traversal-algorithm/blob/master/overview/FastVoxelTraversalOverview.md
std::vector<int> Settings::GetCellsIntersected(float X, float Y, float newX, float newY) {
    auto frac_pos = [](float x) { return 1 - x + floorf(x); };
    auto frac_neg = [](float x) { return x - floorf(x); };
    
    float xstep = newX - X;
    float ystep = newY - Y;
    float tDeltaX = xmx;
    float tDeltaY = ymx;
    int dx = (xstep > 0 ? 1 : (xstep < 0 ? -1 : 0));
    int dy = (ystep > 0 ? 1 : (ystep < 0 ? -1 : 0));
    if (dx != 0) tDeltaX = abs(fmin(dx / xstep, xmx));
    if (dy != 0) tDeltaY = abs(fmin(dy / ystep, ymx));
    float tMaxX, tMaxY;
    
    if (dx > 0) {
        tMaxX = frac_pos(X) * tDeltaX;
    } else {
        tMaxX = frac_neg(X) * tDeltaX;
    }
    
    if (dy > 0) {
        tMaxY = frac_pos(Y) * tDeltaY;
    } else {
        tMaxY = frac_neg(Y) * tDeltaY;
    }
    
    int x = (int) X;
    int y = (int) Y;
    std::vector<int> cellNumbers;
    int j = 0;
    while (true) {
        if (tMaxX < tMaxY) {
            x += dx;
            tMaxX += tDeltaX;
        } else {
            y += dy;
            tMaxY += tDeltaY;
        }
        cellNumbers.push_back(cellFromXY(x, y));
        j++;
        if ((tMaxX > 1 && tMaxY > 1) || j == 200) break;
    }
    return cellNumbers;
}
