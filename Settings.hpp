#ifndef __SETTINGS__
#define __SETTINGS__
#include <Rcpp.h>

#include "misc.hpp"
#include "Vector2d.hpp"
#include "Logger.h"
#include "Timer.hpp"
#include "Block.hpp"
#include "abundanceRegion.hpp"
#include "Fishery.hpp"
#include "FishingEffort.hpp"

// forward declarations
class GridCell;
class Gillnet;
class Block;
class abundanceRegion;
class Logger;
class Porpoise;
extern int DEBUG_LEVEL;

class Settings {
public:
    // simulation data
    Fishery fishery;
    std::list<Gillnet> Gillnets; // for fast random erase
    std::vector<int> TraversableCells;
    std::vector<std::vector<int>> FisheryBlocks;
    std::vector<abundanceRegion> abundanceRegions;
    std::vector<Block> Blocks;
    std::vector<GridCell> Grid;
    std::vector<int> Patches;
    Logger *logger;
    Timer *time;
    int xmn = 0;
    int ymn = 0;
    int xmx, ymx, ncell, block_size;
    int nSurveyBlocks, nFisheryBlocks, nBlocks;
    int nFisheryDataSampleSize, nGillnetTypes;
    bool offGridCellsTraversable;
    float MinimumWaterDepth;
    float FoodGrowthRate; 
    float maxU;
    int QuarterStartingDay[4] = {1, 91, 182, 274};
    std::vector<int> validSurveyBlocks;
    float gillnet_interaction_maxp = 0.0531923;
    int effort_size = 0;
    int stats_size = 0;
    float pinger_effect {0.5};

    // functions
    Settings(Rcpp::List& conf);
    void initData(int x, int y, int nsurv, int nfish, int nblock, int ntrav, int bsize, int npatch,  float meanMaxent1, float meanMaxent2, float meanMaxent3, float meanMaxent4);
    void postProcessData();
    std::vector<int> GetCellsIntersected(float X, float Y, float newX, float newY);
    bool IsPathTraversable(float X, float Y, float newX, float newY);
    int cellFromPoint(Vector2df pt);
    int cellFromXY(float X, float Y);
    int cellFromXY(int xmax, int ymax, float x, float y);
    Vector2df pointFromCell(int cell);
    Vector2df GetXYFromCell(int cellnum);
    bool isCoordValid(Vector2df coords);
    Vector2df randomPoint(int abundanceRegion = -1);
    Vector2df findPathDeepest(Vector2df currentPos, Vector2df mov, float offset, float step, float lookAhead);
    Vector2df findPathFarthestFromShore(Vector2df currentPos, Vector2df mov, float offset, float step);
    Vector2df findPathParallellToCoast(Vector2df currentPos, Vector2df mov, float offset, float step, float min, float max);
    void adjustMoveToAvoidShallowWater(Porpoise* porp, Vector2df& newPos, float& TurningAngle, float& Distance);
    static float subtract_headings(const float origin, const float destination);
    /*void calcBlockAverageFood();
    int blockFromXY(float x, float y);
    int rowFromCell(const float& xmx, const int& cell);
    int colFromCell(const float& xmx, const float& ymx, const int& cell);
    std::vector<int> adjacentCells(float xmx, float ymx, int block);
    std::vector<int> adjacentCells(int cell);
    std::vector<int> adjacentBlocks(int block);
     */

};

#endif // __SETTINGS__