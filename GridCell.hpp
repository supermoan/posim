#ifndef __GRIDCELL__
#define __GRIDCELL__

#include <list>
#include <cmath>

class Settings;
class IO;
class Gillnet;
class Block;

class GridCell {
private:
    static int season;
    static float foodGrowthRate;
    static float meanSeasonalMaxent[4];
public:
    int Id;
    int Block { -1 };
    int abundanceBlock { -1 };
    int fisheryBlock { -1 };
    float Bathymetry { 0.0f }; // Average water depth in cell
    float DistanceToCoast { 0.0f }; // Euclidean distance to coast in units of number of cells
    float DistanceToEdge { 0.0f }; // Distance to nearest edge of landscape
    float CurrentUtility { 0.0f }; // current food level in patch
    float MaximumUtility = { 0.0f }; // max food level in patch, not adjusted for seasonal maxent
    float maxentLevel[4] { 1.0f, 1.0f, 1.0f, 1.0f }; // maxent level per quarter
    float currentMax { 0.0f }; // max adjusted by seasonal maxent
    std::list<Gillnet*> gillnets; // list of gillnets in cell at any given time step
    GridCell(int id, int block, int abundanceBlock, int fisheryBlock, float bathy, float distToCoast, float distToEdge, float food, float maxent1, float maxent2, float maxent3, float maxent4);
    void Regenerate(); // regrows food in cell logistically (if there was food to begin with)
    static void setSeason(int newSeason); // sets the static season variable to the value specified
    void updateMax(); // calculates a new max food level based on current season
    float maxent(int season) const; // returns the current maxent level
    friend class Settings;
    friend class IO;
    friend class Block;
};

#endif // __GRIDCELL__