#include "GridCell.hpp"

float GridCell::foodGrowthRate = 0.2; // overwritten in sim constructor
int GridCell::season = 0; // overwritten in sim constructor
float GridCell::meanSeasonalMaxent[4]; // set in sim::initData()

GridCell::GridCell(int id, int block, int abundanceRegion, int fisheryBlock, float bathy, float distToCoast, float distToEdge, float food, float maxent1, float maxent2, float maxent3, float maxent4) : 
    Id(id), Block(block), abundanceBlock(abundanceRegion), fisheryBlock(fisheryBlock), Bathymetry(bathy), DistanceToCoast(distToCoast), DistanceToEdge(distToEdge), MaximumUtility(food) {
    
    maxentLevel[0] = maxent1;
    maxentLevel[1] = maxent2;
    maxentLevel[2] = maxent3;
    maxentLevel[3] = maxent4;
    
    if (MaximumUtility > 0) {
        updateMax();
        // initialize with extra food to compensate for porpoises starting out with a blank memory
        CurrentUtility = MaximumUtility * 1 / meanSeasonalMaxent[season];
    }
}

void GridCell::Regenerate() {
    if (CurrentUtility < currentMax) {
        //CurrentUtility += foodGrowthRate * CurrentUtility * (1 - CurrentUtility / (currentMax / meanSeasonalMaxent[season]));
        
        float food = CurrentUtility + foodGrowthRate * CurrentUtility * (1 - CurrentUtility / currentMax);

        if (fabs(food - CurrentUtility) > 0.001) {
            for (int i = 0; i < 47; ++i) {
                food += foodGrowthRate * food * (1 - food / currentMax);
            }
        }

        CurrentUtility = food;
    }
}

void GridCell::setSeason(int newSeason) {
    season = newSeason;
}

void GridCell::updateMax() {
    currentMax = MaximumUtility * maxentLevel[season] / meanSeasonalMaxent[season];
}

float GridCell::maxent(int season) const {
    return maxentLevel[season];
}
