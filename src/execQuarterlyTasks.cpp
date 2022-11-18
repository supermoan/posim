#include <algorithm>
#include <string>
#include <random>
#include <iostream>

#include <Rcpp.h>
#include "execHalfhourTasks.h"
#include "Settings.hpp"
#include "Porpoise.hpp"
#include "GridCell.hpp"
#include "misc.hpp"
#include "Vector2d.hpp"
#include "Logger.h"

class GridCell;


void execQuarterlyTasks(Settings& sim) {

    GridCell::setSeason(sim.time->quarter() -1);
    float totfood = 0;
    
    for (const int& patch : sim.Patches) {
        sim.Grid[patch].updateMax();
        totfood += sim.Grid[patch].currentMax;
    }
    
    //for (auto& block : sim.Blocks) {
    //    block.calcDensity();
    //}
    
        /*
    for (auto& block : sim->Blocks) {
        int blockN = std::count_if(Porpoise::Porpoises.begin(), Porpoise::Porpoises.end(), 
                                   [sim, block](const auto& porp) { 
                                       return sim->Grid[porp->currentCell]->Block == block.id();
                                   });
        if (blockN == 0) blockN = 1;
        block.setN(blockN); 
        block.calcValue();
    }
         */

}
