#include <algorithm>
#include <string>
#include <random>
#include <iostream>
#include <cmath>
#include <Rcpp.h>
#include "execHalfhourTasks.h"
#include "Settings.hpp"
#include "Porpoise.hpp"
#include "GridCell.hpp"
#include "misc.hpp"
#include "Vector2d.hpp"
#include "Logger.h"

void execMonthlyTasks(Settings& sim) {

    /*
    // recalculate block values based on abundances
    for (auto& block : sim.Blocks) {
        block.m_N = 0;
    }

    // todo: consider doing this directly as porps move
    for (const auto& porp : Porpoise::Porpoises) {
        if (porp->currentCell != -1) {
            int blockno = sim.Grid[porp->currentCell].Block;
            if (blockno > 0 && blockno < sim.nBlocks) {
                sim.Blocks[blockno].m_N++;
            }
        }
    }
 
    for (auto& block : sim.Blocks) {
        block.calcValue();
    }
    */
    // add up total food in all patches
    float food = 0;
    
    for (const int& i : sim.Patches) {
        food += sim.Grid[i].CurrentUtility;
    }
    
    // add up total energy across all porpoises
    // AND count the number of porpoises of each age class
    // in each abundance area.
    float energy = 0;
    std::vector<int> N(3, 0); // juvenile, adult, old
    std::vector<int> blockN(sim.nSurveyBlocks, 0);
    
    if (Porpoise::Porpoises.size() > 0) {
        for (const auto& porp : Porpoise::Porpoises) {
            
            int cell = porp->currentCell;
            if (cell >= 0 && cell < sim.ncell) {
                int block = sim.Grid[cell].abundanceBlock;
                if (block >= 0 && block < sim.nSurveyBlocks) {
                    blockN[block]++;
                }
            }
            int ageClass = 0;
            if (porp->Age >= 4) {
                if (porp->Age >= 10) {
                    ageClass = 2;
                } else {
                    ageClass = 1;
                }
            }
            N[ageClass]++;
            energy += porp->EnergyLevel;
        }
        if (energy > 0 && Porpoise::Porpoises.size() > 0) {
            energy /= Porpoise::Porpoises.size();
        } else {
            energy = 0;
        }
    
    }
    
    // age structure
    for (int ageClass = 0; ageClass < 3; ++ageClass) {
        sim.logger->log(sim.time->step(), ageClass, N[ageClass], 1);
    }
    
    // abundance in survey blocks
    for (int block = 0; block < sim.nSurveyBlocks; ++block) {
        abundanceRegion& aR = sim.abundanceRegions[block];
        sim.logger->log(sim.time->step(), block, blockN[block], aR._births, aR._deaths, aR._bycatch);
        aR.reset_stats();
    }
    
    sim.logger->log(sim.time->step(), Porpoise::Porpoises.size(), food, energy);
}
