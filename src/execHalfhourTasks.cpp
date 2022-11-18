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
#include "Block.hpp"

extern pcg32 rng; // import from misc.cpp

void execHalfhourTasks(Settings& sim) {
    
    // let porpoises do their thing
    // but randomize the order in which they act
    int N = Porpoise::Porpoises.size();
    std::vector<int> casualties; // holds indices of porpoises that died in this step
    std::vector<int> indices(N);
    std::iota(indices.begin(), indices.end(), 0); // indices from 0 to total number of porps
    std::shuffle(indices.begin(), indices.end(), rng); // randomize indices
    
    #pragma omp parallel for
    for (int j = 0; j < N; ++j) {
        int i = indices[j];
        auto &porp = Porpoise::Porpoises[i];
        
        if (porp->currentCell == -1) {
            Rprintf("porp %d is off-grid!\n", i);
            continue;
        }
        bool entangled = false;
        porp->dispersed = false;
        
        // porpoise dies from old age
        if (porp->Age >= Porpoise::max_age) {
            #pragma omp critical 
            casualties.push_back(i);
            continue;
        }

        // move (correlated random walk + memory)
        porp->intrinsicMove();
        
        // gillnet interaction: if porp is entangled in a gillnet, report and skip to next iteration
        # pragma omp critical
        entangled = porp->Entangled();
        
        if (entangled) {
            int cell = porp->currentCell;
            if (cell != -1) {
                int abundanceRegion = sim.Grid[cell].abundanceBlock;
                if (abundanceRegion != -1) {
                    #pragma omp atomic
                    sim.abundanceRegions[abundanceRegion]._bycatch++;
                }
            }
            #pragma omp critical
            sim.logger->log(sim.time->step(), porp);
            //Logger::debug(0, "Day %d: porp %d got entangled during step %d moving from (%.02f, %.02f) to (%.02f, %.02f)", time.day(), porp->Id, time.step(), porp->X[1], porp->Y[1], porp->X[0], porp->Y[0]);
            #pragma omp critical
            casualties.push_back(i);
            continue;
        }

        // consume food in patch
        #pragma omp critical
        porp->consumeFood();
            
        // dispersal
        
        if (porp->movementMode == Porpoise::directedDispersal || porp->movementMode == Porpoise::returningDispersal) {
            porp->disperseTowardsTarget(); 
        } else if (porp->movementMode == Porpoise::coastalDispersal) {
            porp->disperseAlongCoast();
        }
  
        if (porp->dispersed) {
            porp->dispersalStepCounter++;
        } else {
            porp->dispersalStepCounter = 0;
        }
        
        porp->useEnergy();
        
        if (!porp->checkEnergy()) { // if energy is too low, porp dies.
            #pragma omp critical
            sim.logger->log(sim.time->step(), porp);
            #pragma omp critical 
            casualties.push_back(i);
            int cell = porp->currentCell;
            if (cell != -1) {
                int abundanceRegion = sim.Grid[cell].abundanceBlock;
                if (abundanceRegion != -1 && abundanceRegion < sim.abundanceRegions.size()) {
                    #pragma omp atomic
                    sim.abundanceRegions[abundanceRegion]._deaths++;
                }
            }
            continue;
        }
        #pragma omp critical 
        sim.logger->log(sim.time->step(), porp);
    }
    
    // remove dead porpoises
    if (casualties.size() > 0) {
        std::sort(casualties.begin(), casualties.end());
        int counter = 0;
        for (const int& i : casualties) {
            Porpoise::Porpoises[i + counter] = std::move(Porpoise::Porpoises.back());
            Porpoise::Porpoises.pop_back();
            //Porpoise::Porpoises.erase(Porpoise::Porpoises.begin() + i + counter);
            --counter;
        }
    }
    
    int bycatch = 0;
    int hauled = 0;
    
    // haul (remove) gillnets that have reached their maximum soaktime
    if (sim.Gillnets.size() > 0) {
        auto it = sim.Gillnets.begin();
        while (it != sim.Gillnets.end()) {
            it->soak30m();
            if (it->haulable()) {
                bycatch += it->bycatch(); // tally up number of bycaugh porpoises in net
                ++hauled;
                sim.Gillnets.erase(it++);
            } else {
                ++it;
            }
        }
    }
    sim.logger->bycatch(bycatch);

}
