#include <algorithm>
#include <string>
#include <random>
#include <iostream>

#include <Rcpp.h>
#include "Settings.hpp"
#include "Porpoise.hpp"
#include "GridCell.hpp"
#include "misc.hpp"
#include "Vector2d.hpp"
#include "Logger.h"
#include "execHalfhourTasks.h"

void execDailyTasks(Settings& sim) {

    // regrow food
    for (const int& patch : sim.Patches) {
        sim.Grid[patch].Regenerate();
    }
    const int& season = sim.time->quarter() - 1;
    const int& yday = sim.time->yday() - 1;
    
    // create new gillnet agents
    int netcount = 0; // number of hauls this day

    for (int block = 0; block < sim.fishery.nBlock(); ++block) {
        for (int type = 0; type < sim.fishery.nType(); ++type) {
            int newSets = sim.fishery.sampleNumberOfSets(block, yday, season, type);
            if (newSets > 0) {
                for (int netcounter = 0; netcounter < newSets; ++netcounter) {
                    FishingEffort effort = sim.fishery.sampleEffort(block, yday, season, type);
                    sim.Gillnets.emplace_back(block, type, effort.soaktime, effort.length, effort.pinger);
                }
                netcount += newSets;
            }
        }
    }

    sim.logger->gillnet_set(netcount);

    // iterate using indices, because we'll add new porpoises as calves are weaned
    // since we're adding to the back of the vector, we won't iterate over the newborn calves
    if (sim.time->step() > 0) {
        int N = Porpoise::Porpoises.size();
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            auto &porp = Porpoise::Porpoises[i];
            
            // increase age by 1 day (1/365)
            porp->Age += 0.002739726f;
            
            // dispersal
            if (sim.Blocks.size() > 1) {
                porp->considerDispersing();
            }
    
            // mating
            if (yday == porp->matingDay) {
                porp->Mate();
            }
            
            // giving birth
            if (porp->isPregnant && yday == porp->calfBirthday) {
                porp->giveBirth();
            }
            
            // wean calf (weaningDay is -1 unless porp has a calf, so we can check it like this)
            if (yday == porp->weaningDay) {
                // calf is only added to population if it is female, assuming a sex ratio of 1:1
                if (getRandomFloat(0, 1) < 0.5) {
                    #pragma omp critical
                    porp->weanCalf();
                    
                    // logging. move this code block elsewhere. perhaps to logger class?
                    int cell = porp->currentCell;
                    
                    if (cell >= 0 && cell < sim.ncell) {
                        
                        int abundanceRegion = sim.Grid[cell].abundanceBlock;
                        
                        if (abundanceRegion >= 0 && abundanceRegion < sim.abundanceRegions.size()) {
                            #pragma omp atomic
                            sim.abundanceRegions[abundanceRegion]._births++;
                            
                        }
                    }
                    
                } else {
                    porp->abandonCalf();
                }
            }
            
        }
    }

}
