// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include <cassert>
#include <list>

// for unix-alike machines only
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
#include <unistd.h>
#include <Rinterface.h>
#endif

// load global configuration first, as other classes depend on it.
#include "pcg_random.hpp"
#include "misc.hpp"
#include "Timer.hpp"
#include "Logger.h"
#include "Settings.hpp"
#include "Vector2d.hpp"
#include "GridCell.hpp"
#include "Porpoise.hpp"
#include "Gillnet.h"
#include "IO.hpp"

#include "execHalfhourTasks.h"
#include "execDailyTasks.h"
#include "execMonthlyTasks.h"
#include "execQuarterlyTasks.h"
#include "execYearlyTasks.h"
#include "omp.h"

typedef std::pair<Vector2df, Vector2df> _linestring;

int DEBUG_LEVEL = 0;

// [[Rcpp::export]]
Rcpp::RObject do_sim(Rcpp::List& RSim) {
    // for now, assume settings are ok from the R side - no checks on parameters
    // todo: add handler to validate configuration before doing stuff
    Rcpp::List conf = RSim["conf"];
    DEBUG_LEVEL = Rcpp::as<int>(conf["debug"]);
    Rcpp::IntegerVector N = conf["N"];
    Porpoise::resetId();
    int nThread = Rcpp::as<int>(conf["nThread"]);
    int nMaxThread = omp_get_max_threads();
    
    if (nThread > nMaxThread) {
        nThread = nMaxThread;
        Rcpp::warning("nThread set to %d (maximum allowable)", nThread);
    }
    
    omp_set_num_threads(nThread);
    std::vector<int> follow = Rcpp::as<std::vector<int>>(conf["follow"]);
    int steps = Rcpp::as<int>(conf["steps"]);
    int start = Rcpp::as<int>(conf["start"]);
    
    Logger logger(steps, follow, Rcpp::as<int>(conf["maxAge"]));
    Logger::debug(0, "This is PorpSIM v0.1, using up to %d threads", nThread);
    Logger::debug(0, "Logger set up with debug = %d", DEBUG_LEVEL);
    
    Settings sim{conf};
    
    Timer time(start);
    sim.logger = &logger;
    sim.time = &time;
    
    int ngillnetcells = 0;
    for (const auto& b : sim.FisheryBlocks) ngillnetcells += b.size();
    
    Logger::debug(0, "Landscape is %d x %d cells (%d in total, %d traversable, %d usuable for gillnets, %d food patches)", sim.xmx, sim.ymx, sim.ncell, sim.TraversableCells.size(), ngillnetcells, sim.Patches.size());
    Logger::debug(0, "There are %d traversable blocks, %d abundance region(s) and %d fishery region(s)", sim.Blocks.size(), sim.abundanceRegions.size(), sim.FisheryBlocks.size());
    
    float curfood = 0;
    float totalfood = 0;
    for (auto& i : sim.Patches) {
        curfood += sim.Grid[i].CurrentUtility;
        totalfood += sim.Grid[i].currentMax;
    }
    
    Logger::debug(0, "maxU = %.02f, starting systemic food = %.02f/%.02f", sim.maxU, curfood, totalfood);
    
    Porpoise::Porpoises.reserve(sum(N)*10);
    // Create initial population of porpoise agents
    int Unstructured = N.size() == 1;
    
    for (int j = 0; j < N.size(); ++j) {
        int blockN = N[j];
        if (j - Unstructured == -1) {
            Logger::debug(0, "Created %d porpoises, distributed randomly in landscape", blockN);
        } else {
            Logger::debug(0, "Created %d porpoises in survey block %d", blockN, j-Unstructured);
        }
        
        for (int i = 0; i < blockN; ++i) {
            Porpoise::Porpoises.push_back(std::unique_ptr<Porpoise>(new Porpoise(j-Unstructured)));
            logger.log(0, Porpoise::Porpoises.back());
        }
    }

    Logger::debug(0, "Simulation started on yday %d", time.yday());

    // delegate work to procedures according to increases in step/day/month/year counters
    while (time.step() <= steps) {

        if (time.isNewDay()) {
            
            #if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
                float pct = time.step(); // cast to float
                pct = (pct / steps) * 100;
                REprintf("\rProcessing step %d of %d (%.01f%%)\r", time.step(), steps, pct);
                R_FlushConsole();
            #endif
        
            if (time.isNewMonth()) {
                
                if (time.isNewYear()) execYearlyTasks(sim);
                if (time.isNewQuarter()) execQuarterlyTasks(sim);
                
                execMonthlyTasks(sim);
                Rcpp::checkUserInterrupt(); // allow user to interrupt simulation from R
            }
            execDailyTasks(sim);
        }
        
        execHalfhourTasks(sim);

        if (Porpoise::Porpoises.size() == 0) {
            Logger::debug(0, "Population is extinct!");
            break;
        }
        time.next(); // goto next iteration
    }
    
    execMonthlyTasks(sim);
    
    if (time.step() > steps) {
        Logger::debug(0, "Simulation completed after %d steps", time.step());
    } else {
        Logger::debug(0, "Stopped prematurely on simulation day %d (yday %d, step %d), because population is extinct!", time.day(), time.yday(), time.step());
    }
    

    
    sim.Gillnets.clear(); // important that this is called before sim goes out of scope, since gillnet destructors free memory on grid cells
    RSim["result"] = clone(logger.toList());
    Porpoise::Porpoises.clear();
    //Porpoise::sim = nullptr;
    //Block::sim = nullptr;
    return RSim;
}
