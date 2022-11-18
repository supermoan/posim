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

void execYearlyTasks(Settings& sim) {
    // set mating day
    for (auto &porp : Porpoise::Porpoises) {
        porp->setMatingDay();
    }
}
