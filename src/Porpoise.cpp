#include <Rcpp.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "Settings.hpp"
#include "Porpoise.hpp"
#include "GridCell.hpp"
#include "misc.hpp"
#include "Vector2d.hpp"
#include "Logger.h"
#include "Block.hpp"
#include "Position.hpp"

typedef std::pair<Vector2df, Vector2df> _linestring;

// static porpoise variables
std::vector<std::unique_ptr<Porpoise>> Porpoise::Porpoises;
Settings* Porpoise::sim;
int Porpoise::nextId;
std::shared_ptr<Settings> Porpoise::Config;
float Porpoise::AgeOfMaturity;
float Porpoise::monthlyEnergyMultiplier[12];
float Porpoise::withCalfEnergyMultiplier;
float Porpoise::distEnergyMultiplier;
float Porpoise::stepEnergyMultiplier;
float Porpoise::m_mort_prob;
float Porpoise::x_surv_prob;
int Porpoise::dispersalInertia;
float Porpoise::meanDispersalDistance;
float Porpoise::minDispersalDistance;
float Porpoise::maxDispersalDistance;
float Porpoise::minDispersalDistanceSquared;
float Porpoise::maxDispersalDistanceSquared;
float Porpoise::minDispersalDepth;
float Porpoise::minDispersalDistanceToLand;
float Porpoise::pregnancy_prob;
float Porpoise::max_age;
int Porpoise::maxMemory;
float Porpoise::inertia_const;
float Porpoise::corrLogmov;
float Porpoise::corrAngle;
float Porpoise::m;
float Porpoise::maxLogmov;
std::vector<int> Porpoise::age_dist;
std::vector<float> Porpoise::ref_mem_strength;
std::vector<float> Porpoise::work_mem_strength;

// constructor for all porpoises in starting population
Porpoise::Porpoise(int SurveyBlock) : Id(++nextId) {

    track.reserve(maxMemory);

    // assign a random position to the porpoise (within the specified abundance block, if applicable)
    executeMove(sim->randomPoint(SurveyBlock), Heading, normalMove);

    // set birthday, assuming porpoises are mostly born around june 9
    int birthday = round(getRandomNormal(160, 20));
    
    // pick a random age class
    Age = getRandomDiscrete(&age_dist);

    // increase age by some fraction of a year, corresponding to the time elapsed since the porp's birthday
    if (sim->time->yday() >= birthday) {
        Age += (sim->time->yday() - birthday) * 0.002739726;
    } else {
        Age += (365 - birthday + sim->time->yday()) * 0.002739726;
    }
    
    // if porp was mature last mating season, she may be pregnant with a calf
    if ((Age-1) >= AgeOfMaturity) {
        if (getRandomFloat(0, 1) < pregnancy_prob) {
            float birth = getRandomNormal(160, 20);
            if (birth > sim->time->yday()) {
                isPregnant = true;
                calfBirthday = birth;
            }
        }
        
        // if porp was mature two mating seasons ago, she may be nursing a calf born last summer
        if ((Age - 2) >= AgeOfMaturity && getRandomFloat(0, 1) < pregnancy_prob*0.5) {
            float wean = getRandomNormal(100, 20);
            if (wean > sim->time->yday()) {
                withCalf = true;
                weaningDay = wean;
            }
        }
    }
    
    // set current energy usage according to time of year and calf status
    setEnergyUse();
    setMatingDay();
}

// Calf constructor. Calves start out in the same location as their mother. 
Porpoise::Porpoise(const Porpoise& mother) : Id(++nextId) {
    track.reserve(maxMemory);
    Age = 0.6777778f; // calves are weaned after 8 months
    executeMove(mother.currentPos, mother.Heading, normalMove);
    setEnergyUse();
    setMatingDay();
}

/***********************************************************************************************
 * 
 * 
 *                                   MOVEMENT
 * 
 * 
 ***********************************************************************************************/
// convenience function to always keep heading between 0 and 360
void Porpoise::setHeading(float newHeading) {
    Heading = newHeading;
    
    if (Heading >= 360.0f) {
        Heading -= 360.0f;
    } else if (Heading < 0.0f) {
        Heading += 360.0f;
    }
}
void Porpoise::executeMove(Vector2df newPos, float newHeading, PorpoiseMovementMode mode) {
    
    int newCell = sim->cellFromPoint(newPos);
    
    // catch rare bugs here: don't allow moves to illegal cells
    if (newCell == -1 || sim->Grid[newCell].Bathymetry > sim->MinimumWaterDepth) {
        return; 
    }

    currentCell = newCell;
    lastPos = currentPos;
    currentPos = newPos;
    
    if (mode == normalMove) {
        track.emplace(track.begin(), newPos);
        if (track.size() > maxMemory) track.resize(maxMemory);
    }

    // update heading
    setHeading(newHeading);
    //Logger::debug(0, " set heading to %.02f (moved from %.02f,%.02f to %.02f,%.02f", Heading, lastPos.x, lastPos.y, currentPos.x, currentPos.y);
}

void Porpoise::intrinsicMove() {
    
    float expectedEnergy = calcExpectedEnergy();
    float CRW_contrib = inertia_const + pres_mov * expectedEnergy; // emphasize CRW if food is plentiful
    Vector2df foodAttraction = calcFoodAttractionVector(); // food attraction vector
    Vector2df CRW = calcCRW() * CRW_contrib; // CRW component, scaled by expected foraging success
    Vector2df mov = CRW + foodAttraction; // calculate the sum of CRW and food attraction
    mov = mov.normalize() * (0.25f * pres_mov); // calculate unit vector, and multiply by pres_mov (and scale from 100m to 400m grid)
    
    // calculate new absolute position
    Vector2df newPos = currentPos + mov;
    // calculate turning angle of the new move, relative to current heading (subtract from 90 to get from x-axis to y-axis)
    float newHeading = 90 - atan2(mov.y, mov.x)*180/PI;
    pres_angle = newHeading - Heading;
    if (pres_angle < -180.0f) pres_angle += 360;
    if (pres_angle > 180.f) pres_angle -= 360;
    
    //Logger::debug(0, "mov: %.02f,%.02f, current heading: %.02f, new heading: %.02f, turning angle: %.02f, adjusted turn: %.02f",
     //            mov.x, mov.y, Heading, newHeading, old_angle, pres_angle);
    
    // check if the new course would intersect with any cells that do not
    // have sufficient water depth, and make adjustments as needed
    sim->adjustMoveToAvoidShallowWater(this, newPos, pres_angle, pres_mov);
    this->prev_mov = pres_mov;
    this->prev_logmov = log10(pres_mov);
    this->prev_angle = pres_angle;
    //Logger::debug(0, " final turn = %.02f, final heading = %.02f", pres_angle, Heading + pres_angle);
    
    // execute the move (heading recalculated because it may have changed)
    executeMove(newPos, Heading + pres_angle, normalMove);
}

// CRW: Correlated Random Walk. Calculates turning angle and move length with no
// regard for current energy state or memory, which we'll add in later.
Vector2df Porpoise::calcCRW() {
    
    // R1 = N(0.42, 0.48)       Log10 distance moved per time step (mean +/- 1 SD).
    // R2 = N(0, 38)            Turning angles between steps (mean +/- 1 SD).
    // R3 = N(96, 28)           Parameter controlling relationship between turning angles and step length
    
    // pres_mov is measured in 100m steps
    
    // calculate turning angle. Turning angle should be negatively correlated
    // with the previous turning angle and, when prev_mov < m, with distance moved.
    pres_angle = 999;
    float tmp_angle = 0;
    int j = 0;
    
    if (prev_angle < 0) {
        tmp_angle += 24;
    } else {
        tmp_angle -= 24;
    }
    while (abs(pres_angle) > 180 && j < 200) {
        pres_angle = tmp_angle * -1 * corrAngle + getRandomNormal(0, 38);
        j++;
        if (j == 200) {
            pres_angle = pres_angle * 90 / abs(pres_angle);
        }
    }
    
    int sign = 1;
    if (pres_angle < 0) sign = -1;
    pres_angle = abs(pres_angle);
    // make angle decrease linearly with mov_dist
    bool go_on = true;
    float rnd;
    j = 0;
    while (go_on == true && j < 200) {
        rnd = getRandomNormal(96, 28);
        if (prev_mov <= 5.5) pres_angle += rnd - (rnd*prev_mov/5.5);
        if (pres_angle < 180) go_on = false;
        j++;
        if (j == 200) {
            pres_angle = getRandomInt(0, 20) + 90;
            go_on = false;
        }
    }
    
    pres_angle = pres_angle * sign;
    
    // calculate move distance
    pres_logmov = maxLogmov + 1;
    while (pres_logmov > maxLogmov) {
        pres_logmov = corrLogmov * prev_logmov + getRandomNormal(0.42, 0.48);
    }
    pres_mov = pow(10, pres_logmov);
    
    // slow down if turning sharply
    if (pres_mov > 10 && abs(pres_angle) > 90) {
        pres_mov /= 5; 
    } else if (pres_mov > 7 && abs(pres_angle) > 50) {
        pres_mov /= 2;
    }
    
    float heading = Heading + pres_angle;
    
    if (heading < 0) {
        heading += 360;
    } else if (heading > 360) {
        heading -= 360;
    }
    
    // return unit vector
    float x = sin(heading * PI/180);
    float y = cos(heading * PI/180);
    return Vector2df(x, y);
    
}
Vector2df Porpoise::calcFoodAttractionVector() {
    
    Vector2df attractionVector;
    const int n = std::min(track.size(), ref_mem_strength.size());
    
    // calculate distances between current position and positions in recent past
    for (int i = 1; i < n; ++i) { // start loop at 1 => no attraction to current pos
        if (track[i].food == 0) continue; // if porpoise didn't find any food at this location, skip it
        
        Vector2df av = track[1].pos - currentPos; // vector pointing toward known food location
        float attraction; // attraction strength
        float length = av.length(); // distance to that location
        
        if (length < 0.001) {
            length = 0.001;
            attraction = 9999; // large attraction for close patches
        } else {
            attraction = track[i].food * ref_mem_strength[i] * (1/length);
        }
        
        // calculate unit vector
        av *= (1/length);
    
        // multiply unit vector pointing to previous location by attraction score, 
        // and add the result to the total attraction vector
        attractionVector += av * attraction;
    }
    
    return attractionVector;
}


/***********************************************************************************************
 * 
 * 
 *                                   FOOD AND ENERGY
 * 
 * 
 ***********************************************************************************************/


void Porpoise::consumeFood() {

    float& food = sim->Grid[currentCell].CurrentUtility; // current food level in patch
    track.front().food = food; // porp remembers how much food it found here

    // only eat food if 1) there is food and 2) porp is not already at full energy
    if (food > 0 && EnergyLevel < 20) { 
        float propFoodToEat = (20-EnergyLevel) / 10;
        if (propFoodToEat > 0.99) propFoodToEat = 0.99; // porp never consumes ALL food in a patch
        float foodEaten = propFoodToEat * food;
        EnergyLevel += foodEaten;
        // make sure porpoise doesn't get too energized...
        if (EnergyLevel > 20) {
            foodEaten -= EnergyLevel - 20;
            EnergyLevel = 20.0f;
        }
        food -= foodEaten;
        if (food < 0.01) food = 0.01; // patches don't deplete completely
    }
}

float Porpoise::calcExpectedEnergy() {
    const int n = std::min(track.size(), work_mem_strength.size());
    float expectedEnergy{0};
    for (int i = 0; i < n; i++) {
        expectedEnergy += work_mem_strength[i] * track[i].food;
    }
    return expectedEnergy;
}

void Porpoise::setEnergyUse() {
    E_use = monthlyEnergyMultiplier[sim->time->month()-1];
    if (withCalf) E_use *= withCalfEnergyMultiplier;
}

void Porpoise::useEnergy() {
    // calculate expended energy this step and subtract that amount from the current energy level
    // total distance moved is the sum of prev_mov (food search) and dispersal
    float dist = prev_mov * 2.5; // seems to me this should be * 0.4, not * 2.5
    if (dispersed) dist += meanDispersalDistance; // 
    
    float energyUsed = 0.001 * E_use * (stepEnergyMultiplier + dist * distEnergyMultiplier);
    EnergyLevel -= energyUsed;
    cumulativeEnergy += EnergyLevel;
}

bool Porpoise::checkEnergy() {
    bool survived = false;
    
    if (EnergyLevel > 0) { // porps with zero or negative energy die automatically (no need to check survival)
        float yearly_survival = 1 - (m_mort_prob * exp(-EnergyLevel * x_surv_prob));
        float step_survival = exp(log(yearly_survival) / 17520); // 17520 steps in a year
        
        if (getRandomFloat(0, 1) < step_survival) { // porp survives at current energy
            survived = true;
        } else if (withCalf) { // porp survives by sacrificing calf
            abandonCalf();
            survived = true;
        }
    }
    
    return survived;
}

void Porpoise::calcDailyEnergy() {
    if (track.size() >= 48) {
        DailyEnergy.insert(DailyEnergy.begin(), cumulativeEnergy / 48.0f);
        DailyEnergy.resize(10);
        cumulativeEnergy = 0.0f;
    }
}



/***********************************************************************************************
 * 
 * 
 *                                   DISPERSAL
 * 
 * 
 ***********************************************************************************************/






void Porpoise::considerDispersing() {
    // calculate average energy level for the last 24 hours
    calcDailyEnergy();
    
    dailyPositions.insert(dailyPositions.begin(), currentPos);
    dailyPositions.resize(10);
    
    if (movementMode == normalMove) { // not dispersing
        
        // iterating from 0 to 8, not 0 to 9 (due to the comparison with i+1 inside the loop)
        //for (int i = 0; i < 9; ++i) {
        for (int i = 0; i < dispersalInertia; ++i) {
            if (DailyEnergy[i] >= DailyEnergy[i+1]) {
                return;
            }
        }
        
        // if we get to this point, that means we have decreasing energy for 10 consecutive days -> disperse!
        movementMode = directedDispersal;
        pickDispersalTarget();
        return;
        
    } else {
        
        // if energy today is higher than any day in the previous week, stop dispersing
        if (DailyEnergy[0] == *max_element(DailyEnergy.begin(), DailyEnergy.begin()+7)) {
            dispersalTarget.forget();
            movementMode = normalMove;
            return;
        }
        
        // if energy level was higher last week, disperse towards that area visited 7 days ago.
        if (movementMode == directedDispersal) {

            float energyRecent = std::accumulate(DailyEnergy.begin(), DailyEnergy.begin()+2, 0) / 3;
            float energyLastWeek = std::accumulate(DailyEnergy.begin()+6, DailyEnergy.begin()+9, 0) / 4;
    
            if (energyLastWeek > energyRecent) {
                Vector2df& xy = dailyPositions[7];
                int cell = sim->cellFromPoint(xy);
                if (cell != -1) {
                    int block = sim->Grid[cell].Block;
                    float dist = currentPos.distanceFrom(xy);
                    movementMode = returningDispersal;
                    dispersalTarget.set(block, xy, dist);
                    //Logger::debug(0, "porp %d is returning to area visited 7 days ago (%.02f, %.02f", Id, xy.x, xy.y);
                }
            }
        }
    }
}

void Porpoise::pickDispersalTarget(int excludeBlock) {

    int currentBlock = -1;

    if (currentCell != -1) {
        currentBlock = sim->Grid[currentCell].Block;
        //currentAbundanceBlock = sim->Grid[currentCell].abundanceBlock;
    } 
    
    std::vector<dispersalCandidate> blocks;
    blocks.reserve(sim->nBlocks); // we're never going to need more than this
    
    // select blocks for which the squared distance from the porpoise's current 
    // position to the block center exceeds the minimum squared dispersal distance,
    // but does not exceed maximum squared dispersal distance
    for (auto& block : sim->Blocks) {
        if (block.id() == currentBlock) continue; // skip current block
        float distSquared = currentPos.distanceFromSquared(block.center()); // in map units
        if (distSquared > minDispersalDistanceSquared && distSquared < maxDispersalDistanceSquared && (excludeBlock == -1 || block.id() != excludeBlock)) {
            float distEstimate = distSquared * (1.0f / sqrt(distSquared));
            blocks.emplace_back(block.id(), 
                                distEstimate,
                                block.value(sim->time->quarter()-1) / distEstimate);
        }
    }
    
    // sort blocks in order of decreasing quality, where quality = food density / num porps / distanceSquared
    std::sort(blocks.begin(), blocks.end(), [](const dispersalCandidate& i, const dispersalCandidate& j) {
        return i.value > j.value;
    });
    
    // select 3 best blocks (NB did 12)
    blocks.resize(12);
    
    // select one block at random from the 12 best blocks
    dispersalCandidate& block = sample<dispersalCandidate>(blocks);
    //auto& b = sim->Blocks[43];
    //float dist = pos.distanceFrom(b.center());
    //dispersalCandidate block{ b.id(), dist, b.value(sim->time->quarter()-1) / dist};
    Vector2df target = sim->Blocks[block.id].center();
    dispersalTarget.set(block.id, target, block.distance);

    //Logger::debug(0, "porp %d is dispering towards block %d (distance = %f, value = %.10g, best = %.10g, worst = %.10g, energy = %.3f)", 
    //    Id, block.id, block.distance, block.value, blocks[0].value, blocks[12].value, EnergyLevel);
    //}

}

void Porpoise::disperseTowardsTarget() {
    if ((movementMode != directedDispersal && movementMode != returningDispersal) || !dispersalTarget.isValid()) { return; }
    
    float currentDist = currentPos.distanceFrom(dispersalTarget.pos()); // distance from target

    bool stop = false;
    // check if porp should switch to dispersal mode 2

    if (dispersalStepCounter > 48 && currentPos.distanceFrom(dailyPositions[1]) < 2) { // porp has moved less than 0.8 km since yesterday 
        stop = true;
        //Logger::debug(0, "porp %d switched to coastal dispersal (moved less than 0.8 km in last 24 h)", Id);
    } else if (dispersalStepCounter > 432 &&currentPos.distanceFrom(dailyPositions[8]) < 6) { // porp has moved less than 2.4 km in the last week
        stop = true;
        //Logger::debug(0, "porp %d switched to coastal dispersal (moved less than 2.4 km since last week)", Id);
    } else if (currentDist < 50) {   // porp has crossed into target block
        stop = true;
        //Logger::debug(0, "porp %d arrived in its destination block!", Id);
    //} else if (currentCell != -1 && sim->Grid[cell].DistanceToCoast < minDispersalDistanceToLand) { // porp is moving too close to land
    //    stop = true;
    //    Logger::debug(0, "porp %d switched to coastal dispersal (moving too close to land)", Id);
    }
    
    if (stop) {
        dispersalTarget.forget();
        movementMode = coastalDispersal;
        disperseAlongCoast();
        return;
    }
    
    Vector2df dispStep = dispersalTarget.pos() - currentPos; // vector to destination (total delta in x and y coordinates)
    Vector2df mov = dispStep.normalize() * meanDispersalDistance; // rescale vector length to average dispersal distance

    // adjust dispersal direction by up to 30 degrees to swim towards deeper waters 
    // look 8 cells (2.4 km) ahead to make sure there is enough water 
    mov = sim->findPathDeepest(currentPos, mov, 30, 10, 8);
    
    // next adjust the (already adjusted) dispersal direction by up to 30 degrees to swim towards areas far from land
    mov = sim->findPathFarthestFromShore(currentPos, mov, 30, 10);
    
    // calculate new position
    Vector2df newPos = currentPos + mov;
    
    // check that the elected path is actually traversable, and if it's not, disperse along coast
    // but remember where we want to go!
    if (!sim->IsPathTraversable(currentPos.x, currentPos.y, newPos.x, newPos.y)) {
        //dispersalTarget.forget();
        disperseAlongCoast();
        return;
    }
    
    // execute the move
    float newHeading = 90 - atan2f(mov.y, mov.x) * 180.0 / PI;
    executeMove(newPos, newHeading, directedDispersal);
    
    // set dispersed to true (for energy calculations)
    dispersed = true;
}

// Try to stay 1-4 km from land, travelling away from yesterday's position
void Porpoise::disperseAlongCoast() {

    if (movementMode == normalMove) return;
    
    // vector pointing away from place visited 1 day ago
    Vector2df mov = (currentPos - dailyPositions[1]).normalize() * meanDispersalDistance; 
    
    // turn up to 80 degres in either direction (preferring smaller angles) to find a path
    // between 1 km (2.5 cells) and 4 km (10 cells) from the coast
    mov = sim->findPathParallellToCoast(currentPos, mov, 80, 10, 2.5, 10);
    
    // calculate new position
    Vector2df newPos = currentPos + mov;
    
    // check that the elected path is actually traversable, and if it's not, stop dispersing
    if (!sim->IsPathTraversable(currentPos.x, currentPos.y, newPos.x, newPos.y)) {
        //Logger::debug(0, "porp %d stopped coastal dispersal, because it couldn't find a traversable path (current pos = %.02f, %.02f)", Id, X[0], Y[0]);
        movementMode = normalMove;
        return;
    }
    
    // execute move
    float newHeading = 90 - atan2f(mov.y, mov.x) * 180.0 / PI;
    executeMove(newPos, newHeading, coastalDispersal);
    dispersed = true;
}








/***********************************************************************************************
 * 
 * 
 *                                   MATING AND REPRODUCTION
 * 
 * 
 ***********************************************************************************************/


void Porpoise::setMatingDay() {
    this->matingDay = sanitizeDayNumber(round(getRandomNormal(225, 20)));
}

void Porpoise::Mate() {
    if (Age >= AgeOfMaturity && isPregnant == false && getRandomFloat(0, 1) < pregnancy_prob) {
        isPregnant = true;
        calfBirthday = sanitizeDayNumber(sim->time->yday() - 65); // 10 months pregnancy (wrapping around the year)
    }
}

void Porpoise::giveBirth() {
    isPregnant = false;
    withCalf = true;
    weaningDay = sanitizeDayNumber(sim->time->yday() + 240); // nursing for 8 months
    calfBirthday = { -1 };
    setEnergyUse();
}

void Porpoise::weanCalf() {
    // create a new calf, and set its position to the same as its mother
    auto calf = std::unique_ptr<Porpoise>(new Porpoise(*this)); 
    Porpoises.push_back(std::move(calf));
    withCalf = { false };
    weaningDay = { -1 };
    setEnergyUse();
}

void Porpoise::abandonCalf() {
    withCalf = { false }; 
    weaningDay = { -1 };
    setEnergyUse();
}

/***********************************************************************************************
 * 
 * 
 *                                   GILLNET INTERACTION
 * 
 * 
 ***********************************************************************************************/


bool Porpoise::Entangled() {
    
    // list of gillnets in current cell
    auto& gillnets = sim->Grid[currentCell].gillnets;
    
    if (gillnets.size() == 0) {
        return false;
    }
    // check if the travelled path intersects with any gillnets
    for (Gillnet* gillnet : gillnets) {
        if (gillnet->check4(currentPos, currentCell)) {
            int currentBlock = sim->Grid[currentCell].fisheryBlock;
            sim->logger->log(sim->time->year(), sim->time->month(), sim->time->day(), currentBlock, gillnet->type(), 1, currentPos.x, currentPos.y);
            return true;
        }
    }
    
    return false;
}
