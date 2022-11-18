#ifndef __PORPOISE__
#define __PORPOISE__
#include <algorithm>
#include <iostream>
#include "Settings.hpp"
#include "Vector2d.hpp"
#include "Gillnet.h"
#include "Block.hpp"
#include "Position.hpp"

extern int DEBUG_LEVEL;

class Settings;
class Gillnet;

struct PorpoiseState {
    Vector2df pos;
    float food;
    PorpoiseState() : pos(0, 0), food(0) {};
    PorpoiseState(Vector2df pos) : pos(pos), food(0) {};
    PorpoiseState(Vector2df pos, float food) : pos(pos), food(food) {};
};

class Porpoise {
private:
    static int nextId;
    static Settings* sim;
    friend class Settings;
public:
    enum PorpoiseMovementMode {
        normalMove = 0,
        directedDispersal = 1,
        coastalDispersal = 2,
        returningDispersal = 3
    };
    // static parameters, apply to all porps
    static std::vector<std::unique_ptr<Porpoise>> Porpoises;
    static float AgeOfMaturity;
    static float monthlyEnergyMultiplier[12];
    static float withCalfEnergyMultiplier;
    static float distEnergyMultiplier;
    static float stepEnergyMultiplier;
    static float m_mort_prob;
    static float x_surv_prob;
    static int dispersalInertia;
    static float meanDispersalDistance;
    static float minDispersalDistance;
    static float maxDispersalDistance;
    static float minDispersalDistanceSquared;
    static float maxDispersalDistanceSquared;
    static float minDispersalDepth;
    static float minDispersalDistanceToLand;
    static float pregnancy_prob;
    static float max_age;
    static int maxMemory;
    static float inertia_const;
    static float corrLogmov;
    static float corrAngle;
    static float m;
    static float maxLogmov;
    static std::vector<int> age_dist;
    static std::vector<float> ref_mem_strength; // reference memory decay rate: determines how fast animals forget the location of previously visited food patches
    static std::vector<float> work_mem_strength; // satiation memory decay rate; determines how fast the animals get hungry after eating
    static void resetId() { nextId = 0; } 
    static std::shared_ptr<Settings> Config;
    
    // individual porp properties

    int Id; // incremented automatically
    float Age = getRandomDiscrete(&age_dist); // age in decimal years, incremented daily by 1/365
    float Heading = getRandomFloat(0.0f, 359.9f); // randomly pick an initial direction (0 is north)
    int currentCell{ -1 };
    Vector2df currentPos{}, lastPos{};
    std::vector<PorpoiseState> track;
    std::vector<Vector2df> dailyPositions = std::vector<Vector2df>(10);
    
    bool isPregnant = false; // true for pregnant, false otherwise
    bool withCalf = false; // is the porp nursing a calf?
    int matingDay = { -1 }; // mating day
    int calfBirthday = { -1 }; // calf due day (if pregnant)
    int weaningDay = { -1 }; // calf weaning day (if with calf)
    
    float EnergyLevel = getRandomNormal(10, 1); // 0 - 20
    float cumulativeEnergy = 0;
    std::vector<float> DailyEnergy = std::vector<float>(10, 10);
    float E_use = 1; // rate of energy use

    Position dispersalTarget;    
    int movementMode{normalMove};
    bool dispersed{false}; // did porp disperse during its last action/"turn"? used for calculating energy use.
    int dispersalStepCounter{0};
    
    float prev_mov = { 6.309573 };
    float prev_angle = getRandomFloat(-25.0f, 25.0f);
    float prev_logmov = { 0.8 };
    float pres_mov, pres_logmov, pres_angle, CRW_contrib{ -9999 };
    
    // functions    
    Porpoise(int SurveyBlock); // constructor for initial population of porpoises (calves and non-calves)
    Porpoise(const Porpoise& mother); // constructor for calves that are born as the sim progresses

    Vector2df calcFoodAttractionVector();
    Vector2df calcCRW();
    float calcExpectedEnergy();
    void consumeFood();
    void setEnergyUse();
    void useEnergy();
    bool checkEnergy();
    void calcDailyEnergy();
    void abandonCalf();
    void considerDispersing();
    void pickDispersalTarget(int excludeBlock = -1);
    void disperseTowardsTarget();
    void disperseAlongCoast();
    void intrinsicMove();
    void executeMove(Vector2df newPos, float newHeading, PorpoiseMovementMode mode);
    void setMatingDay();
    bool Entangled();
    void Mate();
    void giveBirth();
    void weanCalf();
    void setHeading(float newHeading);
};

#endif // __PORPOISE__
