#pragma once
#ifndef __MISC__
#define __MISC__

#include <random>
#include "pcg_random.hpp"
#include "Vector2d.hpp"

typedef std::pair<Vector2df, Vector2df> _linestring;

// PCG-random from https://www.pcg-random.org/
extern pcg32 rng;

int getRandomInt(int min, int max);
float getRandomFloat(float min, float max);
int getRandomDiscrete(std::vector<int> *probs);
float getRandomNormal(float mean, float sd);
template<typename T> void shuffle(std::vector<T> const &x);

// adapted from 
// https://stackoverflow.com/questions/6942273/how-to-get-a-random-element-from-a-c-container
template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(rng));
    return start;
}

// randomly picks a single element from a vector and returns a reference to that element
template<typename T>
T& sample(std::vector<T>& vec) {
    auto start = vec.begin();
    std::uniform_int_distribution<> dis(0, std::distance(start, vec.end() - 1));
    std::advance(start, dis(rng));
    return *start;
}

bool intersects(_linestring A, _linestring B);
float distance(_linestring A, Vector2df P);
float distance(_linestring A, _linestring B);
float distanceSquared(_linestring A, Vector2df P);
float distanceSquared(_linestring A, _linestring B);

struct dispersalCandidate {
    int id;
    float distance;
    float value;
    dispersalCandidate() : id(0), distance(0), value(0) {}
    dispersalCandidate(int id, float distance, float value) : id(id), distance(distance), value(value) {}
        
};

// make sure day is between 1 and 365
int sanitizeDayNumber(int day);

#endif // __MISC__
