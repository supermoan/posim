#pragma once

/*
 * Fishing effort is made up of three components :
 *  soaktime: the total duration (in hours) that gillnets fish
 *  length: the total length (in map units) of the string of gillnets
 *  pinger: whether or not pingers are used
 */
struct FishingEffort {
    const int soaktime;
    const float length;
    bool pinger;
    FishingEffort() : soaktime(0), length(0), pinger(false) {};
    FishingEffort(int soaktime, float length) : soaktime(soaktime), length(length), pinger(false) {};
};
