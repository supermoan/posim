#pragma once

#include "Settings.hpp"
#include "FishingEffort.hpp"

/* 
 * abstraction to simplify sampling from historical fishing effort
 */
class Fishery {
private:
    int _nBlock, _nType, _nSample;
    std::vector<std::vector<std::vector<std::vector<FishingEffort>>>> _effort;
    std::vector<std::vector<std::vector<std::vector<int>>>> _hauls;
    std::vector<std::vector<std::vector<bool>>> _pinger;
public:
    void initialize(int nBlock, int nType, int nSample) {
        _nBlock = nBlock;
        _nType = nType;
        _nSample = nSample;
        _effort.resize(_nBlock, std::vector<std::vector<std::vector<FishingEffort>>>(4, std::vector<std::vector<FishingEffort>>(nType)));
        _hauls.resize(_nBlock, std::vector<std::vector<std::vector<int>>>(365, std::vector<std::vector<int>>(nType, std::vector<int>(_nSample, 0))));
        _pinger.resize(_nBlock, std::vector<std::vector<bool>>(365, std::vector<bool>(nType, false)));
    }
    void addHauls(int block, int day, int type, int year, int count, bool pinger) {
        _hauls[block][day][type][year] = count;
        _pinger[block][day][type] = pinger;
    }
    void addEffort(int block, int season, int type, int soaktime, float length) {
        _effort[block][season][type].emplace_back(soaktime, length);
    }
    void removeBlock(int block) {
        _hauls.erase(_hauls.begin()+block);
        _effort.erase(_effort.begin()+block);
        _pinger.erase(_pinger.begin()+block);
        _nBlock--;
    }
    int sampleNumberOfSets(int block, int day, int season, int type) {
        if (block == -1 || block >= _nBlock || day < 0 || day > 364 || type < 0 || type >= _nType || season < 0 || season > 3) { 
            return 0;
        }
        if (_hauls[block][day][type].size() == 0 || _effort[block][season][type].size() == 0) {
            return 0;
        }
        return *select_randomly(_hauls[block][day][type].begin(), _hauls[block][day][type].end());
    }
    FishingEffort sampleEffort(int block, int day, int season, int type) {
        FishingEffort ret = *select_randomly(_effort[block][season][type].begin(), _effort[block][season][type].end());
        ret.pinger = _pinger[block][day][type];
        return ret;
    }
    int nBlock() { return _nBlock; }
    int nType() { return _nType; }
    int nSample() { return _nSample; }
};