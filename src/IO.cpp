#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include "GridCell.hpp"
#include "Block.hpp"
#include "IO.hpp"
#include "Vector2d.hpp"

class GridCell;

void read_gis(const std::string filename, Settings* sim) {

    std::fstream fid (filename, std::ios::in);
    std::string line;
    bool header = false;
    int cellnum = 0;

    if(!fid.is_open()) {
        // throw error in R and stop further execution
        Rcpp::stop("read_gis: Could not open file for reading.");
        return;
    }
    
    while(getline(fid, line)) {
        std::stringstream str(line);
        std::string column;
        std::vector<std::string> data;
        
        if (line.empty()) {
            break;
        }
        
        while (std::getline(str, column, ';')) { 
            if (column != "") data.push_back(column); 
        }
        
        if (header == false) {
            if (data.size() != 12) {
                Rcpp::stop("Bad header in %s: wanted 12 columns, got %d", filename, data.size());
                return;
            }
            int xmx = stoi(data[0]);
            int ymx = stoi(data[1]);
            int nsurv = stoi(data[2]);
            int nfish = stoi(data[3]);
            int nblock = stoi(data[4]);
            int ntrav = stoi(data[5]);
            int bsize = stoi(data[6]);
            int npatch = stoi(data[7]);
            float meanMaxent1 = stof(data[8]);
            float meanMaxent2 = stof(data[9]);
            float meanMaxent3 = stof(data[10]);
            float meanMaxent4 = stof(data[11]);
            sim->initData(xmx, ymx, nsurv, nfish, nblock, ntrav, bsize, npatch, meanMaxent1, meanMaxent2, meanMaxent3, meanMaxent4),
            header = true;
            continue;
        }

        if (data.size() != 10) {
            Rcpp::stop("Bad data on line %d: wanted 10 columns, got %d", cellnum, data.size());
            return;
        }
        
        float averageDepth = stof(data[0]);
        float distToCoast = stof(data[1]);
        int abundanceRegion = stoi(data[2]) - 1;
        
        if (abundanceRegion < 0 || abundanceRegion >= sim->nSurveyBlocks) {
            abundanceRegion = -1;
        }
        int fisheryRegion = stoi(data[3]) - 1;
        int foodBlock = stoi(data[4]) -1;
        float foodLevel = stof(data[5]);
        float maxent1 =  stof(data[6]);
        float maxent2 =  stof(data[7]);
        float maxent3 =  stof(data[8]);
        float maxent4 =  stof(data[9]);
        
        bool cellTraversable = averageDepth <= sim->MinimumWaterDepth;
        bool containsFood = cellTraversable && foodLevel > 0.0f;
        if (containsFood) {
            foodLevel *= sim->maxU;
        } else {
            foodLevel = 0.0f;
        }
        
        Vector2df cellCoords = sim->GetXYFromCell(cellnum);
        float distToEdge = std::min({ cellCoords.x-0.5, 
                                    cellCoords.y-0.5, 
                                    sim->xmx - cellCoords.x + 0.5, 
                                    sim->ymx - cellCoords.y + 0.5 });

        sim->Grid.emplace_back(cellnum, 
                               foodBlock, 
                               -1, //abundanceRegion, 
                               fisheryRegion, 
                               averageDepth, 
                               distToCoast, 
                               distToEdge, 
                               foodLevel, 
                               maxent1, maxent2, maxent3, maxent4);

        if (cellTraversable) {
            sim->TraversableCells.push_back(cellnum);
            sim->Blocks[foodBlock].addCell(cellnum);
            
            if (containsFood) {
                sim->Patches.push_back(cellnum);
                sim->Blocks[foodBlock].addPatch(cellnum);
            }
            
            if (abundanceRegion != -1) {
                sim->abundanceRegions[abundanceRegion].addCell(cellnum);
            }
            
            if (distToCoast > 2.5f && averageDepth >= -400 && fisheryRegion >= 0 && fisheryRegion < sim->nFisheryBlocks) {
                sim->FisheryBlocks[fisheryRegion].push_back(cellnum);
            }
            
        }
        ++cellnum;
    }
    fid.close();
}

int read_fishery_data(const std::string filename, Fishery& fishery) {
    
    int bad_data_count = 0;
    int sample_size = 0;
    int nBlock = fishery.nBlock();
    int nDay = 365;
    int nType = fishery.nType();
    int nSample = fishery.nSample();
    std::fstream fid (filename, std::ios::in);
    std::string line;
    
    if(!fid.is_open()) {
        // throw error in R and stop further execution
        Rcpp::stop("read_fishery_data: Could not open file for reading.");
    }
    
    getline(fid, line); // read and discard first line (skip header)
    
    while(getline(fid, line)) {
        std::stringstream str(line);
        std::string column;
        std::vector<std::string> data;
        
        if (line.empty()) break;
        
        
        while (std::getline(str, column, ';')) { 
            data.push_back(column); 
        }
        
        int block = stoi(data[0]);
        int year = stoi(data[1]);
        int day = stoi(data[2]);
        int type = stoi(data[3]);
        int N = stoi(data[4]);
        bool pinger = stoi(data[5]);

        if (block >= 0 && block < nBlock 
                && day >= 0 && day < nDay 
                && type >= 0 && type < nType 
                && year >= 0 && year < nSample 
                && N > 0) 
        {
            fishery.addHauls(block, day, type, year, N, pinger);
            ++sample_size;
        } else if (N > 0) {
            ++bad_data_count;
        }
    }
    
    fid.close();

    if (bad_data_count > 0) {
        Rcpp::warning("%d rows contained out of range data, please verify data and settings", bad_data_count);
    }
    return sample_size;
}

int read_fishing_effort(const std::string filename, Fishery& fishery) {
    int nBlock = fishery.nBlock();
    int nSeason = 4;
    int nType = fishery.nType();
    int sample_size = 0;
    
    std::fstream fid (filename, std::ios::in);
    std::string line;
    
    if(!fid.is_open()) {
        // throw error in R and stop further execution
        Rcpp::stop("read_fishing_effort: Could not open file for reading.");
    }
    
    getline(fid, line); // read and discard first line (skip header)
    
    while(getline(fid, line)) {
        std::stringstream str(line);
        std::string column;
        std::vector<std::string> data;
        
        if (line.empty()) break;
        
        while (std::getline(str, column, ';')) { 
            data.push_back(column); 
        }
        
        int block = stoi(data[0]);
        int season = stoi(data[1]);
        int type = stoi(data[2]);
        int soaktime = stoi(data[3]);
        float length = stof(data[4]);
        
        if (block >= 0 && block < nBlock 
                && season >= 0 && season < nSeason 
                && type >= 0 && type < nType
                && soaktime >= 0
                && length > 0.0f) 
        {
            fishery.addEffort(block, season, type, soaktime, length);
            ++sample_size;
        }
    }
    
    fid.close();
    return sample_size;
}
