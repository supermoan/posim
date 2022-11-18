#ifndef __LOGGER__
#define __LOGGER__
#include <Rcpp.h>
#include <chrono>
#include <cstring>
#include <iostream>
#include <string>

#include "Porpoise.hpp"
#include "Timer.hpp"

class Settings;
class Porpoise;
class Gillnet;
class Timer;

extern int DEBUG_LEVEL;


class Logger {
private:
    Rcpp::List m_log;
    std::vector<int> m_follow { 0 }; // follow everyone
    std::vector<int> m_step;
    std::vector<int> m_id;
    std::vector<float> m_x;
    std::vector<float> m_y;
    std::vector<int> m_cell;
    std::vector<float> m_age;
    std::vector<float> m_prev_mov;
    std::vector<float> m_energy;
    std::vector<int> g_step;
    std::vector<int> g_gillnet;
    std::vector<float> g_x0;
    std::vector<float> g_y0;
    std::vector<float> g_x1;
    std::vector<float> g_y1;
    std::vector<int> g_type;
    std::vector<int> g_soaktime;
    //std::vector<float> g_length;
    std::vector<int> g_catch;
    std::vector<int> a_step;
    std::vector<int> a_N;
    std::vector<int> a_sets;
    //std::vector<unsigned long> a_effort;
    std::vector<int> a_bycatch;
    std::vector<float> a_food;
    std::vector<float> a_energy;
    std::vector<int> b_step;
    std::vector<int> b_ageClass;
    std::vector<int> b_N;
    std::vector<int> d_step;
    std::vector<int> d_block;
    std::vector<int> d_N;
    std::vector<int> d_births;
    std::vector<int> d_deaths;
    std::vector<int> d_bycatch;
    std::vector<int> z_year, z_month, z_day, z_block, z_type, z_count;
    std::vector<float> z_xpos, z_ypos;
    unsigned int bycatch_count = 0;
    unsigned long effort_count = 0;
    unsigned int set_count = 0;
    Rcpp::List Rlog;
public:
    Logger(int steps, std::vector<int> follow, int max_age) {
        if (follow.size() > 0 && follow[0] != -1) {
            int n = steps * follow.size();
            m_step.reserve(n);
            m_id.reserve(n);
            m_x.reserve(n);
            m_y.reserve(n);
            m_cell.reserve(n);
            m_age.reserve(n);
            m_prev_mov.reserve(n);
            m_energy.reserve(n);
        }
        m_follow = follow; // i think we must pass this by value
        
        int n = std::ceil(steps / 1488) +2;
        int medN = (n+2) * 11;
        int bigN = (n+2) * 3;
        
        a_step.reserve(n);
        a_N.reserve(n);
        a_sets.reserve(n);
        //a_effort.reserve(n);
        a_bycatch.reserve(n);
        a_food.reserve(n);
        a_energy.reserve(n);
        b_step.reserve(bigN);
        b_ageClass.reserve(bigN);
        b_N.reserve(bigN);
        d_step.reserve(medN);
        d_block.reserve(medN);
        d_N.reserve(medN);
    }
    void bycatch(int count);
    void gillnet_set(int count);
    void log(int step, const std::unique_ptr<Porpoise>& porp);
    void log(int step, Gillnet *gn);
    void log(int step, int ageClass, int N, int type);
    void log(int step, int abundanceRegion, int N, int births, int deaths, int bycatch);
    void log(int step, int N, float food, float energy);
    void log(int year, int month, int day, int block, int type, int count, float x, float y);
    // this must be defined here, because otherwise the function will have internal linkage
    template <typename... Args>
    static void debug(int level, const char* format, Args... args) {
        // has to be implement in header file; don't know why.https://open.spotify.com/show/5UTquk3eIW8u9xRV4tnt49
        if (DEBUG_LEVEL < level) return;
        
        std::string formatAsString = "";
        formatAsString += format;
        formatAsString += "\n";
        Rprintf(formatAsString.c_str(), args...);
    };
    
    void clear();
    Rcpp::List toList(); 
};

#endif // __LOGGER__
