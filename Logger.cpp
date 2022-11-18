#include <Rcpp.h>
#include "Logger.h"
#include "Gillnet.h"
#include "Porpoise.hpp"

void Logger::log(int year, int month, int day, int block, int type, int count, float x, float y) {
    z_year.push_back(year);
    z_month.push_back(month);
    z_day.push_back(day);
    z_block.push_back(block);
    z_type.push_back(type);
    z_count.push_back(count);
    z_xpos.push_back(x);
    z_ypos.push_back(y);
}

void Logger::log(int step, const std::unique_ptr<Porpoise>& porp) {
    if ((m_follow.size() == 1 && m_follow[0] == 0) || std::find(m_follow.begin(), m_follow.end(), porp->Id) != m_follow.end()) {
        m_step.push_back(step);
        m_id.push_back(porp->Id);
        m_x.push_back(porp->track.front().pos.x);
        m_y.push_back(porp->track.front().pos.y);
        m_cell.push_back(porp->currentCell);
        m_prev_mov.push_back(pow(10, porp->prev_logmov)*100);
        m_age.push_back(porp->Age);
        m_energy.push_back(porp->EnergyLevel);
    }
}

void Logger::log(int step, Gillnet *gn) {
    auto coords = gn->get();
    g_step.push_back(gn->soaktime());
    g_x0.push_back(coords.first.x);
    g_y0.push_back(coords.first.y);
    g_x1.push_back(coords.second.x);
    g_y1.push_back(coords.second.y);
    g_type.push_back(gn->type());
    g_soaktime.push_back(gn->max_soaktime());
    //g_length.push_back(gn->length());
    g_catch.push_back(gn->bycatch());
}

void Logger::log(int step, int ageClass, int N, int type) {
    if (type == 1) {
        b_step.push_back(step);
        b_ageClass.push_back(ageClass);
        b_N.push_back(N);
    } else {
    }
}

void Logger::log(int step, int abundanceRegion, int N, int births, int deaths, int bycatch) {
    d_step.push_back(step);
    d_block.push_back(abundanceRegion);
    d_N.push_back(N);
    d_births.push_back(births);
    d_deaths.push_back(deaths);
    d_bycatch.push_back(bycatch);
}

void Logger::log(int step, int N, float food, float energy) {
    a_step.push_back(step);
    a_N.push_back(N);
    a_sets.push_back(set_count);
    //a_effort.push_back(effort_count);
    a_bycatch.push_back(bycatch_count);
    a_food.push_back(food);
    a_energy.push_back(energy);
    
    bycatch_count = 0;
    set_count = 0;
    //effort_count = 0;
}

void Logger::bycatch(int count) {
    bycatch_count += count;
}
void Logger::gillnet_set(int count) {
    set_count += count;
}

Rcpp::List Logger::toList() {

    using namespace Rcpp;
    
    return List::create(
        Named("porps") = DataFrame::create(
            Named("step") = wrap(m_step), 
            Named("id") = wrap(m_id),
            Named("age") = wrap(m_age),
            Named("x") = wrap(m_x), 
            Named("y") = wrap(m_y),
            Named("cell") = wrap(m_cell), 
            Named("prev_mov") = wrap(m_prev_mov),
            Named("energy") = wrap(m_energy)),

        Named("bycatch") = DataFrame::create(
            Named("year") = wrap(z_year), 
            Named("month") = wrap(z_month), 
            Named("day") = wrap(z_day),
            Named("block") = wrap(z_block),
            Named("type") = wrap(z_type), 
            Named("x") = wrap(z_xpos), 
            Named("y") = wrap(z_ypos), 
            Named("count") = wrap(z_count)),        
        
        Named("gillnets") = DataFrame::create(
            Named("step") = wrap(g_step),
            Named("x0") = wrap(g_x0),
            Named("y0") = wrap(g_y0),
            Named("x1") = wrap(g_x1),
            Named("y1") = wrap(g_y1),
            Named("type") = wrap(g_type),
            Named("soaktime") = wrap(g_soaktime),
            //Named("length") = wrap(g_length),
            Named("catch") = wrap(g_catch)),
            
        Named("N") = DataFrame::create(
            Named("step") = wrap(a_step), 
            Named("N") = wrap(a_N), 
            Named("sets") = wrap(a_sets),
           // Named("fishing_effort") = wrap(a_effort),
            Named("bycatch") = wrap(a_bycatch),
            Named("food") = wrap(a_food), 
            Named("energy") = wrap(a_energy)),
            
        Named("block") = DataFrame::create(
            Named("step") = wrap(d_step),
            Named("block") = wrap(d_block),
            Named("N") = wrap(d_N),
            Named("births") = wrap(d_births),
            Named("natural_mortality") = wrap(d_deaths),
            Named("bycatch_mortality") = wrap(d_bycatch)),
            
        Named("structure") = DataFrame::create(
            Named("step") = wrap(b_step), 
            Named("ageClass") = wrap(b_ageClass), 
            Named("N") = wrap(b_N))
    );
}