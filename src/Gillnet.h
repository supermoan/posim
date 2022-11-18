#ifndef __GILLNET__
#define __GILLNET__
#include <Rcpp.h>
#include <list>
#include "Vector2d.hpp"
#include "Settings.hpp"
#include "misc.hpp"
#include "Logger.h"

class Settings;

class Gillnet {
    static Settings* sim;
    static float m_catchability_by_type[3];
    friend class Settings;
    std::vector<int> m_cellnums; // which cells are covered by this line segment?
    std::vector<std::list<Gillnet*>::iterator> m_cells;
    std::pair<Vector2df, Vector2df> m_coords;
    static std::vector<float> interaction_probability;
    float m_length { 0.0f }; // total length of gillnet string in units of 400m
    int m_block;
    int m_type; // gillnet mesh size, 0 = small, 1 = medium, 2 = large
    int m_catch = { 0 }; // number of porpoises caught
    float m_soaktime = { 0.0f }; // elapsed soaktime in hours
    int m_max_soaktime; // soaktime in hours
    float m_catchability; // harbour porpoise catchability
    bool m_pingered { false }; // does the gillnet have pingers?
public:
    Gillnet(int block, int type, int soaktime, float length, bool pingered);
    ~Gillnet();
    std::pair<Vector2df, Vector2df> get() { return m_coords; }
    int soaktime() { return m_soaktime; }
    int max_soaktime() { return m_max_soaktime; }
    void soak24h() { ++m_soaktime; }
    void soak30m() { m_soaktime += 0.021; }
    bool haulable() { return m_soaktime >= m_max_soaktime; }
    int bycatch() { return m_catch; }
    bool check(std::pair<Vector2df, Vector2df> path);
    bool check2(std::pair<Vector2df, Vector2df> path);
    bool check3(Vector2df xy, int cell);
    bool check4(Vector2df xy, int cell);
    int type() { return m_type; }
    float length() { return m_length; }
};

#endif // __GILLNET__
