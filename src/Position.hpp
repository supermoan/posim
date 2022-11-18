#ifndef __POSITION__
#define __POSITION__
#include "Vector2d.hpp"

class Position {
private:
    int m_block; // index of block in Blocks vector
    Vector2df m_pos;
    std::vector<float> m_distance;
public:
    Position() : m_block(-1), m_pos(Vector2df(0,0)) {
        m_distance.reserve(9);
    };
    Vector2df pos() { return m_pos; };
    int block() { return m_block; };
    bool isValid() { return m_block != -1; }
    float x() { return m_pos.x; }
    float y() { return m_pos.y; }
    void forget() {
        m_block = -1;
        m_distance.clear();
        m_pos.x = 0;
        m_pos.y = 0;
    }
    void set(int block, Vector2df pos, float dist) {
        m_block = block;
        m_pos = pos;
        m_distance.clear();
        m_distance.push_back(dist);
    }
    void set(Vector2df pos, float dist) {
        m_block = 0;
        m_pos = pos;
        m_distance.clear();
        m_distance.push_back(dist);
    }
    void updateDistance(float dist) {
        m_distance.insert(m_distance.begin(), dist);
        m_distance.resize(9);
    }
    void updateDistance(Vector2df pos) {
        int dist = pos.distanceFrom(m_pos);
        updateDistance(dist);
    }
    bool gettingCloser() {
        if (m_distance.size() < 9) return true;
        
        float dist[3];
        for (int i = 0; i < 3; ++i) {
            int j = 3*(i-1);
            dist[i] = std::accumulate(m_distance.begin() + j, m_distance.begin() + j + 2, 0);
        }
        return dist[0] < dist[1] && dist[1] < dist[2];
    }
};


#endif // __POSITION__

