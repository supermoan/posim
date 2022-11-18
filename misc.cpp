
#include <random>
#include "pcg_random.hpp"
#include "Vector2d.hpp"
#include "misc.hpp"

typedef std::pair<Vector2df, Vector2df> _linestring;

// PCG-random from https://www.pcg-random.org/
pcg_extras::seed_seq_from<std::random_device> seed_source;
pcg32 rng(seed_source);

int getRandomInt(int min, int max) {
    return std::uniform_int_distribution<int>{min, max}(rng);
}

float getRandomFloat(float min, float max) {
    return std::uniform_real_distribution<float>{min, max}(rng);
}

int getRandomDiscrete(std::vector<int> *probs) {
    std::discrete_distribution<int> dist(probs->begin(), probs->end());
    return dist(rng);
}

float getRandomNormal(float mean, float sd) {
    return std::normal_distribution<float>{mean, sd}(rng);
}

template<typename T> void shuffle(std::vector<T> const &x) {
    std::shuffle (x.begin(), x.end(), rng);
}


bool intersects(_linestring A, _linestring B) {
    float dx1 = A.second.x - A.first.x; // delta x, A 
    float dy1 = A.second.y - A.first.y; // delta y, A
    float dx2 = B.second.x - B.first.x; // delta x, B
    float dy2 = B.second.y - B.first.y; // delta y, B
    float det = dx2 * dy1 - dy2 * dx1;
    if (det < 0.001f) return false;  // segments are parallel
    float s = (dx1 * (B.first.y - A.first.y) + dy1 * (A.first.x - B.first.x)) / det;
    float t = (dx2 * (A.first.y - B.first.y) + dy2 * (B.first.x - A.first.x)) / (-det);
    return (s >= 0 && s <= 1) && (t >= 0 && t <= 1);
}


float distance(_linestring A, _linestring B) {
    if (intersects(A, B)) {
        return 0;
    }
    
    std::vector<float> dist;
    dist.push_back(distanceSquared(A, B.first));
    dist.push_back(distanceSquared(A, B.second));
    dist.push_back(distanceSquared(B, A.first));
    dist.push_back(distanceSquared(B, A.second));
    return sqrt(*min_element(dist.begin(), dist.end()));
}
float distanceSquared(_linestring A, _linestring B) {
    if (intersects(A, B)) {
        return 0;
    }
    
    std::vector<float> dist;
    dist.push_back(distanceSquared(A, B.first));
    dist.push_back(distanceSquared(A, B.second));
    dist.push_back(distanceSquared(B, A.first));
    dist.push_back(distanceSquared(B, A.second));
    return *min_element(dist.begin(), dist.end());
}
/*
 float distance(_linestring A, Vector2df P) {
 return sqrt(distanceSquared(A, P));
 }
 
 */



float distance(_linestring A, Vector2df P) {
    
    float dx = A.second.x - A.first.x;
    float dy = A.second.y - A.first.y;
    
    float t = ((P.x - A.first.x) * dx + (P.y - A.first.y) * dy) / (dx * dx + dy * dy);
    if (t < 0) {
        dx = P.x - A.first.x;
        dy = P.y - A.first.y;
    } else if (t > 1) {
        dx = P.x - A.second.x;
        dy = P.y - A.second.y;
    } else {
        dx = P.x - (A.first.x + t*dx);
        dy = P.y - (A.first.y + t*dy);
    }
    
    return sqrt(dx*dx + dy*dy);
}
float distanceSquared(_linestring A, Vector2df P) {
    
    float dx = A.second.x - A.first.x;
    float dy = A.second.y - A.first.y;
    
    float t = ((P.x - A.first.x) * dx + (P.y - A.first.y) * dy) / (dx * dx + dy * dy);
    if (t < 0) {
        dx = P.x - A.first.x;
        dy = P.y - A.first.y;
    } else if (t > 1) {
        dx = P.x - A.second.x;
        dy = P.y - A.second.y;
    } else {
        dx = P.x - (A.first.x + t*dx);
        dy = P.y - (A.first.y + t*dy);
    }
    
    return dx*dx + dy*dy;
}

int sanitizeDayNumber(int day) {
    if (day > 365) return day - 365;
    if (day <= 0) return day + 365;
    return day;
}
