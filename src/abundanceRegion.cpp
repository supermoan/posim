
#include "abundanceRegion.hpp"

Settings* abundanceRegion::_sim {};

void abundanceRegion::addCell(const int cell) {
    _cells.push_back(cell);
}
bool abundanceRegion::empty() {
    return _cells.empty();
}
int abundanceRegion::randomCell() {
    return *select_randomly(_cells.begin(), _cells.end());
}
Vector2df abundanceRegion::randomPoint() {
    Vector2df pos{_sim->pointFromCell(randomCell())};
    pos.x += getRandomFloat(-0.49, 0.49);
    pos.y += getRandomFloat(-0.49, 0.49);
    return pos;
}
void abundanceRegion::reset_stats() {
    _births = 0;
    _deaths = 0;
    _bycatch = 0;
}

int abundanceRegion::id() {
    return _id;
}
void abundanceRegion::setId(int id) {
    _id = id;
}

