#ifndef __execHalfhourTasks__
#define __execHalfhourTasks__
#include <Rcpp.h>
#include <string>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <chrono>

// load global configuration first, as other classes depend on it.
#include "Timer.hpp"
#include "Logger.h"
#include "Settings.hpp"
#include "Vector2d.hpp"
#include "misc.hpp"
#include "GridCell.hpp"
#include "Porpoise.hpp"
#include "Gillnet.h"
#include "IO.hpp"

class Porpoise;
class Config;
class Logger;
class Timer;
class Settings;

void execHalfhourTasks(Settings& sim);

#endif // __execHalfhourTasks__
