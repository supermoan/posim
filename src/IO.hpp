#ifndef __SimIO__
#define __SimIO__

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "GridCell.hpp"
#include "Settings.hpp"
#include "Logger.h"
#include "Vector2d.hpp"
#include "Fishery.hpp"

class Settings;
class Fishery;

void read_gis(const std::string filename, Settings* sim);
int read_fishery_data(const std::string filename, Fishery& fishery);
int read_fishing_effort(const std::string filename, Fishery& fishery);
#endif // __SimIO__
