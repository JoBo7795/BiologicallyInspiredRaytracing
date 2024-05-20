#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>

#include "../ShapeData/ShapeData.h"
#include "../Mathematics/Mathematics.h"

namespace DataReader
{
	ShapeData readPLY(std::string path);
};

