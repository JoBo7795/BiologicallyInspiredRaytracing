#include "ProgramParams.h"

float ProgramParams::windowWidth = 1920 / 2, ProgramParams::windowHeight = 1080 / 2;
float ProgramParams::posX = 0.0f, ProgramParams::posY = 0.0f, ProgramParams::posZ = 0.0f;
double ProgramParams::g;

std::string ProgramParams::path = (char*)"C:\\Users\\Johannes\\Desktop\\Renderings\\";
std::string ProgramParams::paramPaths;
bool ProgramParams::picOnly;
bool ProgramParams::CPURender;
float ProgramParams::lenseScale;
float ProgramParams::lenseDistance;
float ProgramParams::sceneBrightness=1.0f;
int ProgramParams::rDepth;
std::string ProgramParams::params;
bool ProgramParams::debugMode = false;
Shader ProgramParams::shaderProgram1;
Shader ProgramParams::shaderProgram2;
std::string ProgramParams::shaderName;
Shader ProgramParams::FragmentHUD, ProgramParams::tmpFragmentHUD;
Shader ProgramParams::currentHUDFragment = FragmentHUD;
Window ProgramParams::window;
double ProgramParams::b, ProgramParams::dpt, ProgramParams::f;