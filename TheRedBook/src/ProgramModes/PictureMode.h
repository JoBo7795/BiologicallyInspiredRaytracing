#pragma once
#include "GL/glew.h"

#include "../Shader/ShaderManager.h"
#include "../Renderer/Renderer.h"
#include "ProgramParams.h"
#include "RealtimeMode.h"
#include "../Camera/Camera.h"
#include "../Objects/Player.h"
#include "../Window/CreateWindow.h"
#include "../Window/InputControl.h"
#include "../Objects/LensData.h"
#include <iomanip>
#include <chrono>

class PictureMode
{
public:
	static void InitPictureMode(float in_lenseDistance, bool in_renderImageOnly,bool debugMode = false, int in_rDepth = 4);
	static const float infinity;
	static glm::vec3 imagePlanePos, imagePlaneDir;
	static bool debugMode;
	static Camera* cam,debugCam;
	static LensData lens;

	static int imageWidth;
	static int imageHeight;
	static int numRenders;
	static int numLoops;
	static int stepX;
	static int stepY;

	static void Render();

	static void SetDebugParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir);

};

