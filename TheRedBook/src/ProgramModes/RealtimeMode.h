#pragma once

#include "../Shader/ShaderManager.h"
#include "../Renderer/Renderer.h"
#include "ProgramParams.h"
#include "../Camera/Camera.h"
#include "../Objects/Player.h"
#include "../Window/CreateWindow.h"
#include "../Window/InputControl.h"
#include "../Objects/LensData.h"

class RealtimeMode
{
public:

	static glm::vec3 imagePlanePos, imagePlaneDir;
	static bool debugMode;
	static Camera* cam, debugCam;
	static LensData lens;

	static int imageWidth;
	static int imageHeight;

	static void InitPictureMode(float in_lenseDistance, bool in_debugMode, int in_rDepth = 4);
	static void SetDebugParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir);
	static void SetCameraParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir);
	static void Render();



};

