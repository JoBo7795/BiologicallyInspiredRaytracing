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
	void InitPictureMode(float in_lenseDistance, bool in_renderImageOnly,bool debugMode = false, int in_rDepth = 4);
	void Render();
	void SetDebugParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir);	
	static PictureMode* GetInstance();


private:
	static PictureMode* instance;
	int imageWidth;
	int imageHeight;
	int numRenders;
	int numLoops;
	int stepX;
	int stepY;

	const float infinity;
	glm::vec3 imagePlanePos, imagePlaneDir;
	bool debugMode;
	Camera* cam, debugCam;
	Renderer* rendererInstanceRef;
	LensData lens;
	PictureMode();
};

