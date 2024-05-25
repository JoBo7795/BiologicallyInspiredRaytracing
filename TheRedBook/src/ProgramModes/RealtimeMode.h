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


	int imageWidth;
	int imageHeight;

	void InitPictureMode(float in_lenseDistance, bool in_debugMode, int in_rDepth = 4);
	void SetDebugParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir);
	void SetCameraParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir);
	void Render();

	static RealtimeMode* GetInstance();
	
private:

	static RealtimeMode* instance;
	glm::vec3 imagePlanePos, imagePlaneDir;
	bool debugMode;
	Camera* cam, debugCam;
	LensData lens;
	Renderer* rendererInstanceRef;
	RealtimeMode();


};

