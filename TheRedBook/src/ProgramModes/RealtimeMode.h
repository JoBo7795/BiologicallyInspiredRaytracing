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


	void SetDebugParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir);
	void SetCameraParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir);
	void Render();

	static RealtimeMode* GetInstance();

	glm::vec3 GetImagePlanePos();
	void SetImagePlanePos(glm::vec3& newPos);

	glm::vec3 GetImagePlaneDir();
	void SetImagePlaneDir(glm::vec3& newDir);

	bool IsDebugMode();
	void SetDebugMode(bool newDebugMode);

	Camera* GetCamera();
	void SetCamera(Camera* newCam);

	Camera GetDebugCamera();
	void SetDebugCamera(Camera newDebugCam);

	LensData GetLens();
	void SetLens(LensData& newLens);

	Renderer* GetRendererInstanceRef();
	void SetRendererInstanceRef(Renderer* newRenderer);

	// Getter and Setter for Renderer class parameters
	float RendererGetLensDistance();
	void RendererSetLensDistance(float lensDistance);

	bool RendererGetRenderImageOnly();
	void RendererSetRenderImageOnly(bool renderImageOnly);

	int RendererGetRenderDepth();
	void RendererSetRenderDepth(int renderDepth);

	
private:

	static RealtimeMode* instance;
	glm::vec3 imagePlanePos, imagePlaneDir;
	bool debugMode;
	Camera* cam, debugCam;
	LensData lens;
	Renderer* rendererInstanceRef;
	RealtimeMode();


};

