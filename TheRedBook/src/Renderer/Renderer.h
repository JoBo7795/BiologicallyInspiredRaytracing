#pragma once

#include "GL/glew.h"
#include "../Shader/Shader.h"
#include <vector>
#include <iostream>

#include "../../Dependencies/glm/glm/glm.hpp"
#include "../../Dependencies/glm/glm/gtc/matrix_transform.hpp"
#include "../../Dependencies/glm/glm/gtc/type_ptr.hpp"
#include "GLFW/glfw3.h"

#include "../Objects/GameObject.h"
#include "../Camera/Camera.h"

#include "../Objects/ObjectCreator.h"
#include "../Texture/TextureManager.h"
#include "../Window/Window.h"



class Renderer
{

public: 

	static Renderer* GetInstance();

	void RenderRayTraceShader(Shader shaderProgram);
	void RenderHUD(Shader shaderProgram);
	void Render(Shader shaderProgram);	
	GameObject* GetQueueObject(GLuint gameObjectID);
	void DeleteObjectFromQueue(GLuint id);
	void DeleteChunkObjectFromQueue(GLuint chunkObjectID);
	GLuint AppendToRenderQueue(GameObject& gameObject);
	GLuint AppendToRTRenderQueue(GLuint& gameObject);
	GLuint AppendToHUDRenderQueue(GameObject gameObject);
	std::vector<GameObject> GetRenderQueue();
	void FlushQueue();
	void FlushQueueWithVAO();
	void FlushRTQueue();
	void FlushHUDQueue();
	void SetCamera(Camera* camera);
	Camera* GetCamera();
	void SetRenderMode(GLuint renderMode);
	GLuint GetRenderMode();	
	GLfloat DeltaTime();




	int GetTranslationModelLoc();
	void SetTranslationModelLoc(int in_translationModelLoc);

	int GetRotationModelLoc();
	void SetRotationModelLoc(int in_rotationModelLoc);

	int GetScaleModelLoc();
	void SetScaleModelLoc(int in_scaleModelLoc);

	int GetCustomColorModelLoc();
	void SetCustomColorModelLoc(int in_customColorModelLoc);

	int GetTextureSetModelLoc();
	void SetTextureSetModelLoc(int in_textureSetModelLoc);

	int GetCounter();
	void SetCounter(int in_counter);

	int GetChunkRenderCounter();
	void SetChunkRenderCounter(int in_chunkRenderCounter);

	int GetRenderQueueSize();
	void SetRenderQueueSize(int in_renderQueueSize);

	int GetLightPosModelLoc();
	void SetLightPosModelLoc(int in_lightPosModelLoc);

	int GetViewModelLoc();
	void SetViewModelLoc(int in_viewModelLoc);

	int GetProjectionModelLoc();
	void SetProjectionModelLoc(int in_projectionModelLoc);

	int GetViewPosModelLoc();
	void SetViewPosModelLoc(int in_viewPosModelLoc);

	glm::vec3 GetLightPos();
	void SetLightPos(glm::vec3 in_lightPos);

	bool GetShadowsActive();
	void SetShadowsActive(bool in_shadows_active);

	bool GetDebugMode();
	void SetDebugMode(bool in_debugMode);

	float GetGlobalRefractIndex();
	void SetGlobalRefractIndex(float in_global_refract_index);

	glm::vec2 GetFactXY();
	void SetFactXY(glm::vec2 in_factXY);

	bool GetRenderImageOnly();
	void SetRenderImageOnly(bool in_renderImageOnly);

	float GetLenseDistance();
	void SetLenseDistance(float in_lenseDistance);

	int GetLenseSizeX();
	void SetLenseSizeX(float in_lenseSizeX);

	float GetLenseSizeY();
	void SetLenseSizeY(float in_lenseSizeY);

	float GetLenseSizeZ();
	void SetLenseSizeZ(float in_lenseSizeZ);

	float GetSceneBrightness();
	void SetSceneBrightness(float in_sceneBrightness);

	float GetHeight();
	void SetHeight(float in_height);

	float GetWidth();
	void SetWidth(float in_width);

	int GetRenderDepth();
	void SetRenderDepth(int in_renderDepth);

	glm::vec3 GetInitCamPos();
	void SetInitCamPos(glm::vec3 in_initCamPos);

	glm::vec3 GetInitlowerLeftCorner();
	void SetInitlowerLeftCorner(glm::vec3 in_initlowerLeftCorner);

	glm::vec3 GetInitHorizontal();
	void SetInitHorizontal(glm::vec3 in_initHorizontal);

	glm::vec3 GetInitVertical();
	void SetInitVertical(glm::vec3 in_initVertical);

	glm::vec3 GetLensPos();
	void SetLensPos(glm::vec3 in_lensPos);

private:

	static Renderer* instance;

	static Camera* camera;
	static GLuint renderMode, queueObjectID;
	static std::vector<GameObject> RenderQueue, HUDRenderQueue;
	static std::vector<GLuint> RayTraceRenderQueue;
	static GLfloat lastFrame, deltaTime;

	// GPU model locations
	int r_translationModelLoc;
	int r_rotationModelLoc;
	int r_scaleModelLoc;
	int r_customColorModelLoc;
	int r_textureSetModelLoc;
	int r_counter, r_chunkRenderCounter;
	int r_renderQueueSize;
	int r_lightPosModelLoc;
	int r_viewModelLoc;
	int r_projectionModelLoc;
	int r_viewPosModelLoc;


	glm::vec3 r_lightPos;
	bool r_shadows_active;
	bool r_debugMode;
	float r_global_refract_index;
	glm::vec2 r_factXY;
	bool r_renderImageOnly;
	float r_lenseDistance;
	float r_lenseSizeX;
	float r_lenseSizeY;
	float r_lenseSizeZ;
	float r_sceneBrightness;
	float r_height;
	float r_width;
	int r_renderDepth;
	glm::vec3 r_initCamPos;
	glm::vec3 r_initlowerLeftCorner;
	glm::vec3 r_initHorizontal;
	glm::vec3 r_initVertical;
	glm::vec3 r_lensPos;

	void CalcDeltaTime();
};

