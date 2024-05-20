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
private:

	static Camera* camera;
	static GLuint renderMode, queueObjectID;
	static std::vector<GameObject> RenderQueue, HUDRenderQueue;
	static std::vector<GLuint> RayTraceRenderQueue;
	static GLfloat lastFrame, deltaTime;	

	// variables for chunkRenderer
	static int r_translationModelLoc;
	static int r_rotationModelLoc;
	static int r_scaleModelLoc;
	static int r_customColorModelLoc;
	static int r_textureSetModelLoc;
	static int r_counter, r_chunkRenderCounter;
	static int r_renderQueueSize;
	static int lightPosModelLoc;
	static int viewModelLoc;
	static int projectionModelLoc;
	static int viewPosModelLoc;


	static void CalcDeltaTime();

public: 
	static void RenderHUD(Shader shaderProgram);
	static void RenderRayTrace(Shader shaderProgram);
	static void InitChunkRenderer(Shader shaderProgram);
	static void Render(Shader shaderProgram);	
	static GameObject* GetQueueObject(GLuint gameObjectID);
	static void DeleteObjectFromQueue(GLuint id);
	static void DeleteChunkObjectFromQueue(GLuint chunkObjectID);
	static GLuint AppendToRenderQueue(GameObject& gameObject);
	static GLuint AppendToRTRenderQueue(GLuint& gameObject);
	static GLuint AppendToHUDRenderQueue(GameObject gameObject);
	static std::vector<GameObject> GetRenderQueue();
	static void FlushQueue();
	static void FlushQueueWithVAO();
	static void FlushRTQueue();
	static void FlushHUDQueue();
	static void SetCamera(Camera* camera);
	static Camera* GetCamera();
	static void SetRenderMode(GLuint renderMode);
	static GLuint GetRenderMode();	
	static GLfloat DeltaTime();

	static glm::vec3 lightPos;
	static bool shadows_active;
	static bool debugMode;
	static float global_refract_index;
	static glm::vec2 factXY;
	static bool renderImageOnly;
	static float lenseDistance;
	static float lenseSizeX;
	static float lenseSizeY;
	static float lenseSizeZ;
	static float sceneBrightness;
	static float height;
	static float width;
	static int renderDepth;
	static glm::vec3 initCamPos;
	static glm::vec3 initlowerLeftCorner;
	static glm::vec3 initHorizontal;
	static glm::vec3 initVertical;
	static glm::vec3 lensPos;

};

