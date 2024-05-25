#include "Renderer.h"


Camera* Renderer::camera = new Camera(glm::vec3(0.0f, 0.0f, 0.0f));
GLuint Renderer::renderMode = GL_FILL, Renderer::queueObjectID = 0;
std::vector<GameObject> Renderer::RenderQueue;
std::vector<GameObject> Renderer::HUDRenderQueue;
std::vector<GLuint> Renderer::RayTraceRenderQueue;
GLfloat Renderer::lastFrame = 0.0f, Renderer::deltaTime = 0.0f;

// TODO Need to refactor class

 // int Renderer::r_translationModelLoc;
 // int Renderer::r_rotationModelLoc;
 // int Renderer::r_scaleModelLoc;
 // int Renderer::r_customColorModelLoc;
 // int Renderer::r_textureSetModelLoc;
 // int Renderer::r_counter, Renderer::r_chunkRenderCounter;
 // int Renderer::r_renderQueueSize;
 // int Renderer::lightPosModelLoc;
 // int Renderer::viewModelLoc;
 // int Renderer::projectionModelLoc;
 // int Renderer::viewPosModelLoc;

// float Renderer::lenseDistance = 0;//291.68f;
// float Renderer::lenseSizeX = 1.0f;
// float Renderer::lenseSizeY = 1.0f;
// float Renderer::lenseSizeZ = 1.0f;
// float Renderer::sceneBrightness = 1.0f;
// glm::vec3 Renderer::lightPos = glm::vec3(1000.0f, -1000.0f, 0.0f);
// bool Renderer::shadows_active = false;
// bool Renderer::debugMode = false;
// float Renderer::global_refract_index = 1.59f;
// glm::vec2 Renderer::factXY;
// bool Renderer::renderImageOnly = false;
// float Renderer::height = 1920;
// float Renderer::width = 1080;
// int Renderer::renderDepth = 3;
// glm::vec3 Renderer::initCamPos;
// glm::vec3 Renderer::initlowerLeftCorner;
// glm::vec3 Renderer::initHorizontal;
// glm::vec3 Renderer::initVertical;
// glm::vec3 Renderer::lensPos = glm::vec3(0.0);
 Renderer* Renderer::instance;


Renderer* Renderer::GetInstance() {
	if (!instance) {
		instance = new Renderer();
	}
	return instance;
}


void Renderer::RenderRayTraceShader(Shader shaderProgram) {

	shaderProgram.Use();
	int shaderID = shaderProgram.GetProgramID();
	glm::vec3 view = camera->GetPosition() + camera->GetDirection();
	glm::vec3 viewPos = camera->GetPosition();
	glm::vec3 up = camera->GetUp();

	for (GameObject& gameObject : HUDRenderQueue) {

		DrawObjectManager::UseDrawObject(gameObject.GetDrawObjectID());

		for (GLuint id : gameObject.GetTextureIDs()) {
			TextureManager::UseTexture(id);
		}

		int modelLoc = glGetUniformLocation(shaderID, "view");
		glUniform3f(modelLoc, view.x, view.y, view.z);

		modelLoc = glGetUniformLocation(shaderID, "viewPos");
		glUniform3f(modelLoc, viewPos.x, viewPos.y, viewPos.z);

		modelLoc = glGetUniformLocation(shaderID, "up");
		glUniform3f(modelLoc, up.x, up.y, up.z);

		modelLoc = glGetUniformLocation(shaderID, "lensPos");
		glUniform3f(modelLoc, r_lensPos.x, r_lensPos.y, r_lensPos.z);

		modelLoc = glGetUniformLocation(shaderID, "lightPos");
		glUniform3f(modelLoc, r_lightPos.x, r_lightPos.y, r_lightPos.z);

		modelLoc = glGetUniformLocation(shaderID, "shadows_active");
		glUniform1i(modelLoc, r_shadows_active);

		modelLoc = glGetUniformLocation(shaderID, "debugMode");
		glUniform1i(modelLoc, r_debugMode);

		modelLoc = glGetUniformLocation(shaderID, "global_refract_index");
		glUniform1f(modelLoc, r_global_refract_index);

		modelLoc = glGetUniformLocation(shaderID, "factXY");
		glUniform2f(modelLoc, r_factXY.x, r_factXY.y);

		modelLoc = glGetUniformLocation(shaderID, "widthHeight");
		glUniform2f(modelLoc, r_width, r_height);

		modelLoc = glGetUniformLocation(shaderID, "renderImageOnly");
		glUniform1f(modelLoc, r_renderImageOnly);

		modelLoc = glGetUniformLocation(shaderID, "renderDepth");
		glUniform1i(modelLoc, r_renderDepth);
		
		glm::vec3 transView = glm::vec3(0.0f,0.0f,0.0f) * glm::normalize(camera->GetDirection()) * r_lenseDistance;
		auto lenseTransform = camera->GetViewMat();

		glm::vec4 posL = (camera->GetViewMat() * glm::vec4(0.0, 0.0, 0.0, 1.0));

		glm::vec3 normalizedDir = -posL;
		glm::vec3 lensOrig = camera->GetPosition() + normalizedDir;

		modelLoc = glGetUniformLocation(shaderID, "lenseTransform");
		glUniformMatrix4fv(modelLoc,1,false, &lenseTransform[0][0]);	

		modelLoc = glGetUniformLocation(shaderID, "initCamPos");
		glUniform3f(modelLoc, r_initCamPos.x, r_initCamPos.y, r_initCamPos.z);

		modelLoc = glGetUniformLocation(shaderID, "initlowerLeftCorner");
		glUniform3f(modelLoc, r_initlowerLeftCorner.x , r_initlowerLeftCorner.y, r_initlowerLeftCorner.z);

		modelLoc = glGetUniformLocation(shaderID, "initHorizontal");
		glUniform3f(modelLoc, r_initHorizontal.x, r_initHorizontal.y, r_initHorizontal.z);

		modelLoc = glGetUniformLocation(shaderID, "initVertical");
		glUniform3f(modelLoc, r_initVertical.x, r_initVertical.y, r_initVertical.z);
				
		modelLoc = glGetUniformLocation(shaderID, "sceneBrightness");
		glUniform1f(modelLoc, r_sceneBrightness);

		glDrawElements(GL_TRIANGLES, gameObject.GetIndiSize(), GL_UNSIGNED_INT, 0);
	}
}

void Renderer::RenderHUD(Shader shaderProgram) {

	shaderProgram.Use();

	for (GameObject& gameObject : HUDRenderQueue) {

		glDrawElements(GL_TRIANGLES, gameObject.GetIndiSize(), GL_UNSIGNED_INT, 0);

	}
}

void Renderer::Render(Shader shaderProgram) {

	shaderProgram.Use();

	glPolygonMode(GL_FRONT_AND_BACK, renderMode);

	glm::mat4 proj = glm::perspective(glm::radians(camera->GetZoom()), 1920.0f / 1080.0f, 0.1f, 100.0f);
	glm::mat4 view = camera->LookAt();
	glm::vec3 viewPos = camera->GetPosition();
	glm::vec3 viewDir = camera->GetDirection();

	int shaderID = shaderProgram.GetProgramID();

	int modelLoc = glGetUniformLocation(shaderID, "view");
	glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(view));

	modelLoc = glGetUniformLocation(shaderID, "projection");
	glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(proj));

	modelLoc = glGetUniformLocation(shaderID, "viewPos");
	glUniform3f(modelLoc, viewPos.x, viewPos.y, viewPos.z);

	modelLoc = glGetUniformLocation(shaderID, "dir");
	glUniform3f(modelLoc, viewDir.x, viewDir.y, viewDir.z);

	int renderQueueSize = RenderQueue.size();

	int translationModelLoc = glGetUniformLocation(shaderID, "translation");
	int rotationModelLoc = glGetUniformLocation(shaderID, "rotation");
	int scaleModelLoc = glGetUniformLocation(shaderID, "scale");
	int customColorModelLoc = glGetUniformLocation(shaderID, "customObjectColor");
	int textureSetModelLoc = glGetUniformLocation(shaderID, "textureSet");

	for (int i = 0; i < renderQueueSize; i++) {

		DrawObjectManager::UseDrawObject(RenderQueue[i].GetDrawObjectID());

		if (RenderQueue[i].GetTextureSet()) {

			for (GLuint id : RenderQueue[i].GetTextureIDs()) {				
				TextureManager::UseTexture(id);
			}			
			glUniform1i(textureSetModelLoc, true);
		}
		else {
			glUniform1i(textureSetModelLoc, false);
		}

		glm::vec3 objectColor = RenderQueue[i].GetColor();
		glUniform3f(customColorModelLoc, objectColor.x, objectColor.y, objectColor.z);

		glUniformMatrix4fv(translationModelLoc, 1, GL_FALSE, glm::value_ptr(RenderQueue[i].GetTranslation()));

		glUniformMatrix4fv(rotationModelLoc, 1, GL_FALSE, glm::value_ptr(RenderQueue[i].GetRotation()));

		glUniformMatrix4fv(scaleModelLoc, 1, GL_FALSE, glm::value_ptr(RenderQueue[i].GetScale()));

		glDrawElements(GL_TRIANGLES, DrawObjectManager::GetCurrentBound().shapeData.indexSize, GL_UNSIGNED_INT, 0);

	}

}

GameObject* Renderer::GetQueueObject(GLuint gameObjectID) {

	for (int i = 0; i < RenderQueue.size(); i++) {
		// -1 possible source of failure
		if (RenderQueue[i].GetObjectID() == gameObjectID-1) {
			return &RenderQueue[i];
		}
	}
	return nullptr;
}

void Renderer::DeleteObjectFromQueue(GLuint gameObjectID) {

	for (int i = 0; i < RenderQueue.size(); i++) {
		if (RenderQueue[i].GetObjectID() == gameObjectID) {
			RenderQueue.erase(RenderQueue.begin() + i);
		}
	}
}


void Renderer::DeleteChunkObjectFromQueue(GLuint chunkObjectID) {

	for (int i = 0; i < RenderQueue.size(); i++) {
		if (RenderQueue[i].GetChunkObjectID() == chunkObjectID) {
			RenderQueue.erase(RenderQueue.begin() + i);
		}
	}
}

GLuint Renderer::AppendToRenderQueue(GameObject& gameObject) {

	RenderQueue.push_back(gameObject);
	RenderQueue[RenderQueue.size() - 1].SetObjectID(queueObjectID++);
	return queueObjectID;
}

GLuint Renderer::AppendToRTRenderQueue(GLuint& gameObject) {

	RayTraceRenderQueue.push_back(gameObject);
	
	return queueObjectID;
}

GLuint Renderer::AppendToHUDRenderQueue(GameObject gameObject) {
	HUDRenderQueue.push_back(gameObject);
	return HUDRenderQueue.size();
}

std::vector<GameObject> Renderer::GetRenderQueue() {
	return RenderQueue;
}

void Renderer::FlushQueue() {

	RenderQueue.clear();
}

void Renderer::FlushQueueWithVAO() {

	//for (GameObject gameObject : RenderQueue)
		//if (glIsBuffer(gameObject.GetBuffer().GetVAO()))
			//gameObject.GetBuffer().DeleteBuffer();

	//RenderQueue.clear();
}

void Renderer::FlushRTQueue() {

	RayTraceRenderQueue.clear();
}

void Renderer::FlushHUDQueue() {

	HUDRenderQueue.clear();
}




void Renderer::SetCamera(Camera* camera) {
	Renderer::camera = camera;
}

Camera* Renderer::GetCamera() {
	return Renderer::camera;
}

void Renderer::SetRenderMode(GLuint renderMode) {
	Renderer::renderMode = renderMode;
}

GLuint Renderer::GetRenderMode() {
	return Renderer::renderMode;
}

void Renderer::CalcDeltaTime() {

	float currentFrame = (float)glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;
}

GLfloat Renderer::DeltaTime() {
	return deltaTime;
}


int Renderer::GetTranslationModelLoc() {
	return r_translationModelLoc;
}

void Renderer::SetTranslationModelLoc(int in_translationModelLoc) {
	this->r_translationModelLoc = in_translationModelLoc;
}



int Renderer::GetRotationModelLoc() {
	return r_rotationModelLoc;
}

void Renderer::SetRotationModelLoc(int in_rotationModelLoc) {
	this->r_rotationModelLoc = in_rotationModelLoc;
}




int Renderer::GetScaleModelLoc() {
	return r_scaleModelLoc;
}

void Renderer::SetScaleModelLoc(int in_scaleModelLoc) {
	this->r_scaleModelLoc = in_scaleModelLoc;
}



int Renderer::GetCustomColorModelLoc() {
	return r_customColorModelLoc;
}

void Renderer::SetCustomColorModelLoc(int in_customColorModelLoc) {
	this->r_customColorModelLoc = in_customColorModelLoc;
}



int Renderer::GetTextureSetModelLoc() {
	return r_textureSetModelLoc;
}

void Renderer::SetTextureSetModelLoc(int in_textureSetModelLoc) {
	this->r_textureSetModelLoc = in_textureSetModelLoc;
}



int Renderer::GetCounter() {
	return r_counter;
}

void Renderer::SetCounter(int in_counter) {
	this->r_counter = in_counter;
}



int Renderer::GetChunkRenderCounter () {
	return r_chunkRenderCounter;
}

void Renderer::SetChunkRenderCounter(int in_chunkRenderCounter) {
	this->r_chunkRenderCounter = in_chunkRenderCounter;
}



int Renderer::GetRenderQueueSize() {
	return r_renderQueueSize;
}

void Renderer::SetRenderQueueSize(int in_renderQueueSize) {
	this->r_renderQueueSize = in_renderQueueSize;
}




int Renderer::GetLightPosModelLoc() {
	return r_lightPosModelLoc;
}

void Renderer::SetLightPosModelLoc(int in_lightPosModelLoc) {
	this->r_lightPosModelLoc = in_lightPosModelLoc;
}



int Renderer::GetViewModelLoc() {
	return r_viewModelLoc;
}

void Renderer::SetViewModelLoc(int in_viewModelLoc) {
	this->r_viewModelLoc = in_viewModelLoc;
}



int Renderer::GetProjectionModelLoc() {
	return r_projectionModelLoc;
}

void Renderer::SetProjectionModelLoc(int in_projectionModelLoc) {
	this->r_projectionModelLoc = in_projectionModelLoc;
}



int Renderer::GetViewPosModelLoc() {
	return r_viewPosModelLoc;
}

void Renderer::SetViewPosModelLoc(int in_viewPosModelLoc) {
	this->r_viewPosModelLoc = in_viewPosModelLoc;
}



glm::vec3 Renderer::GetLightPos() {
	return r_lightPos;
}

void Renderer::SetLightPos(glm::vec3 in_lightPos) {
	this->r_lightPos = in_lightPos;
}



bool Renderer::GetShadowsActive() {
	return r_shadows_active;
}

void Renderer::SetShadowsActive(bool in_shadows_active) {
	this->r_shadows_active = in_shadows_active;
}



bool Renderer::GetDebugMode() {
	return r_viewPosModelLoc;
}

void Renderer::SetDebugMode(bool in_debugMode) {
	this->r_debugMode = in_debugMode;
}



float Renderer::GetGlobalRefractIndex() {
	return r_global_refract_index;
}

void Renderer::SetGlobalRefractIndex(float in_global_refract_index) {
	this->r_global_refract_index = in_global_refract_index;
}



glm::vec2 Renderer::GetFactXY() {
	return r_factXY;
}

void Renderer::SetFactXY(glm::vec2 in_factXY) {
	this->r_factXY = in_factXY;
}



bool Renderer::GetRenderImageOnly() {
	return r_renderImageOnly;
}

void Renderer::SetRenderImageOnly(bool in_renderImageOnly) {
	this->r_renderImageOnly = in_renderImageOnly;
}


float Renderer::GetLenseDistance() {
	return r_lenseDistance;
}

void Renderer::SetLenseDistance(float in_lenseDistance) {
	this->r_lenseDistance = in_lenseDistance;
}



int Renderer::GetLenseSizeX() {
	return r_lenseSizeX;
}

void Renderer::SetLenseSizeX(float in_lenseSizeX) {
	this->r_lenseSizeX = in_lenseSizeX;
}



float Renderer::GetLenseSizeY() {
	return r_lenseSizeY;
}

void Renderer::SetLenseSizeY(float in_lenseSizeY) {
	this->r_lenseSizeY = in_lenseSizeY;
}



float Renderer::GetLenseSizeZ() {
	return r_lenseSizeZ;
}

void Renderer::SetLenseSizeZ(float in_lenseSizeZ) {
	this->r_lenseSizeZ = in_lenseSizeZ;
}



float Renderer::GetSceneBrightness() {
	return r_sceneBrightness;
}

void Renderer::SetSceneBrightness(float in_sceneBrightness) {
	this->r_sceneBrightness = in_sceneBrightness;
}



float Renderer::GetHeight() {
	return r_height;
}

void Renderer::SetHeight(float in_height) {
	this->r_height = in_height;
}



float Renderer::GetWidth() {
	return r_width;
}

void Renderer::SetWidth(float in_width) {
	this->r_width = in_width;
}



int Renderer::GetRenderDepth() {
	return r_renderDepth;
}

void Renderer::SetRenderDepth(int in_renderDepth) {
	this->r_renderDepth = in_renderDepth;
}



glm::vec3 Renderer::GetInitCamPos() {
	return r_initCamPos;
}

void Renderer::SetInitCamPos(glm::vec3 in_initCamPos) {
	this->r_initCamPos = in_initCamPos;
}



glm::vec3 Renderer::GetInitlowerLeftCorner() {
	return r_initlowerLeftCorner;
}

void Renderer::SetInitlowerLeftCorner(glm::vec3 in_initlowerLeftCorner) {
	this->r_initlowerLeftCorner = in_initlowerLeftCorner;
}



glm::vec3 Renderer::GetInitHorizontal() {
	return r_initHorizontal;
}

void Renderer::SetInitHorizontal(glm::vec3 in_initHorizontal) {
	this->r_initHorizontal = in_initHorizontal;
}



glm::vec3 Renderer::GetInitVertical() {
	return r_initVertical;
}

void Renderer::SetInitVertical(glm::vec3 in_initVertical) {
	this->r_initVertical = in_initVertical;
}



glm::vec3 Renderer::GetLensPos() {
	return r_lensPos;
}

void Renderer::SetLensPos(glm::vec3 in_lensPos) {
	this->r_lensPos = in_lensPos;
}