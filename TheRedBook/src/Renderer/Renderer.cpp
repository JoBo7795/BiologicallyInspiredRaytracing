#include "Renderer.h"


Camera* Renderer::camera = new Camera(glm::vec3(0.0f, 0.0f, 0.0f));
GLuint Renderer::renderMode = GL_FILL, Renderer::queueObjectID = 0;
std::vector<GameObject> Renderer::RenderQueue;
std::vector<GameObject> Renderer::HUDRenderQueue;
std::vector<GLuint> Renderer::RayTraceRenderQueue;
GLfloat Renderer::lastFrame = 0.0f, Renderer::deltaTime = 0.0f;

// TODO Need to refactor class

int Renderer::r_translationModelLoc;
int Renderer::r_rotationModelLoc;
int Renderer::r_scaleModelLoc;
int Renderer::r_customColorModelLoc;
int Renderer::r_textureSetModelLoc;
int Renderer::r_counter, Renderer::r_chunkRenderCounter;
int Renderer::r_renderQueueSize;
int Renderer::lightPosModelLoc;
int Renderer::viewModelLoc;
int Renderer::projectionModelLoc;
int Renderer::viewPosModelLoc;

float Renderer::lenseDistance = 0;//291.68f;
float Renderer::lenseSizeX = 1.0f;
float Renderer::lenseSizeY = 1.0f;
float Renderer::lenseSizeZ = 1.0f;
float Renderer::sceneBrightness = 1.0f;
glm::vec3 Renderer::lightPos = glm::vec3(1000.0f, -1000.0f, 0.0f);
bool Renderer::shadows_active = false;
bool Renderer::debugMode = false;
float Renderer::global_refract_index = 1.59f;
glm::vec2 Renderer::factXY;
bool Renderer::renderImageOnly = false;
float Renderer::height = 1920;
float Renderer::width = 1080;
int Renderer::renderDepth = 3;
glm::vec3 Renderer::initCamPos;
glm::vec3 Renderer::initlowerLeftCorner;
glm::vec3 Renderer::initHorizontal;
glm::vec3 Renderer::initVertical;
glm::vec3 Renderer::lensPos = glm::vec3(0.0);




void Renderer::InitChunkRenderer(Shader shaderProgram) {

}

void Renderer::RenderHUD(Shader shaderProgram) {

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
		glUniform3f(modelLoc, lensPos.x, lensPos.y, lensPos.z);

		modelLoc = glGetUniformLocation(shaderID, "lightPos");
		glUniform3f(modelLoc, lightPos.x, lightPos.y, lightPos.z);

		modelLoc = glGetUniformLocation(shaderID, "shadows_active");
		glUniform1i(modelLoc, shadows_active);

		modelLoc = glGetUniformLocation(shaderID, "debugMode");
		glUniform1i(modelLoc, debugMode);

		modelLoc = glGetUniformLocation(shaderID, "global_refract_index");
		glUniform1f(modelLoc, global_refract_index);

		modelLoc = glGetUniformLocation(shaderID, "factXY");
		glUniform2f(modelLoc, factXY.x, factXY.y);

		modelLoc = glGetUniformLocation(shaderID, "widthHeight");
		glUniform2f(modelLoc, width, height);

		modelLoc = glGetUniformLocation(shaderID, "renderImageOnly");
		glUniform1f(modelLoc, renderImageOnly);

		modelLoc = glGetUniformLocation(shaderID, "renderDepth");
		glUniform1i(modelLoc, renderDepth);
		
		glm::vec3 transView = glm::vec3(0.0f,0.0f,0.0f) * glm::normalize(camera->GetDirection()) * lenseDistance;
		auto lenseTransform = camera->GetViewMat();

		glm::vec4 posL = (camera->GetViewMat() * glm::vec4(0.0, 0.0, 0.0, 1.0));


		glm::vec3 normalizedDir = -posL;// 0.9 * 2 * -1;
		glm::vec3 lensOrig = camera->GetPosition() + normalizedDir;

		modelLoc = glGetUniformLocation(shaderID, "lenseTransform");
		glUniformMatrix4fv(modelLoc,1,false, &lenseTransform[0][0]);	

		modelLoc = glGetUniformLocation(shaderID, "initCamPos");
		glUniform3f(modelLoc, initCamPos.x, initCamPos.y, initCamPos.z);

		modelLoc = glGetUniformLocation(shaderID, "initlowerLeftCorner");
		glUniform3f(modelLoc, initlowerLeftCorner.x , initlowerLeftCorner.y, initlowerLeftCorner.z);

		modelLoc = glGetUniformLocation(shaderID, "initHorizontal");
		glUniform3f(modelLoc, initHorizontal.x, initHorizontal.y, initHorizontal.z);

		modelLoc = glGetUniformLocation(shaderID, "initVertical");
		glUniform3f(modelLoc, initVertical.x, initVertical.y, initVertical.z);
				
		modelLoc = glGetUniformLocation(shaderID, "sceneBrightness");
		glUniform1f(modelLoc, sceneBrightness);

		glDrawElements(GL_TRIANGLES, gameObject.GetIndiSize(), GL_UNSIGNED_INT, 0);
	}
}

void Renderer::RenderRayTrace(Shader shaderProgram) {

	shaderProgram.Use();

	int renderQueueSize = RayTraceRenderQueue.size();

	for (int i = 0; i < renderQueueSize; i++) {

		glBindBuffer(GL_SHADER_STORAGE_BUFFER,RayTraceRenderQueue[i]);

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