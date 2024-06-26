#include "SceneDescription.h"

#define RENDER_IMAGE_ONLY true
#define CPU_RENDER false
#define DEBUG_MODE true
void Scene::SceneDescription() {

	ProgramParams::picOnly = RENDER_IMAGE_ONLY;
	ProgramParams::CPURender = CPU_RENDER;
	ProgramParams::debugMode = DEBUG_MODE;

	float lenseScale = ProgramParams::lenseScale;
	float lenseDistance = ProgramParams::lenseDistance;
	float g = ProgramParams::g;

	float posX = ProgramParams::posX, posY =ProgramParams::posY, posZ = ProgramParams::posZ;

	lenseDistance = 2.136;
	ProgramParams::f = lenseDistance;


	glm::vec3 imagePlanePos(-lenseDistance, 0, 0);
	glm::vec3 imagePlaneDir(normalize(glm::vec3(1, 0, 0)));


	// Sphere Objects //
	 Sphere sphere;	 
	 
	 sphere.radius = .2;
	 sphere.position = glm::vec3(5.2, 1.0f, 2.0f);
	 sphere.color = glm::vec4(.3f, .5f, .7f, 1.0f);
	 sphere.albedo = glm::vec4(1.,0,0, 1.0f);
	 sphere.fuzz = 1.0f;
	 sphere.mat_ptr = PLAIN;
	 DrawObjectManager::RegisterRTSphere(sphere);
	 
	 sphere.radius = .2;
	 sphere.position = glm::vec3(105.0, -.3f, -.3f);	 
	 sphere.color = glm::vec4(.7f, .5f, .3f, 1.0f);
	 sphere.albedo = glm::vec4(0., .0, 1.0, 1.0f);
	 sphere.fuzz = 1.0f;
	 sphere.mat_ptr = PLAIN;
	 DrawObjectManager::RegisterRTSphere(sphere);

	 sphere.radius = .2;
	 sphere.position = glm::vec3(10.0, -10.75f, 0.0f);	 
	 sphere.color = glm::vec4(.7f, .5f, .3f, 1.0f);
	 sphere.albedo = glm::vec4(0., 0.0, 1.0, 1.0f);
	 sphere.fuzz = 1.0f;
	 sphere.mat_ptr = PLAIN;
	 DrawObjectManager::RegisterRTSphere(sphere);
	 
	 sphere.radius = .05f;
	 sphere.position = glm::vec3(-300.0f, -100.5f, 0.0f);
	 sphere.color = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
	 sphere.albedo = glm::vec4(0.0, 0.0, 1.0, 1.0f);
	 sphere.fuzz = 1.0f;
	 sphere.mat_ptr = DIELECTRIC;
	 DrawObjectManager::RegisterRTSphere(sphere);
	 

	 // Tropfen

	 sphere.radius = .1f;
	 sphere.position = glm::vec3(5.5f, 10.0f, 0.4f);
	 sphere.color = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
	 sphere.albedo = glm::vec4(1.0, 0.0, 0.0, 1.0f);
	 sphere.fuzz = 1.0f;
	 sphere.mat_ptr = DIELECTRIC;
	 DrawObjectManager::RegisterRTSphere(sphere);
	 
	 sphere.radius = .01f;
	 sphere.position = glm::vec3(5.5f, -10.0f, -0.4f);
	 sphere.color = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
	 sphere.albedo = glm::vec4(0.0, 1.0, 0.0, 1.0f);
	 sphere.fuzz = 1.0f;
	 sphere.mat_ptr = DIELECTRIC;
	 DrawObjectManager::RegisterRTSphere(sphere);


	 // Triangle-Mesh Objects //

	 char buffer[FILENAME_MAX];
	
	 std::string(_getcwd(buffer, FILENAME_MAX)) + std::string(OBJECTPATH);
	 ShapeData shape = DataReader::readPLY(std::string(_getcwd(buffer, FILENAME_MAX)) + std::string(OBJECTPATH) + std::string("paintedCube.ply"));

	 shape.matPointer = PLAIN;
	 shape.albedo = glm::vec4(0.0,0.0,0.0,1.0);
	 DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(shape));

	 GameObject go = GameObject(DrawObjectManager::GetCurrentBound().id, 0);
	 go.SetTranslation(glm::translate(glm::vec3(10.5, 25.0, 0.0f))* glm::scale(glm::vec3(.15f, .15f, .15f)));
	   
	  GameObjectManager::AppendGameObject(go);



	// Camera Options //
	auto cam = Renderer::GetInstance()->GetCamera();

	cam->SetPosition(imagePlanePos);
	cam->SetDirection(imagePlaneDir);
	cam->Update();

	RealtimeMode* realTimeMode = RealtimeMode::GetInstance();
	PictureMode* pictureMode = PictureMode::GetInstance();

	// Debug Position Parameters //
	if (DEBUG_MODE) {
		
		glm::vec3 debugImagePlanePos(0.0, posY, posX);
		glm::vec3 debugImagePlaneDir(glm::normalize(glm::vec3(1.0,0.0,0.0)));

		pictureMode->SetDebugParams(debugImagePlanePos, debugImagePlaneDir);
		realTimeMode->SetDebugParams(debugImagePlanePos, debugImagePlaneDir);
		CPUMode::SetDebugParams(debugImagePlanePos, debugImagePlaneDir);

		cam->SetPosition(glm::vec3(0.0, 0.0, -40.0));
		cam->SetDirection(normalize(glm::vec3(0.0, 0.0, 0.0) - cam->GetPosition()));
		cam->Update();
	}



	if (RENDER_IMAGE_ONLY) {
		if (CPU_RENDER) {
			CPUMode::InitCPUMode(imagePlanePos, imagePlaneDir);			
		}
		else {
			pictureMode->SetDebugMode(DEBUG_MODE);
			pictureMode->RendererSetLensDistance(lenseDistance);
			pictureMode->RendererSetRenderImageOnly(true);

		}
		
	}
	else {

		realTimeMode->SetDebugMode(DEBUG_MODE);
		realTimeMode->RendererSetLensDistance(lenseDistance);
		realTimeMode->RendererSetRenderImageOnly(true);
		realTimeMode->SetCameraParams(imagePlanePos, imagePlaneDir);
		
	}
}