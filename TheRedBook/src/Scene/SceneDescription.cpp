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

	// f = 2.52
	lenseDistance = 2.136;
	//lenseDistance = 0.847457627118644;
	ProgramParams::f = lenseDistance;
	g = 10.0f;

	//glm::vec3 imagePlanePos(-lenseDistance, posY, posZ);
	glm::vec3 imagePlanePos(-lenseDistance, 0, 0);
	//glm::vec3 imagePlaneDir(normalize(glm::vec3(-lenseDistance, posY, posZ)) - normalize(glm::vec3(g, 0.0f, 0.0f)));
	//glm::vec3 imagePlaneDir(normalize(glm::vec3((-lenseDistance), posY, posZ) - glm::vec3(0.0f, 0.0f, 0.0f)));
	//glm::vec3 imagePlaneDir(normalize(glm::vec3(1, posY, posZ)));
	glm::vec3 imagePlaneDir(normalize(glm::vec3(1, 0, 0)));


	// Sphere Objects //
	 Sphere sphere;
	 
	 
	 sphere.radius = .2;
	 sphere.position = glm::vec3(5.2, 1.0f, 2.0f);
	 //sphere.position = glm::vec3(g, 10.5f, 0.0f);
	 sphere.color = glm::vec4(.3f, .5f, .7f, 1.0f);
	 sphere.albedo = glm::vec4(1.,0,0, 1.0f);
	 sphere.fuzz = 1.0f;
	 sphere.mat_ptr = PLAIN;
	 DrawObjectManager::RegisterRTSphere(sphere);
	 
	 sphere.radius = .2;
	 sphere.position = glm::vec3(105.0, -.3f, -.3f);
	 //sphere.position = glm::vec3(g, 10.5f, 0.0f);
	 sphere.color = glm::vec4(.7f, .5f, .3f, 1.0f);
	 sphere.albedo = glm::vec4(0., .0, 1.0, 1.0f);
	 sphere.fuzz = 1.0f;
	 sphere.mat_ptr = PLAIN;
	 DrawObjectManager::RegisterRTSphere(sphere);

	 sphere.radius = .2;
	 sphere.position = glm::vec3(10.0, -10.75f, 0.0f);
	 //sphere.position = glm::vec3(g, 10.5f, 0.0f);
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
	 
	  //phere.radius = .05f;
	  //phere.position = glm::vec3(-(2.23642172523962+.05), -100.25f, 0.0f);
	  //phere.color = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
	  //phere.albedo = glm::vec4(0.0, 1.0, 0.0, 1.0f);
	  //phere.fuzz = 1.0f;
	  //phere.mat_ptr = PLAIN;
	  //rawObjectManager::RegisterRTSphere(sphere);
	  //

	  //
	  //sphere.radius = .05f;
	  //sphere.position = glm::vec3(0.0, -100.25f, 0.0f);
	  //sphere.color = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
	  //sphere.albedo = glm::vec4(1.0, 0.0, 0.0, 1.0f);
	  //sphere.fuzz = 1.0f;
	  //sphere.mat_ptr = PLAIN;
	  //DrawObjectManager::RegisterRTSphere(sphere);







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



	   
	   // sphere.radius = .025f;
	   // //sphere.position = glm::vec3(-lenseDistance, 1.5f, 0.0f);
	   // sphere.position = glm::vec3((3.0f), .1f, -.15f);
	   // sphere.color = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
	   // sphere.albedo = glm::vec4(0.0, 0.0, 1.0, 1.0f);
	   // sphere.fuzz = 1.0f;
	   // sphere.mat_ptr = DIELECTRIC;
	   // DrawObjectManager::RegisterRTSphere(sphere);
	   // 
	   // 
	   // sphere.radius = .025f;
	   // sphere.position = glm::vec3(3.0f, -.025f, .0f);
	   // sphere.color = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
	   // sphere.albedo = glm::vec4(1.0, 0.0, 0.0, 1.0f);
	   // sphere.fuzz = 1.0f;
	   // sphere.mat_ptr = DIELECTRIC;
	   // DrawObjectManager::RegisterRTSphere(sphere);
	   // 
	   // sphere.radius = .025f;
	   // sphere.position = glm::vec3(3.0f, .1f, -.25f);
	   // sphere.color = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
	   // sphere.albedo = glm::vec4(0.0, 1.0, 0.0, 1.0f);
	   // sphere.fuzz = 1.0f;
	   // sphere.mat_ptr = DIELECTRIC;
	   // DrawObjectManager::RegisterRTSphere(sphere);
	   // 
	   // sphere.radius = .02f;
	   // //sphere.position = glm::vec3(-lenseDistance, 1.5f, 0.0f);
	   // sphere.position = glm::vec3((3.0f), 0.1f, .30f);
	   // sphere.color = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
	   // sphere.albedo = glm::vec4(0.0, 0.0, 1.0, 1.0f);
	   // sphere.fuzz = 1.0f;
	   // sphere.mat_ptr = DIELECTRIC;
	   // DrawObjectManager::RegisterRTSphere(sphere);
	   // 
	   // 
	   // sphere.radius = .025f;
	   // sphere.position = glm::vec3(3.0f, -.2f, .20f);
	   // sphere.color = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
	   // sphere.albedo = glm::vec4(1.0, 0.0, 0.0, 1.0f);
	   // sphere.fuzz = 1.0f;
	   // sphere.mat_ptr = DIELECTRIC;
	   // DrawObjectManager::RegisterRTSphere(sphere);
	   // 
	   // sphere.radius = .025f;
	   // sphere.position = glm::vec3(3.0f, -.0f, -0.2f);
	   // sphere.color = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
	   // sphere.albedo = glm::vec4(0.0, 1.0, 0.0, 1.0f);
	   // sphere.fuzz = 1.0f;
	   // sphere.mat_ptr = DIELECTRIC;
	   // DrawObjectManager::RegisterRTSphere(sphere);
	   // 
	   // sphere.radius = .02f;
	   // //sphere.position = glm::vec3(-lenseDistance, 1.5f, 0.0f);
	   // sphere.position = glm::vec3((3.0f), 0.2f, .20f);
	   // sphere.color = glm::vec4(0.0f, 1.0f, 0.0f, 1.0f);
	   // sphere.albedo = glm::vec4(0.0, 0.0, 1.0, 1.0f);
	   // sphere.fuzz = 1.0f;
	   // sphere.mat_ptr = DIELECTRIC;
	   // DrawObjectManager::RegisterRTSphere(sphere);









	// Triangle-Mesh Objects //
	 //ShapeData shape = DataReader::readPLY(OBJECTPATH + std::string("particle.ply"));
	 //ShapeData shape = DataReader::readPLY(OBJECTPATH + std::string("drop.ply"));
	 
	   ShapeData shape = DataReader::readPLY(OBJECTPATH + std::string("paintedCube.ply"));
	   //ShapeData shape = DataReader::readPLY(OBJECTPATH + std::string("lens.ply"));

	 //ShapeData shape = DataReader::readPLY(OBJECTPATH + std::string("affekopp.ply"));
	 //ShapeData shape = DataReader::readPLY(OBJECTPATH + std::string("weihnachtenRotated.ply"));
	 shape.matPointer = PLAIN;
	 shape.albedo = glm::vec4(0.0,0.0,0.0,1.0);
	 DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(shape));

//	//DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(DataReader::readPLY(OBJECTPATH+std::string("NeueLinseBlender.ply"))));
//	//DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(DataReader::readPLY(OBJECTPATH+std::string("concave_lens_blender.ply"))));
//	//DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(DataReader::readPLY(OBJECTPATH+std::string("testkugelBlender.ply"))));
//	//DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(DataReader::readPLY(OBJECTPATH+std::string("EbenePlatteBlender.ply"))));
//	//DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(DataReader::readPLY(OBJECTPATH+std::string("NeueLinseGroßBlender.ply"))));
//	//DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(DataReader::readPLY(OBJECTPATH+std::string("lens.ply"))));
//	//DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(DataReader::readPLY(OBJECTPATH+std::string("refractionSceneA.ply"))));
//	//DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(DataReader::readPLY(OBJECTPATH+std::string("plane.ply"))));
//	//DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(DataReader::readPLY(OBJECTPATH+std::string("Platteplatte.ply"))));
//	//DrawObjectManager::UseDrawObject(DrawObjectManager::RegisterRTDrawObject(DataReader::readPLY(OBJECTPATH+std::string("planeBig.ply"))));




	 GameObject go = GameObject(DrawObjectManager::GetCurrentBound().id, 0);




//	//go.SetTranslation(glm::scale(glm::vec3(0.38977829817f, 0.38977829817f, 0.38977829817f)));
//	//go.SetTranslation(glm::scale(glm::vec3(10.0f, 1.0f, 1.0f)));
	//posY = -8;

	//go.SetTranslation(glm::translate(glm::vec3(10.3, 0.0, 0.0f))* glm::scale(glm::vec3(.15f, .15f, .15f)));// *glm::rotate(45.0f, glm::vec3(1.0f, 1.0f, 0.0f)));
	 go.SetTranslation(glm::translate(glm::vec3(10.5, 25.0, 0.0f))* glm::scale(glm::vec3(.15f, .15f, .15f)));// *glm::rotate(180.0f, glm::vec3(0.0f, 1.0f, 0.0f))




//	//go.SetTranslation(glm::rotate(180.0f, glm::vec3(0.0f, 1.0f, 0.0f)) );
//	//go.SetTranslation(glm::translate(glm::vec3(10.0f,2.0f,0.0f)) *glm::scale(glm::mat4(1.0),glm::vec3(5.0)) );
//	go.SetTranslation(glm::translate(glm::vec3(3.6f,.0f,0.0f))* glm::scale(glm::mat4(1.0), glm::vec3(.05)));
	

	   
	  GameObjectManager::AppendGameObject(go);



	// Camera Options //
	auto cam = Renderer::GetCamera();
	//cam->SetPosition(glm::vec3(-60.0f, 0.0f, -250.0f));	
	//cam->SetDirection(normalize(glm::vec3(30.0f, 0.0f, 0.0f) - cam->GetPosition()));

	cam->SetPosition(imagePlanePos);
	cam->SetDirection(imagePlaneDir);
	cam->Update();



	// Debug Position Parameters //
	if (DEBUG_MODE) {
		//glm::vec3 debugImagePlanePos(-10.0, posY, posX);
		glm::vec3 debugImagePlanePos(0.0, posY, posX);
		//glm::vec3 debugImagePlaneDir(normalize(glm::vec3(-lenseDistance, posY, posZ)) - normalize(glm::vec3(g, 0.0f, 0.0f)));
		glm::vec3 debugImagePlaneDir(glm::normalize(glm::vec3(1.0,0.0,0.0)));
		//glm::vec3 debugImagePlaneDir(normalize(glm::vec3(lenseDistance, posY, posZ)));
		//glm::vec3 debugImagePlaneDir(normalize(glm::vec3(-lenseDistance, posY, posZ) - glm::vec3(0.0f, 0.0f, 0.0f)));


		PictureMode::SetDebugParams(debugImagePlanePos, debugImagePlaneDir);
		RealtimeMode::SetDebugParams(debugImagePlanePos, debugImagePlaneDir);
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
			PictureMode::InitPictureMode(lenseDistance, true, DEBUG_MODE);
		}
		
	}
	else {
		RealtimeMode::InitPictureMode(lenseDistance, true, DEBUG_MODE);
		RealtimeMode::SetCameraParams(imagePlanePos, imagePlaneDir);
		
	}
}