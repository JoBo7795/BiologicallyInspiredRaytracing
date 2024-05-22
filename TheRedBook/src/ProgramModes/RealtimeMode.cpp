#include "RealtimeMode.h"

glm::vec3 RealtimeMode::imagePlanePos, RealtimeMode::imagePlaneDir;
bool RealtimeMode::debugMode;
Camera* RealtimeMode::cam;
Camera RealtimeMode::debugCam;
int RealtimeMode::imageWidth;
int RealtimeMode::imageHeight;
LensData RealtimeMode::lens;



void RealtimeMode::InitPictureMode(float in_lenseDistance, bool in_debugMode, int in_rDepth) {
	cam = Renderer::GetCamera();
	Renderer::lenseDistance = in_lenseDistance;
	Renderer::renderDepth = in_rDepth;
	debugMode = in_debugMode;

	imageWidth = ProgramParams::windowWidth;
	imageHeight = ProgramParams::windowHeight;

}

void RealtimeMode::SetDebugParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir) {
	imagePlanePos = in_imagePlanePos;
	imagePlaneDir = in_imagePlaneDir;

	debugCam = Camera(imagePlanePos);

	debugCam.SetPosition(imagePlanePos);
	debugCam.SetDirection(imagePlaneDir);
	debugCam.InitRTCamera();
	debugCam.Update();


	Renderer::initHorizontal = debugCam.horizontal;
	Renderer::initVertical = debugCam.vertical;
	Renderer::initlowerLeftCorner = debugCam.lower_left_corner;
	Renderer::initCamPos = debugCam.GetPosition();
}

void RealtimeMode::SetCameraParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir) {
	imagePlanePos = in_imagePlanePos;
	imagePlaneDir = in_imagePlaneDir;

	cam->SetPosition(imagePlanePos);
	cam->SetDirection(imagePlaneDir);
	cam->InitRTCamera();
	cam->Update();

}


void RealtimeMode::Render() {

	ObjSelect currentGameObjectID;
	currentGameObjectID.objectID = -1;
	currentGameObjectID.objectType = -1;
	glm::vec3 currentMouseDirection;

	glm::vec3 pos;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_SCISSOR_TEST);
	
	Player player = Player(cam,cam->GetPosition());

	bool first = true;

	lens = LensData(player.GetPlayerPosition());

	while (!ProgramParams::window.WindowShouldClose()) {

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		
		InputControl::processInput(ProgramParams::window.GetWindowRef(), &player);

		if (ShaderManager::GetRegisteredShader(ProgramParams::shaderName, ProgramParams::tmpFragmentHUD))
			ProgramParams::FragmentHUD = ProgramParams::tmpFragmentHUD;

		auto posBefore = player.GetPlayerPosition();
		lens.lensOrigin = player.GetPlayerPosition();
		ProgramParams::f = lens.calcF();
		ProgramParams::dpt = lens.calcOneByF();

		pos = cam->GetPosition();

		ImGui::Begin("Player Position");
		ImGui::InputFloat("X: ", &pos.x, .005f, 1.0f);
		ImGui::InputFloat("Y: ", &pos.y, .05f, 1.0f);
		ImGui::InputFloat("Z: ", &pos.z, .05f, 1.0f);


		ImGui::Text(glm::to_string(currentMouseDirection).c_str());
		ImGui::End();

		ImGui::Begin("Refraction");		
		ImGui::Text(std::to_string(Renderer::global_refract_index).c_str());
		ImGui::InputFloat("Index(G): ", &Renderer::global_refract_index, .05f, 1.0f);
		ImGui::SliderFloat("Lense Distance ", &Renderer::lenseDistance, 0.0f, 3000.0f);
		ImGui::SliderFloat("Lense SizeX ", &Renderer::lenseSizeX, 0.0f, 3.0f);
		
		ImGui::SliderFloat("Lense SizeY ", &Renderer::lenseSizeY, 0.0f, 100.0f);
		ImGui::SliderFloat("Lense SizeZ ", &Renderer::lenseSizeZ, 0.0f, 100.0f);
		ImGui::SliderFloat("Retinal Image FormX", &Renderer::factXY.x, -100.0f, 100.0f);
		ImGui::SliderFloat("RaysPerPixel", &Renderer::factXY.y, 1.0f, 100.0f);

		ImGui::InputFloat("Index(G): ", &Renderer::factXY.x, .05f, 1.0f);
		ImGui::InputFloat("Index(G): ", &Renderer::factXY.y, .05f, 1.0f);

		ImGui::End();


		ProgramParams::g = 14.0;
		ProgramParams::calcImageDist();


		
		Renderer::lensPos = cam->GetPosition();
		cam->SetPosition(lens.lensOrigin + glm::normalize(lens.lensOriginLeft) * glm::vec3(ProgramParams::b));
		if (cam->updateData)
			cam->Update();
		if (first) {
			cam->SetDirection(glm::vec3(1.0f, 0.0f, 0.0f));
			
			cam->Update();
			Renderer::initHorizontal = cam->horizontal;
			Renderer::initVertical = cam->vertical;
			Renderer::initlowerLeftCorner = cam->lower_left_corner;
			Renderer::initCamPos = cam->GetPosition();
			first = false;

			// Activate to debug
			//cam->SetDirection(imagePlaneDir);

			
			cam->Update();

		}

		Renderer::RenderHUD(ProgramParams::FragmentHUD);

		currentMouseDirection = cam->GetPosition() + cam->GetDirection();

		ImGui::Begin("Framerate");
		ImGui::Text(std::to_string(ImGui::GetIO().Framerate).c_str());
		ImGui::End();

		ImGui::Begin("Selected GameObject");
		ImGui::Text("ID:");
		currentGameObjectID = Player::GetSelectedObject();
		if (currentGameObjectID.objectID == -1) {
			ImGui::Text("NONE");

		}
		else {
			ImGui::Text(std::to_string(currentGameObjectID.objectID).c_str());

			if (currentGameObjectID.objectType == SPHERE) {
				ImGui::Text("SPHERE");

			}
			if (currentGameObjectID.objectType == TRIANGLE) {
				ImGui::Text("TRIANGLE");

			}
		}

		ImGui::End();

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		ProgramParams::window.SwapChain();

	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

}