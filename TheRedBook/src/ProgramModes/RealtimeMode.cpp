#include "RealtimeMode.h"

RealtimeMode* RealtimeMode::instance = nullptr;

RealtimeMode* RealtimeMode::GetInstance() {

	if (instance == nullptr) {
		instance = new RealtimeMode();
	}

	return instance;
}

RealtimeMode::RealtimeMode() {
	rendererInstanceRef = Renderer::GetInstance();
	cam = rendererInstanceRef->GetCamera();
	rendererInstanceRef->SetLenseDistance(5.0f);
	rendererInstanceRef->SetRenderDepth(4);
	debugMode = false;

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


	rendererInstanceRef->SetInitHorizontal(debugCam.horizontal);
	rendererInstanceRef->SetInitVertical(debugCam.vertical);
	rendererInstanceRef->SetInitlowerLeftCorner(debugCam.lower_left_corner);
	rendererInstanceRef->SetInitCamPos(debugCam.GetPosition());
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
		float refractIndex  = rendererInstanceRef->GetGlobalRefractIndex();
		float lensDistance = rendererInstanceRef->GetLenseDistance();
		float lensSizeX	= rendererInstanceRef->GetLenseSizeX();
		float lensSizeY	= rendererInstanceRef->GetLenseSizeY();
		float lensSizeZ	= rendererInstanceRef->GetLenseSizeZ();
		glm::vec2 factXY = rendererInstanceRef->GetFactXY();

		ImGui::Text(std::to_string(refractIndex).c_str());

		ImGui::InputFloat("Index(G): ", &refractIndex, .05f, 1.0f);
		ImGui::SliderFloat("Lense Distance ", &lensDistance, 0.0f, 3000.0f);

		ImGui::SliderFloat("Lense SizeX ", &lensSizeX, 0.0f, 3.0f);		
		ImGui::SliderFloat("Lense SizeY ", &lensSizeY, 0.0f, 100.0f);
		ImGui::SliderFloat("Lense SizeZ ", &lensSizeZ, 0.0f, 100.0f);

		ImGui::SliderFloat("Retinal Image FormX", &factXY.x, -100.0f, 100.0f);
		ImGui::SliderFloat("RaysPerPixel", &factXY.y, 1.0f, 100.0f);

		ImGui::InputFloat("Index(G): ", &factXY.x, .05f, 1.0f);
		ImGui::InputFloat("Index(G): ", &factXY.y, .05f, 1.0f);

		ImGui::End();

		rendererInstanceRef->SetGlobalRefractIndex(refractIndex);
		rendererInstanceRef->SetLenseDistance(lensDistance);
		rendererInstanceRef->SetLenseSizeX(lensSizeX);
		rendererInstanceRef->SetLenseSizeY(lensSizeY);
		rendererInstanceRef->SetLenseSizeZ(lensSizeZ);
		rendererInstanceRef->SetFactXY(factXY);



		ProgramParams::g = 14.0;
		ProgramParams::calcImageDist();


		
		rendererInstanceRef->SetLensPos(cam->GetPosition());
		cam->SetPosition(lens.lensOrigin + glm::normalize(lens.lensOriginLeft) * glm::vec3(ProgramParams::b));
		if (cam->updateData)
			cam->Update();
		if (first) {
			cam->SetDirection(glm::vec3(1.0f, 0.0f, 0.0f));
			
			cam->Update();
			rendererInstanceRef->SetInitHorizontal(cam->horizontal);
			rendererInstanceRef->SetInitVertical(cam->vertical);
			rendererInstanceRef->SetInitlowerLeftCorner(cam->lower_left_corner);
			rendererInstanceRef->SetInitCamPos(cam->GetPosition());
			first = false;

			// Activate to debug
			//cam->SetDirection(imagePlaneDir);

			
			cam->Update();

		}

		rendererInstanceRef->RenderRayTraceShader(ProgramParams::FragmentHUD);

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

glm::vec3 RealtimeMode::GetImagePlanePos() {
	return imagePlanePos;
}


void RealtimeMode::SetImagePlanePos(glm::vec3& newPos) {
	imagePlanePos = newPos;
}

glm::vec3 RealtimeMode::GetImagePlaneDir() {
	return imagePlaneDir;
}


void RealtimeMode::SetImagePlaneDir(glm::vec3& newDir) {
	imagePlaneDir = newDir;
}


bool RealtimeMode::IsDebugMode() {
	return debugMode;
}


void RealtimeMode::SetDebugMode(bool newDebugMode) {
	debugMode = newDebugMode;
}


Camera* RealtimeMode::GetCamera() {
	return cam;
}


void RealtimeMode::SetCamera(Camera* newCam) {
	cam = newCam;
}


Camera RealtimeMode::GetDebugCamera() {
	return debugCam;
}


void RealtimeMode::SetDebugCamera(Camera newDebugCam) {
	debugCam = newDebugCam;
}


LensData RealtimeMode::GetLens() {
	return lens;
}


void RealtimeMode::SetLens(LensData& newLens) {
	lens = newLens;
}


Renderer* RealtimeMode::GetRendererInstanceRef() {
	return rendererInstanceRef;
}


void RealtimeMode::SetRendererInstanceRef(Renderer* newRenderer) {
	rendererInstanceRef = newRenderer;
}

float RealtimeMode::RendererGetLensDistance() {
	return Renderer::GetInstance()->GetLenseDistance();
}

void RealtimeMode::RendererSetLensDistance(float lensDistance) {
	Renderer::GetInstance()->SetLenseDistance(lensDistance);
}

bool RealtimeMode::RendererGetRenderImageOnly() {

	return Renderer::GetInstance()->GetRenderImageOnly();
}

void RealtimeMode::RendererSetRenderImageOnly(bool renderImageOnly) {
	Renderer::GetInstance()->SetRenderImageOnly(renderImageOnly);
}

int RealtimeMode::RendererGetRenderDepth() {
	return Renderer::GetInstance()->GetRenderDepth();

}
void RealtimeMode::RendererSetRenderDepth(int renderDepth) {
	Renderer::GetInstance()->SetLenseDistance(renderDepth);
}