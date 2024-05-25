#include "PictureMode.h"


PictureMode* PictureMode::instance = nullptr;

PictureMode* PictureMode::GetInstance() {

	if (instance == nullptr) {
		instance = new PictureMode();
	}

	return instance;
}

PictureMode::PictureMode(): infinity(std::numeric_limits<float>::infinity())
{
	rendererInstanceRef = Renderer::GetInstance();

	cam = Renderer::GetInstance()->GetCamera();
	Renderer::GetInstance()->SetLenseDistance(5.0f);
	Renderer::GetInstance()->SetRenderImageOnly(true);
	Renderer::GetInstance()->SetRenderDepth(4);
	debugMode = false;

	imageWidth = ProgramParams::windowWidth;
	imageHeight = ProgramParams::windowHeight;
	numRenders = 10000; // number of tiles per picture
	numLoops = sqrt(numRenders);
	stepX = imageWidth / numLoops;
	stepY = imageHeight / numLoops;
}


void PictureMode::SetDebugParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir) {
	imagePlanePos = in_imagePlanePos;
	imagePlaneDir = in_imagePlaneDir;

	debugCam = Camera(imagePlanePos);

	debugCam.SetPosition(imagePlanePos);
	debugCam.SetDirection(imagePlaneDir);
	debugCam.InitRTCamera();
	debugCam.Update();


	Renderer::GetInstance()->SetInitHorizontal(debugCam.horizontal);
	Renderer::GetInstance()->SetInitVertical (debugCam.vertical);
	Renderer::GetInstance()->SetInitlowerLeftCorner(debugCam.lower_left_corner);
	Renderer::GetInstance()->SetInitCamPos(debugCam.GetPosition());
}



void PictureMode::Render() {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_SCISSOR_TEST);


	std::chrono::time_point<std::chrono::high_resolution_clock> start, current, end;
	std::chrono::milliseconds duration,elapsed;
	int64_t hours, minutes, seconds,gHours, gMinutes, gSeconds, gDuration;

	bool first = true;

	lens = LensData(glm::vec3(0.0, 0.0, 0.0));
	ProgramParams::f = lens.calcF();
	ProgramParams::dpt = lens.calcOneByF();	

	auto sphereQueue = DrawObjectManager::GetSphereQueue();
	glm::vec3 rayPos = cam->GetPosition();
	
	rayPos = glm::vec3(0.0, 0.0, 0.0);

	double closestSoFar = infinity;
	double dist = 0.0;
	for (size_t i = 0; i < sphereQueue.size(); i++)
	{
		Sphere sphere = sphereQueue[i];
		glm::vec3 rayDir = glm::vec3(1.0, 0.0, 0.0);

		if (CollisionDetection::RaySphereIntersection(sphere.position, sphere.radius, rayPos, rayDir, 0.0, closestSoFar)) {
			dist = glm::length(glm::vec3(sphere.position - rayPos));
		}
	}
			

	if(dist > 1.0)
		ProgramParams::g = dist;

	ProgramParams::g = 5.0;

	ProgramParams::calcImageDist();

	Renderer::GetInstance()->SetSceneBrightness(ProgramParams::sceneBrightness/100);

	if (debugMode) {

		 SetDebugParams(glm::vec3(0.0, 0.0, 0.0), glm::vec3(1.0, 0.0, 0.0));

		 lens.lensOrigin = debugCam.GetPosition();
		 Renderer::GetInstance()->SetLensPos(lens.lensOrigin);
		 debugCam.SetPosition(lens.lensOrigin + glm::normalize(lens.lensOriginLeft) * glm::vec3(ProgramParams::b));
		 debugCam.Update();

		 SetDebugParams(debugCam.GetPosition(), glm::vec3(1.0, 0.0, 0.0));

	}
	else {

		lens.lensOrigin = debugCam.GetPosition();
		Renderer::GetInstance()->SetLensPos(lens.lensOrigin);
		cam->SetPosition(lens.lensOrigin + glm::normalize(lens.lensOriginLeft) * glm::vec3(ProgramParams::b));
		cam->Update();	

	}
	int nextProg = 10;
	start = std::chrono::high_resolution_clock::now();
	std::cout << std::flush << "rendering progress: " << std::setprecision(2) << static_cast<float>(0) * 100 << "%" << std::endl;
	for (int i = 0; i < numLoops; i++) {
		for (int j = 0; j < numLoops; j++) {
			
			if (ShaderManager::GetRegisteredShader(ProgramParams::shaderName, ProgramParams::tmpFragmentHUD))
				ProgramParams::FragmentHUD = ProgramParams::tmpFragmentHUD;			

			if (cam->updateData)
				cam->Update();
					

			glScissor(i * stepX, j * stepY, stepX, stepY);

			Renderer::GetInstance()->SetDebugMode(debugMode);
			Renderer::GetInstance()->SetRenderImageOnly (true);
			Renderer::GetInstance()->SetRenderDepth(ProgramParams::rDepth);
			Renderer::GetInstance()->RenderRayTraceShader(ProgramParams::FragmentHUD);

			double progress = static_cast<double>(i * numLoops + j) / numRenders;


			if (progress > 0.01 && first) {
				end = std::chrono::high_resolution_clock::now();
				duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
				gDuration = duration.count()* 100;

				gHours = gDuration / 3600000;
				gDuration -= gHours * 3600000;
				gMinutes = gDuration / 60000;
				gDuration -= gMinutes * 60000;
				gSeconds = gDuration / 1000;

				first = false;

				std::cout << std::flush << "rendering progress: " << std::setprecision(2) << static_cast<float>(progress) * 100 << "%" << std::endl << std::endl;
				std::cout << "Guessed time " << gHours << ":" << gMinutes << ":" << gSeconds << "." << std::endl;
				
			}
			else if (!first){				

				current = std::chrono::high_resolution_clock::now();
				elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(current - start);

				hours = elapsed.count() / 3600000;
				elapsed -= std::chrono::hours(hours);
				minutes = elapsed.count() / 60000;
				elapsed -= std::chrono::minutes(minutes);
				seconds = elapsed.count() / 1000;

				float calc = std::round((std::fmod(static_cast<float>(progress * 100) , nextProg)));
				if ( calc == 0) {
					std::cout << std::flush << "rendering progress: " << std::setprecision(2) << static_cast<float>(progress) * 100 << "%" << std::endl;
					std::cout << "Elapsed time: " << hours << ":" << minutes << ":" << seconds << "." << std::endl << std::endl;
					nextProg += 10;
				}

			}


		}
	}
	
}

glm::vec3 PictureMode::GetImagePlanePos()  {
	return imagePlanePos;
}


void PictureMode::SetImagePlanePos(glm::vec3& newPos) {
	imagePlanePos = newPos;
}

glm::vec3 PictureMode::GetImagePlaneDir() {
	return imagePlaneDir;
}


void PictureMode::SetImagePlaneDir( glm::vec3& newDir) {
	imagePlaneDir = newDir;
}


bool PictureMode::IsDebugMode() {
	return debugMode;
}


void PictureMode::SetDebugMode(bool newDebugMode) {
	debugMode = newDebugMode;
}


Camera* PictureMode::GetCamera()  {
	return cam;
}


void PictureMode::SetCamera(Camera* newCam) {
	cam = newCam;
}


Camera PictureMode::GetDebugCamera()  {
	return debugCam;
}


void PictureMode::SetDebugCamera(Camera newDebugCam) {
	debugCam = newDebugCam;
}


LensData PictureMode::GetLens() {
	return lens;
}


void PictureMode::SetLens(LensData& newLens) {
	lens = newLens;
}


Renderer* PictureMode::GetRendererInstanceRef() {
	return rendererInstanceRef;
}


void PictureMode::SetRendererInstanceRef(Renderer* newRenderer) {
	rendererInstanceRef = newRenderer;
}


float PictureMode::RendererGetLensDistance() {
	return Renderer::GetInstance()->GetLenseDistance();
}

void PictureMode::RendererSetLensDistance(float lensDistance) {
	Renderer::GetInstance()->SetLenseDistance(lensDistance);
}

bool PictureMode::RendererGetRenderImageOnly() {

	return Renderer::GetInstance()->GetRenderImageOnly();
}

void PictureMode::RendererSetRenderImageOnly(bool renderImageOnly) {
	Renderer::GetInstance()->SetRenderImageOnly(renderImageOnly);
}

int PictureMode::RendererGetRenderDepth() {
	return Renderer::GetInstance()->GetRenderDepth();

}
void PictureMode::RendererSetRenderDepth(int renderDepth) {

	Renderer::GetInstance()->SetLenseDistance(renderDepth);
}