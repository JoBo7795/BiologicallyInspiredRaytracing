#include "PictureMode.h"

glm::vec3 PictureMode::imagePlanePos, PictureMode::imagePlaneDir;
bool PictureMode::debugMode;
Camera* PictureMode::cam;
Camera PictureMode::debugCam;

int PictureMode::imageWidth;
int PictureMode::imageHeight;
int PictureMode::numRenders;
int PictureMode::numLoops;
int PictureMode::stepX;
int PictureMode::stepY;
LensData PictureMode::lens;
const float PictureMode::infinity = std::numeric_limits<float>::infinity();


void PictureMode::InitPictureMode(float in_lenseDistance, bool in_renderImageOnly,bool in_debugMode, int in_rDepth) {
	cam = Renderer::GetCamera();
	Renderer::lenseDistance = in_lenseDistance;
	Renderer::renderImageOnly = in_renderImageOnly;
	Renderer::renderDepth = in_rDepth;
	debugMode = in_debugMode;

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


	Renderer::initHorizontal = debugCam.horizontal;
	Renderer::initVertical = debugCam.vertical;
	Renderer::initlowerLeftCorner = debugCam.lower_left_corner;
	Renderer::initCamPos = debugCam.GetPosition();
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

	Renderer::sceneBrightness = ProgramParams::sceneBrightness/100;

	if (debugMode) {

		 SetDebugParams(glm::vec3(0.0, 0.0, 0.0), glm::vec3(1.0, 0.0, 0.0));

		 lens.lensOrigin = debugCam.GetPosition();
		 Renderer::lensPos = lens.lensOrigin;
		 debugCam.SetPosition(lens.lensOrigin + glm::normalize(lens.lensOriginLeft) * glm::vec3(ProgramParams::b));
		 debugCam.Update();

		 SetDebugParams(debugCam.GetPosition(), glm::vec3(1.0, 0.0, 0.0));

	}
	else {

		lens.lensOrigin = debugCam.GetPosition();
		Renderer::lensPos = lens.lensOrigin;
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

			Renderer::debugMode = debugMode;
			Renderer::renderImageOnly = true;
			Renderer::renderDepth = ProgramParams::rDepth;
			Renderer::RenderHUD(ProgramParams::FragmentHUD);

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