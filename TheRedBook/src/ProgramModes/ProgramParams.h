#pragma once
#include <string>
#include "../Definitions.h"
#include "../Shader/ShaderManager.h"
#include "../Window/CreateWindow.h"

#include "../../Dependencies/imgui/imgui.h"
#include "../../Dependencies/imgui/imgui_impl_glfw.h"
#include "../../Dependencies/imgui/imgui_impl_opengl3.h"

class ProgramParams{
public:
	static float windowWidth, windowHeight;
	static float posX, posY, posZ;
	static double g,b,dpt,f;

	static std::string path;
	static std::string paramPaths;
	static bool picOnly;
	static bool CPURender;
	static float lenseScale;
	static float lenseDistance;
	static float sceneBrightness;
	static int rDepth;
	static std::string params;
	static bool debugMode;
	static Shader shaderProgram1;
	static Shader shaderProgram2;
	static std::string shaderName;
	static Shader FragmentHUD, tmpFragmentHUD;
	static Shader currentHUDFragment;
	static Window window;

	// object and image wide


	static void ComputeProgramArgs(int argc, char* argv[]) {

		if (argc > 1) {
			for (int i = 1; i < argc; i++) {
				std::string arg = argv[i];
				size_t pos = arg.find("=");
				if (pos != std::string::npos) {
					std::string key = arg.substr(0, pos);
					std::string value_str = arg.substr(pos + 1);

					if (key == "posX") {
						posX = std::stof(value_str);
						params += ("_posX" + value_str);
						std::cout << "posX: " << posX << std::endl;
					}
					else if (key == "posY") {
						posY = std::stof(value_str);
						params += ("_posY" + value_str);
						std::cout << "posY: " << posY << std::endl;
					}
					else if (key == "posZ") {
						posZ = std::stof(value_str);
						params += ("_posZ" + value_str);
						std::cout << "posZ: " << posZ << std::endl;
					}
					else if (key == "sbright") {
						sceneBrightness = std::stof(value_str);
						params += ("_sceneBrigthness" + value_str);
						std::cout << "sceneBrigthness: " << posZ << std::endl;
					}
					else if (key == "lScale") {
						lenseScale = std::stof(value_str);
						params += ("_lenseScale" + value_str);
						std::cout << "lenseScale: " << lenseScale << std::endl;
					}
					else if (key == "lDist") {
						lenseDistance = std::stof(value_str);
						params += ("_lenseDistance" + value_str);
						std::cout << "lenseDistance: " << lenseDistance << std::endl;
					}
					else if (key == "picWidth") {
						windowWidth = std::stof(value_str);
						params += ("_picWidth" + value_str);
						std::cout << "picWidth: " << windowWidth << std::endl;
					}
					else if (key == "picHeight") {
						windowHeight = std::stof(value_str);
						params += ("_picHeight" + value_str);
						std::cout << "picHeight: " << windowHeight << std::endl;
					}
					else if (key == "rDepth") {
						rDepth = std::stoi(value_str);
						params += ("_rDepth" + value_str);
						std::cout << "rDepth: " << rDepth << std::endl;
					}
					else if (key == "path") {
						paramPaths += value_str;
						std::cout << "path: " << path << std::endl;
					}
					else if (key == "bWidth") {
						g = std::stof(value_str);
						params += ("_bWidth" + value_str);
						std::cout << "bWidth: " << g << std::endl;
					}
					else if (key == "debug") {

						if (value_str == "true" || value_str == "1") {
							debugMode = true;
						}
						else
						{
							debugMode = false;
						}
						paramPaths += value_str;
						std::cout << "Debug mode with rays activated" << std::endl;
					}
					else if (key == "pic") {
						if (std::stoi(value_str)) {
							picOnly = true;
							std::cout << "Only a picture will be rendered into folder: " << path << std::endl;
						}
						else {
							std::cout << "Rendering realtime" << std::endl;
						}
					}
					else {
						std::cout << "Unbekannter Schlüssel: " << key << std::endl;
					}
				}
				else {
					std::cout << "Kein gültiges Argument angegeben" << std::endl;
				}
			}
		}
		else {
			std::cout << "Kein Argument angegeben" << std::endl;
		}
	}

	static void ProgramInit() {

		window = CreateWindow::CreateWindow(windowWidth, windowHeight);

		Renderer::GetInstance()->SetHeight(windowHeight);
		Renderer::GetInstance()->SetWidth(windowWidth);
		Renderer::GetInstance()->SetInitCamPos(glm::vec3(-lenseDistance, posY, posZ));

		window.SetCallback(BUFFER_SIZE);
		window.SetCallback(MOUSE_INPUT);
		window.SetCallback(SCROLL_INPUT);

		glewInit();

		std::vector<ShaderInfo> shaderInfo1 = {
			{GL_VERTEX_SHADER, VERTEX_SHADER_PATH, NULL },
			{GL_FRAGMENT_SHADER, FRAGMENT_SHADER_PATH, NULL}
		};

		std::vector<ShaderInfo> shaderInfo2 = {
			{GL_VERTEX_SHADER, VERTEX_SHADER_PATH_RAYTRACE, NULL },
			{GL_FRAGMENT_SHADER, FRAGMENT_SHADER_PATH_RAYTRACE, NULL}
		};

		shaderProgram1 = CreateShader::GenerateShader(shaderInfo1);
		shaderProgram2 = CreateShader::GenerateShader(shaderInfo2);
		shaderName = std::string(HUD_VERTEX_SHADER);
		ShaderManager::RegisterShader(shaderName, shaderProgram1);
		shaderName = std::string(HUD_FRAGMENT_SHADER);
		ShaderManager::RegisterShader(shaderName, shaderProgram2);

		ShaderManager::GetRegisteredShader(shaderName, FragmentHUD);
		currentHUDFragment = FragmentHUD;

		IMGUI_CHECKVERSION();
		ImGui::CreateContext();
		auto io = ImGui::GetIO();
		ImGui::StyleColorsDark();
		ImGui_ImplGlfw_InitForOpenGL(window.GetWindowRef(), true);
		ImGui_ImplOpenGL3_Init("#version 460");

		glClearColor(1.0f, 0.0f, 0.0f, 1.0f);


		// Initial camera setup necessary
		auto cam = Renderer::GetInstance()->GetCamera();
		cam->SetPosition(glm::vec3(1.0f, 0.0f, 0.0f));
		cam->SetDirection(normalize(glm::vec3(0.0f, 0.0f, 0.0f) - cam->GetPosition()));
		cam->Update();
		cam->InitRTCamera();

	}

	static void calcImageDist() {
		
		ProgramParams::b = 1.0 / ((ProgramParams::dpt) - (1.0 / ProgramParams::g));

	}
};


