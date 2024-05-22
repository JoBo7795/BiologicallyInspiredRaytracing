#include "GL/glew.h"
#include "GLFW/glfw3.h"

#include <iostream>
#include <vector>

#include <time.h>  
#include <array>

#include <string>
#include <chrono>
#include <iomanip>


#include "Objects/GameObjectManager.h"
#include "Buffer/Buffer.h"
#include "Buffer/CreateBuffer.h"
#include "Window/CreateWindow.h"
#include "Objects/Player.h"
#include "Window/InputControl.h"
#include "Shader/CreateShader.h"
#include "Shader/ShaderInfo.h"
#include "DataReader/DataReader.h"
#include "Texture/CreateTexture.h"
#include "Shader/ShaderManager.h"
#include "Definitions.h"

#include "ProgramModes/ProgramParams.h"
#include "ProgramModes/PictureMode.h"
#include "ProgramModes/RealtimeMode.h"
#include "ProgramModes/CPUMode.h"
#include "Scene/SceneDescription.h"
#include "DataReader/DataWriter.h"


#include "../Dependencies/imgui/imgui.h"
#include "../Dependencies/imgui/imgui_impl_glfw.h"
#include "../Dependencies/imgui/imgui_impl_opengl3.h"


int main(int argc, char* argv[])
{
		
	ProgramParams::ComputeProgramArgs(argc,argv);

	ProgramParams::ProgramInit();

	Scene::SceneDescription();
	
	ObjectCreator::SerializeToSSBO(DrawObjectManager::GetSphereQueue());

	for (GameObject& go : GameObjectManager::GetGameObjectList()) {

		ObjectCreator::SerializeToSSBO(go, 16);

	}

	GameObject Canvas = ObjectCreator::Generate2DSquare();

	Renderer::AppendToHUDRenderQueue(Canvas);
	
	if (ProgramParams::picOnly) {
		
		if (ProgramParams::CPURender) {
			CPUMode::Render();
		}
		else {
			PictureMode::Render();
		}

		ProgramParams::window.SwapChain();

		DataWriter::WriteToProgramPathFile();
	}
	else {
		
		RealtimeMode::Render();

	}

	ProgramParams::window.DeleteWindow();
	
	return 0;
}