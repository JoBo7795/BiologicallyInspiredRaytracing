#pragma once

#include "../Definitions.h"
#include "../Renderer/Renderer.h"
#include "Window.h"
#include "../Objects/ObjectCreator.h"
#include "../Objects/Player.h"
#include "../Shader/CreateShader.h"

#include "../Shader/ShaderManager.h"


//static Shader shaderProgram1, shaderProgram2;

namespace InputControl
{
	
	void processInput(GLFWwindow* window, Player* player);	

};

