#pragma once

#include "GL/glew.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "ShaderInfo.h"
#include "Shader.h"

namespace CreateShader
{
	GLuint CompileShader(std::string program, GLuint type);
	GLuint InstallShader(std::vector<GLuint> program);
	Shader GenerateShader(std::vector<ShaderInfo> shaderInfo);
	const char* ReadShaderCode(const char* filePath);
};

