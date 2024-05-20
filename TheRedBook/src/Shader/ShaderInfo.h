#pragma once
#include "GL/glew.h"

struct ShaderInfo {
	GLuint shader_type;	
	const char* filePath;
	const char* shaderCode;
};