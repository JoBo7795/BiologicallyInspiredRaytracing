#pragma once

#include "GL/glew.h"

class Shader
{
	GLuint vertexID,fragmentID,shaderProgramID;

public:

	GLuint GetProgramID();
	void SetProgramID(GLuint programID);

	GLuint GetVertexID();
	void SetVertexID(GLuint vertexID);

	GLuint GetFragmentID();
	void SetFragmentID(GLuint fragmentID);

	void Use();

	void SetInt(const char* varName, GLint val);

	void DeleteShaderProgram();
	void DeleteShaders();
};

