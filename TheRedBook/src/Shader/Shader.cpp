#include "Shader.h"

GLuint Shader::GetProgramID() {

	return this->shaderProgramID;
}
void Shader::SetProgramID(GLuint programID) {

	this->shaderProgramID = programID;
}

GLuint Shader::GetVertexID() {

	return this->vertexID;
}

void Shader::SetVertexID(GLuint vertexID) {

	this->vertexID = vertexID;
}

GLuint Shader::GetFragmentID() {

	return this->fragmentID;
}

void Shader::SetFragmentID(GLuint fragmentID) {

	this->fragmentID = fragmentID;
}

void Shader::Use() {

	glUseProgram(this->shaderProgramID);
}

void Shader::SetInt(const char* varName,GLint val) {

	glUniform1i(glGetUniformLocation(this->shaderProgramID, varName), val);

}

void Shader::DeleteShaders() {
	glDeleteShader(this->vertexID);
	glDeleteShader(this->fragmentID);
}

void Shader::DeleteShaderProgram() {

	glDeleteProgram(this->shaderProgramID);
}