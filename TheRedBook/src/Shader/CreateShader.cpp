#include "CreateShader.h"

GLuint CreateShader::CompileShader(std::string program, GLuint type) {

	const char* ShaderSource = program.c_str();

	GLuint shaderID;
	GLint  success;
	char infoLog[512];
	shaderID = glCreateShader(type);

	glShaderSource(shaderID, 1, &ShaderSource, NULL);
	glCompileShader(shaderID);
	
	glGetShaderiv(shaderID, GL_COMPILE_STATUS, &success);

	if (!success)
	{
		glGetShaderInfoLog(shaderID, 512, NULL, infoLog);
		if(shaderID == GL_VERTEX_SHADER)
			std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
		if (shaderID == GL_FRAGMENT_SHADER)
			std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
	}	

	return shaderID;
}

GLuint CreateShader::InstallShader(std::vector<GLuint> program) {

	GLuint shaderProgram;	
	GLint  success;
	char infoLog[512];

	shaderProgram = glCreateProgram();

	for(int i = 0; i < program.size(); i++)
		glAttachShader(shaderProgram, program[i]);
	
	glLinkProgram(shaderProgram);	

	glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
	if (!success) {
		glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::PROGRAM::COMPILATION_FAILED\n" << infoLog << std::endl;
	}

	for (int i = 0; i < program.size(); i++)
		glDeleteShader(program[i]);	

	return shaderProgram;
}

const char* CreateShader::ReadShaderCode(const char* filePath){

	std::ifstream mystring(filePath);

	if (!mystring.good()) {
		std::cout << " Shader file is not good! \n";
	}

	std::string tmpString = std::string(std::istreambuf_iterator<char>(mystring), std::istreambuf_iterator<char>());

	char* shaderString = new char[tmpString.size() + 1];

	std::copy(tmpString.begin(), tmpString.end(), shaderString);

	shaderString[tmpString.size()] = '\0';

	return shaderString;
}

Shader CreateShader::GenerateShader(std::vector<ShaderInfo> shaderInfo) {

	std::vector<GLuint> programIDs;
	Shader shader;

	for (ShaderInfo elem : shaderInfo) {
		if (elem.shaderCode == NULL && elem.filePath != NULL) {
			elem.shaderCode = ReadShaderCode(elem.filePath);			
		}

		if (elem.shader_type == GL_VERTEX_SHADER) {
			shader.SetVertexID(CompileShader(elem.shaderCode, elem.shader_type));
			programIDs.push_back(shader.GetVertexID());			
		}

		if (elem.shader_type == GL_FRAGMENT_SHADER) {
			shader.SetFragmentID(CompileShader(elem.shaderCode, elem.shader_type));
			programIDs.push_back(shader.GetFragmentID());			
		}

		delete(elem.shaderCode);
	}

	shader.SetProgramID(InstallShader(programIDs));

	shader.DeleteShaders();

	return shader;
}