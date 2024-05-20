#pragma once

#include "Shader.h"
#include <vector>
#include <map>
#include <string>

class ShaderManager
{

private:
	static std::vector<Shader> ShaderStorage;	
	static std::map<std::string,Shader> ShaderStorageMap;	
	

public:

	static void RegisterShader(std::string& shaderName,Shader& shader);
	static bool GetRegisteredShader(std::string& shaderName, Shader& shader);

	static void RegisterShader(Shader& shader);
	static bool GetRegisteredShader(int programID, Shader& shader);

};

