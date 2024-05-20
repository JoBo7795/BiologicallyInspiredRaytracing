#include "ShaderManager.h"

std::vector<Shader> ShaderManager::ShaderStorage;
std::map<std::string, Shader> ShaderManager::ShaderStorageMap;


void ShaderManager::RegisterShader(std::string& shaderName, Shader& shader) {

	ShaderStorageMap[shaderName] = shader;
}

bool ShaderManager::GetRegisteredShader(std::string& shaderName, Shader& shader) {
	if (ShaderStorageMap.find(shaderName) != ShaderStorageMap.end()) {
		shader = ShaderStorageMap[shaderName];
		return true;
	}
	return false;

}

void ShaderManager::RegisterShader(Shader& shader) {
	ShaderStorage.push_back(shader);
}

bool ShaderManager::GetRegisteredShader(int programID, Shader& shader){

	for (auto& shaderIt : ShaderStorage) {
		if (shader.GetProgramID() == programID) {
			shader = shaderIt;
			return true;
		}
	}

	return false;
}