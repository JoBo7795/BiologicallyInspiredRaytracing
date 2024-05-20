#pragma once

#include "GameObject.h"

#include <vector>
#include "../Buffer/CreateBuffer.h"
#include "../Mathematics/Mathematics.h"
#include "DrawObjectManager.h"

namespace ObjectCreator
{
	GameObject GenerateSquare();
	BoundBox GenerateClickPoint(glm::vec3 position);
	BoundBox GeneratePlayerBB(glm::vec3 position);
	GameObject CustomObject(GLuint drawObjectID, GLuint textureID);
	RTGameObject CustomRTObject(GLuint drawObjectID, GLuint textureID);
	GameObject Generate2DSquare();
	void SerializeToSSBO(GameObject& go, int objectSize);
	void SerializeToSSBO(std::vector<float>& data);
	void SerializeToSSBO(std::vector<Sphere>& data);
	ShaderStorageBuffer SerializeToUBO(std::vector<float>& data);


};

