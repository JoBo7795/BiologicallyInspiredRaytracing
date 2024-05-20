#pragma once

#include "GameObject.h"

class GameObjectManager
{
private:
	static std::vector<GameObject> gameObjectList;

public:
	GameObjectManager();

	static GLuint AppendGameObject(GameObject gameObject);
	static GameObject& GetGameObjectByID(GLuint gameObjectID);

	static std::vector<GameObject> GetGameObjectList();

};

