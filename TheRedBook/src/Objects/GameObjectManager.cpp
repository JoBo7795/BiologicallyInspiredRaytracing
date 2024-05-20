#include "GameObjectManager.h"

std::vector<GameObject> GameObjectManager::gameObjectList;

GLuint GameObjectManager::AppendGameObject(GameObject gameObject) {

	gameObjectList.push_back(gameObject);
	return gameObjectList.size() - 1;

}

GameObject& GameObjectManager::GetGameObjectByID(GLuint gameObjectID) {

	for (int i = 0; i < gameObjectList.size(); i++) {
		if (i == gameObjectID) {
			return gameObjectList[gameObjectID];
		}
	}
}

std::vector<GameObject> GameObjectManager::GetGameObjectList() {
	return gameObjectList;
}