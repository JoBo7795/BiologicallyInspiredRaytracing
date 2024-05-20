#pragma once

#include "GameObject.h"
#include "ObjectCreator.h"
#include "../Camera/Camera.h"
#include "../Physics/CollisionDetection.h"

struct ObjSelect {


	int objectID, objectType;


};

class Player :
	public GameObject

{
private:
	Camera* playerCamera;
	float playerSpeed;
	glm::vec3 playerPosition, jumpHeight, jumpMaxHeight;
	bool jumpActive, canJump;
	static ObjSelect selectedObject;

public:

	Player(Camera* camera, glm::vec3 playerSpawn = glm::vec3(0.0f, 0.0f, 0.0f), glm::mat4 scale = glm::mat4(1.0f), glm::mat4 translation = glm::mat4(1.0f), glm::mat4 rotation = glm::mat4(1.0f));

	Camera* GetCamera();
	void SetCamera(Camera* camera);

	void MoveForward();
	void MoveBackwarts();
	void MoveLeft();
	void MoveRight();
	void MoveDown();
	void MoveUp();
	void CustomMove(glm::vec3 direction, float speed);

	void HandleJump();
	void ActivateJump();

	bool MoveForwardPossible();
	bool MoveBackwartsPossible();
	bool MoveLeftPossible();
	bool MoveRightPossible();

	void SetCanJump(bool canJump);
	bool GetCanJump();

	static void SetSelectedObject(ObjSelect selectedObject);
	static ObjSelect GetSelectedObject();

	glm::vec3 GetPlayerPosition();
	void SetPlayerPosition(glm::vec3 position);
};

