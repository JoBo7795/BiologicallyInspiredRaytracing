#include "Player.h"

ObjSelect Player::selectedObject;// = ObjSelect(-1, -1);


Player::Player(Camera* camera, glm::vec3 playerSpawn, glm::mat4 scale, glm::mat4 translation, glm::mat4 rotation)
	:GameObject(scale, translation, rotation) {
	this->playerCamera = camera;
	this->playerPosition = playerSpawn;
	this->playerCamera->SetPosition(playerSpawn);
	this->boundingbox = ObjectCreator::GeneratePlayerBB(playerSpawn);
	this->playerSpeed = 0.05f;
	this->jumpActive = false;
	this->jumpHeight = glm::vec3(0.0f, 0.2f, 0.0f);
	Player::selectedObject.objectID = -1;
	Player::selectedObject.objectType = -1;
}

Camera* Player::GetCamera() {
	return this->playerCamera;
}

void Player::SetCamera(Camera* camera) {
	this->playerCamera = camera;
}

void Player::MoveForward() {
	this->playerPosition = playerCamera->GetPosition() + (playerSpeed * glm::vec3(playerCamera->GetDirection().x, 0.0f, playerCamera->GetDirection().z));
	this->playerCamera->SetPosition(playerPosition);
	this->boundingbox = ObjectCreator::GeneratePlayerBB(playerCamera->GetPosition());
}

void Player::MoveBackwarts() {
	this->playerPosition = playerCamera->GetPosition() - (playerSpeed * glm::vec3(playerCamera->GetDirection().x, 0.0f, playerCamera->GetDirection().z));
	this->playerCamera->SetPosition(playerPosition);
	this->boundingbox = ObjectCreator::GeneratePlayerBB(playerCamera->GetPosition());
}

void Player::MoveLeft() {
	this->playerPosition = playerCamera->GetPosition() - glm::normalize(glm::cross(playerCamera->GetDirection(), playerCamera->GetUp())) * playerSpeed;
	this->playerCamera->SetPosition(playerPosition);
	this->boundingbox = ObjectCreator::GeneratePlayerBB(playerCamera->GetPosition());
}

void Player::MoveRight() {
	this->playerPosition = playerCamera->GetPosition() + glm::normalize(glm::cross(playerCamera->GetDirection(), playerCamera->GetUp())) * playerSpeed;
	this->playerCamera->SetPosition(playerPosition);
	this->boundingbox = ObjectCreator::GeneratePlayerBB(playerCamera->GetPosition());
}

void Player::MoveDown() {
	this->playerPosition = playerCamera->GetPosition() - playerCamera->GetUp() * playerSpeed;
	this->playerCamera->SetPosition(this->playerPosition);
	this->boundingbox = ObjectCreator::GeneratePlayerBB(playerCamera->GetPosition());
}

void Player::MoveUp() {
	this->playerPosition = playerCamera->GetPosition() + playerCamera->GetUp() * playerSpeed;
	this->playerCamera->SetPosition(this->playerPosition);
	this->boundingbox = ObjectCreator::GeneratePlayerBB(playerCamera->GetPosition());
}

void Player::CustomMove(glm::vec3 direction, float speed) {
	this->playerPosition = playerCamera->GetPosition() + direction * speed;
	this->playerCamera->SetPosition(this->playerPosition);
	this->boundingbox = ObjectCreator::GeneratePlayerBB(playerCamera->GetPosition());
}

void Player::ActivateJump() {
	if (!this->jumpActive && this->canJump) {
		this->jumpActive = true;
		this->jumpMaxHeight = this->playerPosition + this->jumpHeight;
	}
}

void Player::HandleJump() {

	if (this->jumpActive && this->jumpMaxHeight.y > this->playerPosition.y) {
		this->playerPosition = playerCamera->GetPosition() + playerCamera->GetUp() * playerSpeed;
		this->playerCamera->SetPosition(this->playerPosition);
		this->boundingbox = ObjectCreator::GeneratePlayerBB(playerCamera->GetPosition());
	}
	else if (this->jumpMaxHeight.y < this->playerPosition.y) {
		this->jumpActive = false;
	}
}

bool Player::MoveForwardPossible() {
	BoundBox clickPoint = ObjectCreator::GenerateClickPoint(playerCamera->GetPosition() + (0.2f * glm::vec3(playerCamera->GetDirection().x, 0.0f, Renderer::GetCamera()->GetDirection().z)));
	//return !CollisionDetection::PlayerWorldCollisionDetection(&clickPoint);
	return true;
}

bool Player::MoveBackwartsPossible() {
	BoundBox clickPoint = ObjectCreator::GenerateClickPoint(playerCamera->GetPosition() - (0.2f * glm::vec3(playerCamera->GetDirection().x, 0.0f, Renderer::GetCamera()->GetDirection().z)));
	//return !CollisionDetection::PlayerWorldCollisionDetection(&clickPoint);
	return true;
}

bool Player::MoveLeftPossible() {
	BoundBox clickPoint = ObjectCreator::GenerateClickPoint(Renderer::GetCamera()->GetPosition() - glm::normalize(glm::cross(Renderer::GetCamera()->GetDirection(), Renderer::GetCamera()->GetUp())) * 0.2f);
	//return !CollisionDetection::PlayerWorldCollisionDetection(&clickPoint);
	return true;
}

bool Player::MoveRightPossible() {
	BoundBox clickPoint = ObjectCreator::GenerateClickPoint(Renderer::GetCamera()->GetPosition() + glm::normalize(glm::cross(Renderer::GetCamera()->GetDirection(), Renderer::GetCamera()->GetUp())) * 0.2f);
	//return !CollisionDetection::PlayerWorldCollisionDetection(&clickPoint);
	return true;
}

void Player::SetCanJump(bool canJump) {
	this->canJump = canJump;
}

bool Player::GetCanJump() {
	return this->canJump;
}

glm::vec3 Player::GetPlayerPosition() {
	return this->playerPosition;
}

void Player::SetPlayerPosition(glm::vec3 position) {
	this->playerPosition = position;
}

void Player::SetSelectedObject(ObjSelect selectedObject_in) {
	selectedObject = selectedObject_in;
}

ObjSelect Player::GetSelectedObject() {
	return selectedObject;
}