#include "Physics.h"

void Physics::HandlePhysics(Player& player) {
	BoundBox bb = player.GetBoundingBox();
	/*if (!CollisionDetection::PlayerWorldCollisionDetection(&bb)) {

		player.CustomMove(glm::vec3(0.0f,-1.0f,0.0f),GRAVITY);
		player.SetCanJump(false);
	}
	else {
		player.SetCanJump(true);
	}*/
}