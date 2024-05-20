#include "WindowCallback.h"

#define MOVESPEED  1.0f
#define NO_CAM_MOVEMENT false

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
	if (NO_CAM_MOVEMENT)
		return;

	glm::vec3 direction = Renderer::GetCamera()->GetDirection();

	Renderer::GetCamera()->rawXposMouse = xpos;
	Renderer::GetCamera()->rawYposMouse = ypos;

	direction.x = cos(glm::radians(xpos)) * cos(glm::radians(-ypos));
	direction.y = sin(glm::radians(ypos));
	direction.z = sin(glm::radians(xpos)) * cos(glm::radians(ypos));	

	Renderer::GetCamera()->SetDirection(glm::normalize(direction));
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	Renderer::GetCamera()->SetZoom(Renderer::GetCamera()->GetZoom() - (float)yoffset);
	if (Renderer::GetCamera()->GetZoom() < 1.0f)
		Renderer::GetCamera()->SetZoom(1.0f);
	if (Renderer::GetCamera()->GetZoom() > 45.0f)
		Renderer::GetCamera()->SetZoom(45.0f);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {


		Camera* camera = Renderer::GetCamera();
		glm::vec3 viewPos = camera->GetPosition();
		float M_PI = 3.14159265359;
		float x = camera->rawXposMouse, y = camera->rawYposMouse;

		float s = x / (camera->imageWidth-1);

		float t = 1-y / (camera->imageHeight-1);
   
		glm::vec3 direction = camera->lower_left_corner + s * camera->horizontal + t * camera->vertical - viewPos;

        double t_max = 100.0;
        bool hit_anything = false;
        double closest_so_far = t_max;
		auto sphereQueue = DrawObjectManager::GetSphereQueue();
		int size = sphereQueue.size();


        for (int i = 0; i < sphereQueue.size(); i++) {
			Sphere currSphere = sphereQueue[i];
            if (CollisionDetection::RaySphereIntersection(currSphere.position,currSphere.radius,viewPos,direction,0.0, closest_so_far)) {
				std::cout << "hit Sphere: " << currSphere.ID << std::endl;
				ObjSelect obj;
				obj.objectID = currSphere.ID;
				obj.objectType = SPHERE;
				Player::SetSelectedObject(obj);

            }
        }
		auto list = GameObjectManager::GetGameObjectList();
		auto queue = DrawObjectManager::GetObjectSlotQueue();

		//for (int i = 0; i < queue.size();i++) {
		for (int i = 0; i < list.size();i++) {
			glm::vec3 intersectPoint;
			auto drawObj = DrawObjectManager::GetDrawObjectByID(list[i].GetDrawObjectID());
			glm::vec3 spherePos = list[i].GetTranslation() * glm::vec4(glm::vec3(0.0f, 0.0f, 0.0f), 1.0f);
			if (!CollisionDetection::RaySphereIntersectionSimple(spherePos, drawObj.shapeData.maxRad, viewPos,direction)) break;
			for (auto& mesh : drawObj.shapeData.meshVerts) {
				mesh[0] = list[i].GetTranslation() * glm::vec4(mesh[0],1.0f);
				mesh[1] = list[i].GetTranslation() * glm::vec4(mesh[1],1.0f);
				mesh[2] = list[i].GetTranslation() * glm::vec4(mesh[2],1.0f);
				if (CollisionDetection::RayTriangleIntersect(viewPos, direction,0.0f , closest_so_far, mesh, intersectPoint)) {
					std::cout << "hit at obj: " << i << "at point: " << glm::to_string(intersectPoint) << std::endl;
					ObjSelect obj;
					obj.objectID = i;
					obj.objectType = TRIANGLE;
					Player::SetSelectedObject(obj);

				}
			}

		}

        // if (hit_anything)
          //  return hit_anything;



        /*Triangle tri;
        Sphere boundingSphere;
        vec3 boundingSpherePos;
        int offsetI = 0;
        for (int o = 0; o < drawObjectDataList.length(); o++) {
            vec3 boundingSpherePos = vec3(drawObjectDataList[o].translation * vec4(0.0f, 0.0f, 0.0f, 1.0f));
            boundingSphere.positionAndradius = vec4(boundingSpherePos, vertArr[0].fuzzAndmat_ptr.a);
            //for (int i = 0; i < indiArr.length(); i += 3) {
            for (int i = int(drawObjectDataList[o].IndiLeftRightVertiLeftRight.x);
                i < int(drawObjectDataList[o].IndiLeftRightVertiLeftRight.y); i += 3) {

                // if (!sphere_hit_simple(r, boundingSphere)) break;



                if (o > 0) {
                    offsetI = int(drawObjectDataList[o].IDType.z);
                }

                tri.position[0] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i]].position;
                tri.position[1] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i + 1]].position;
                tri.position[2] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i + 2]].position;

                tri.normal = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i]].normal;

                tri.albedo = vertArr[indiArr[i]].albedo;
                tri.fuzzAndmat_ptr = vertArr[indiArr[i]].fuzzAndmat_ptr;


                //  tri.fuzzAndmat_ptr = vertArr[indiArr[i]].fuzzAndmat_ptr;

                if (RayTriangleIntersect(r, temp_rec, t_min, closest_so_far, tri, p)) {
                    hit_anything = true;
                    closest_so_far = temp_rec.t;
                    rec = temp_rec;
                    rec.arrIndex = i;
                    rec.geoType = TRIANGLE;
                    rec.tri = tri;
                }
            }
        }
        //return kd_tree_search(t_min, r, t_max, rec);*/
        //return hit_anything;


	}
	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {

	
	}

}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {

	/*const float cameraSpeed = 0.05f; // adjust accordingly
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
	if (glfwGetKey(window, GLFW_KEY_F1) == GLFW_PRESS)
		Renderer::SetRenderMode(GL_FILL);
	if (glfwGetKey(window, GLFW_KEY_F2) == GLFW_PRESS)
		Renderer::SetRenderMode(GL_LINE);
	if (glfwGetKey(window, GLFW_KEY_F3) == GLFW_PRESS)
		Renderer::shadows_active = !Renderer::shadows_active;
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS && player->MoveForwardPossible())
		player->MoveForward();
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS && player->MoveBackwartsPossible())
		player->MoveBackwarts();
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS && player->MoveLeftPossible())
		player->MoveLeft();
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS && player->MoveRightPossible())
		player->MoveRight();
	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
		player->ActivateJump();
	if (glfwGetKey(window, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS) {
		auto cam = Renderer::GetCamera();
		cam->updateData = !cam->updateData;
	}
	if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
		player->MoveUp();
	if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
		player->MoveDown();
	if (glfwGetKey(window, GLFW_KEY_F5) == GLFW_PRESS) {

		std::vector<ShaderInfo> shaderInfo2 = {
			{GL_VERTEX_SHADER, VERTEX_SHADER_PATH_CROSS, NULL },
			{GL_FRAGMENT_SHADER, FRAGMENT_SHADER_PATH_CROSS, NULL}
		};

		std::string shaderName = std::string(HUD_FRAGMENT_SHADER);
		auto shader = CreateShader::GenerateShader(shaderInfo2);

		ShaderManager::RegisterShader(shaderName, shader);

		std::cout << HUD_FRAGMENT_SHADER << " Reloaded" << std::endl;
	}
	if (glfwGetKey(window, GLFW_KEY_PAGE_UP) == GLFW_PRESS)
		Renderer::global_refract_index += .1f;
	if (glfwGetKey(window, GLFW_KEY_PAGE_DOWN) == GLFW_PRESS)
		Renderer::global_refract_index -= .1f;*/
}