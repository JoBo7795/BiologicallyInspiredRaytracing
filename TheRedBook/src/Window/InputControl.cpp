#include "InputControl.h"
bool rAltReleased = true;
bool ffReleased = true;
void InputControl::processInput(GLFWwindow* window, Player* player)
{
	auto rendererReference = Renderer::GetInstance();

	const float cameraSpeed = 0.05f; // adjust accordingly
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
	if (glfwGetKey(window, GLFW_KEY_F1) == GLFW_PRESS)
		rendererReference->SetRenderMode(GL_FILL);
	if (glfwGetKey(window, GLFW_KEY_F2) == GLFW_PRESS)
		rendererReference->SetRenderMode(GL_LINE);
	if (glfwGetKey(window, GLFW_KEY_F3) == GLFW_PRESS)
		rendererReference->SetShadowsActive(!rendererReference->GetShadowsActive());
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
	if (glfwGetKey(window, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS && rAltReleased) {
		auto cam = rendererReference->GetCamera();
		cam->updateData = !cam->updateData;
		rAltReleased = false;
	}
	if (glfwGetKey(window, GLFW_KEY_RIGHT_ALT) == GLFW_RELEASE) {
		rAltReleased = true;
	}
	if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
		player->MoveUp();
	if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS)
		player->MoveDown();
	if (glfwGetKey(window, GLFW_KEY_F5) == GLFW_PRESS &&ffReleased) {

		std::vector<ShaderInfo> shaderInfo2 = {
			{GL_VERTEX_SHADER, VERTEX_SHADER_PATH_RAYTRACE, NULL },
			{GL_FRAGMENT_SHADER, FRAGMENT_SHADER_PATH_RAYTRACE, NULL}
		};

		std::string shaderName = std::string(HUD_FRAGMENT_SHADER);
		auto shader = CreateShader::GenerateShader(shaderInfo2);
		
		ShaderManager::RegisterShader( shaderName,shader);
		ffReleased = false;
		std::cout << HUD_FRAGMENT_SHADER << " Reloaded" << std::endl;		
	}
	if (glfwGetKey(window, GLFW_KEY_F5) == GLFW_RELEASE && !ffReleased) {

		ffReleased = true;
	}
	if (glfwGetKey(window, GLFW_KEY_PAGE_UP) == GLFW_PRESS) {
		float currVal = rendererReference->GetGlobalRefractIndex();
		rendererReference->SetGlobalRefractIndex((currVal += .1f));
	}

	
	if (glfwGetKey(window, GLFW_KEY_PAGE_DOWN) == GLFW_PRESS) {
		float currVal = rendererReference->GetGlobalRefractIndex();
		rendererReference->SetGlobalRefractIndex((currVal -= .1f));
	}

	if (rendererReference->GetGlobalRefractIndex() < 0.0f)
		rendererReference->SetGlobalRefractIndex(0.0f);
}