#pragma once

#include "Window.h"
#include "../Objects/ObjectCreator.h"
#include "../Physics/CollisionDetection.h"
#include "../../glm/glm/gtx/vector_angle.hpp"
#include "../Physics/CollisionDetection.h"
#include "../Objects/DrawObjectManager.h"

#include "../Definitions.h"
#include "../Renderer/Renderer.h"
#include "Window.h"
#include "../Objects/ObjectCreator.h"
#include "../Objects/Player.h"
#include "../Shader/CreateShader.h"

#include "../Shader/ShaderManager.h"
#include "../Objects/Player.h"
#include "../Objects/GameObjectManager.h"


void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);