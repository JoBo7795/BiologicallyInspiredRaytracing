#include "CreateWindow.h"

Window CreateWindow::CreateWindow(float windowWidth,float windowHeight) {

	Window windowObject;

	glfwInit();

	glfwWindowHint(GLFW_VISIBLE, true);
	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, true);

	windowObject.SetWindowRef(glfwCreateWindow(windowWidth, windowHeight, "MeinCraft_RTXOn_DLSS5.0_und_GTA6", NULL, NULL));
	if (windowObject.GetWindowRef() == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
	}

	glfwMakeContextCurrent(windowObject.GetWindowRef());	
	glViewport(0, 0, windowWidth, windowHeight);

	return windowObject;
}