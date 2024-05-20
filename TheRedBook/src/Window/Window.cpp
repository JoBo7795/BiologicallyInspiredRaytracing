#include "Window.h"

GLFWwindow* Window::GetWindowRef() {

	return this->window;
}

void Window::SetWindowRef(GLFWwindow* windowRef) {
	this->window = windowRef;
}

void Window::DeleteWindow() {

	glfwTerminate();
}

void Window::SwapChain() {

	glfwSwapBuffers(this->GetWindowRef());
	glfwPollEvents();
}

bool Window::WindowShouldClose() {
	return glfwWindowShouldClose(this->GetWindowRef());
}

void Window::SetCallback(GLuint callbackID) {

	switch (callbackID)
	{
	case MOUSE_INPUT:
		glfwSetCursorPosCallback(window, mouse_callback);
		 
	case SCROLL_INPUT:
		glfwSetScrollCallback(window, scroll_callback);

	case MOUSE_CLICK:
		glfwSetMouseButtonCallback(window, mouse_button_callback);

	case KEY_INPUT:
		glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	case BUFFER_SIZE:
		glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	default:
		break;
	}

}

void Window::SetWidth(GLfloat width) {
	this->width = width;
}

GLfloat Window::GetWidth() {
	return this->width;
}

void Window::SetHeigth(GLfloat height) {
	this->height = height;
}

GLfloat Window::GetHeigth() {
	return this->height;
}