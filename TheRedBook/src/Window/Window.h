#pragma once

#include<direct.h>
#include "../Renderer/Renderer.h"
#include "WindowCallback.h"
#include "../../Dependencies/stb/stb_image_write.h"


class Window
{
private:
	GLFWwindow* window;	
	GLfloat height, width;
	

public:	
	GLFWwindow* GetWindowRef();
	void SetWindowRef(GLFWwindow* windowRef);

	void DeleteWindow();
	void SwapChain();
	bool WindowShouldClose();
	void SetCallback(GLuint callbackID);

	void SetWidth(GLfloat width);
	GLfloat GetWidth();

	void SetHeigth(GLfloat height);
	GLfloat GetHeigth();

	inline void saveImage(char* filepath, GLFWwindow* w) {
		int width, height;
		glfwGetFramebufferSize(w, &width, &height);
		GLsizei nrChannels = 3;
		GLsizei stride = nrChannels * width;
		stride += (stride % 4) ? (4 - stride % 4) : 0;
		GLsizei bufferSize = stride * height;
		std::vector<char> buffer(bufferSize);

		// read pixel from gpu
		glPixelStorei(GL_PACK_ALIGNMENT, 4);
		glReadBuffer(GL_FRONT);
		glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());
		stbi_flip_vertically_on_write(true);

		auto path = std::string(filepath);
		struct stat info;
		size_t pos = 0;
		while ((pos = path.find_first_of("\\/", pos + 1)) != std::string::npos) {
			std::string subPath = path.substr(0, pos);
			if (stat(subPath.c_str(), &info) != 0) {
				std::cout << "Creating directory: " << subPath << std::endl;
				_mkdir(subPath.c_str());
			}
		}


		if (stbi_write_png(filepath, width, height, nrChannels, buffer.data(), stride) == 0)
			assert("image not written");
	}

	inline void saveImageCPU(char* filepath, uint8_t* data) {
		
		GLsizei nrChannels = 3;
		GLsizei stride = nrChannels * width;
		stride += (stride % 4) ? (4 - stride % 4) : 0;
		GLsizei bufferSize = stride * height;
		std::vector<char> buffer(bufferSize);

		stbi_flip_vertically_on_write(true);

		auto path = std::string(filepath);
		struct stat info;
		size_t pos = 0;
		while ((pos = path.find_first_of("\\/", pos + 1)) != std::string::npos) {
			std::string subPath = path.substr(0, pos);
			if (stat(subPath.c_str(), &info) != 0) {
				std::cout << "Creating directory: " << subPath << std::endl;
				_mkdir(subPath.c_str());
			}
		}


		if (stbi_write_png(filepath, this->width, this->height, nrChannels, data, stride) == 0)
			assert("image not written");
	}

	float mouseRawXpos, mouseRawYpos;

};

