#include "CreateTexture.h"

Texture CreateTexture::CreateTexture(std::string texPath, GLuint programID,bool alpha, bool flip,int textureMode,int wrappingMethod) {	

	Texture texture;
	int width, height, colorChannels, texID = -1;
	const char* tmpPath = texPath.c_str();
	std::cout << "texturePath: " << tmpPath << std::endl;
	stbi_set_flip_vertically_on_load(flip);

	unsigned char* data;

	if(!alpha)
		data = stbi_load(tmpPath, &width, &height, &colorChannels, 3);
	else
		data = stbi_load(tmpPath, &width, &height, &colorChannels, 4);

	texture.SetWidth(width);
	texture.SetHeight(height);

	if (data) {
		
		if(!alpha)
			texID = CreateBuffer::CreateTextureBuffer(data, width, height,GL_RGB);
		if (alpha)
			texID = CreateBuffer::CreateTextureBuffer(data, width, height,GL_RGBA);

		glTexParameteri(textureMode, GL_TEXTURE_WRAP_S, wrappingMethod);
		glTexParameteri(textureMode, GL_TEXTURE_WRAP_T, wrappingMethod);

		glTexParameteri(textureMode, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(textureMode, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	}
	else {
		std::cout << "Failed to load texture" << std::endl;
	}

	stbi_image_free(data);

	texture.SetTexID(texID);

	return texture;
}
