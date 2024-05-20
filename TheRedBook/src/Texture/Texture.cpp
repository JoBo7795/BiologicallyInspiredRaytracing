#include "Texture.h"

void Texture::SetBorderColor(GLfloat colors) {

	this->borderColor = new GLfloat(colors);
	glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, this->borderColor);
}

GLfloat* Texture::GetBorderColor() {

	return this->borderColor;
}

GLint Texture::GetTexID() {

	return this->texID;
}

void Texture::SetTexID(GLint texID) {

	this->texID = texID;
}

GLint Texture::GetWidth() {

	return this->width;
}

void Texture::SetWidth(GLint width) {

	this->width = width;
}

GLint Texture::GetHeight() {

	return this->height;
}

void Texture::SetHeight(GLint height) {

	this->height = height;
}

void Texture::Use() {

	glBindTexture(GL_TEXTURE_2D, this->texID);
}

void Texture::Unbind() {

	glBindTexture(GL_TEXTURE_2D, 0);
}