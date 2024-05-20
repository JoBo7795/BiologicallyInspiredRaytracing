#pragma once

#include "GL/glew.h"

class Texture
{

private:

	GLfloat* borderColor;
	GLint width, height, texID;
public:
	void SetBorderColor(GLfloat colors);
	GLfloat* GetBorderColor();

	GLint GetTexID();
	void SetTexID(GLint texID);

	GLint GetWidth();
	void SetWidth(GLint texID);

	GLint GetHeight();
	void SetHeight(GLint texID);

	void Use();
	void Unbind();

};

