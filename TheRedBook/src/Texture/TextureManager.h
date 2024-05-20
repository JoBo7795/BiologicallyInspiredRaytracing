#pragma once

#include <vector>
#include "Texture.h"

class TextureManager
{
private:
	static std::vector<Texture> TextureSlotQueue;
	static std::vector<GLuint> BoundTextures;
	static GLuint boundID;
public:
	static int RegisterTexture(Texture texture);
	static void UseTexture(GLuint TexID);
	static void UnbindTexture();
};

