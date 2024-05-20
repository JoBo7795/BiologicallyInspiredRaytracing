#pragma once

#include "Texture.h"
#include <string>
#include <iostream>

#include "../Buffer/CreateBuffer.h"
#include "../../Dependencies/stb/stb_image.h"

namespace CreateTexture
{
	Texture CreateTexture(std::string texPath, GLuint programID,bool alpha = false, bool flip = false, int textureMode = GL_TEXTURE_2D, int wrappingMethod = GL_MIRRORED_REPEAT);
};

