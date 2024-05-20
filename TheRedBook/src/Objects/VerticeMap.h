#pragma once

struct VerticeMap
{
	unsigned int positionSize, normalSize, colorSize, textureSize, stride;

	VerticeMap(unsigned int posSize = 0, unsigned int normSize = 0, unsigned int colorSize = 0, unsigned int texSize = 0) : positionSize(posSize), normalSize(normSize), colorSize(colorSize), textureSize(texSize) {
		stride = positionSize + normalSize + colorSize + textureSize;
	}
};