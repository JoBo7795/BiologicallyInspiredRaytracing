#include "TextureManager.h"

std::vector<Texture> TextureManager::TextureSlotQueue;
std::vector<GLuint> TextureManager::BoundTextures;
GLuint TextureManager::boundID = -1;

int TextureManager::RegisterTexture(Texture texture) {
	TextureSlotQueue.push_back(texture);
	return TextureSlotQueue.size() - 1;
}

void TextureManager::UseTexture(GLuint texID) {

	for (int i = 0; i < TextureSlotQueue.size(); i++) {
		if (i == texID && texID != boundID) {
			TextureSlotQueue[texID].Use();
			boundID = texID;
		}
	}
}

void TextureManager::UnbindTexture() {

	glBindTexture(GL_TEXTURE_2D, 0);
}
