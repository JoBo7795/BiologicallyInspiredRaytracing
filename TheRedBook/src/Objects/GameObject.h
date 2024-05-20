#pragma once

#include <vector>
#include <iostream>
#include "../../Dependencies/glm/glm/glm.hpp"
#include "../../Dependencies/glm/glm/gtc/matrix_transform.hpp"
#include "../../Dependencies/glm/glm/gtc/type_ptr.hpp"
#include "../../Dependencies/glm/glm/gtx/string_cast.hpp""

#include "../Texture/Texture.h"
#include "../Buffer/Buffer.h"
#include "BBStruct.h"
#include "../Mathematics/Mathematics.h"
#include "DrawObjectManager.h"

class GameObject
{

private:

	GLuint objectID, DrawObjectID, chunkObjectID;

	unsigned int vertSize, indiSize;

	std::vector<GLuint> textureIDs;
	bool textureSet;
	glm::vec3 position, color;
	bool collision;

protected:
	glm::mat4 translation, rotation, scale;
	BoundBox boundingbox;

public:

	GameObject();
	GameObject(glm::mat4 scale, glm::mat4 translation, glm::mat4 rotation);
	GameObject(GLuint drawObjectID, GLuint textureID);

	~GameObject();

	GLuint GetObjectID();
	void SetObjectID(GLuint objectID);

	GLuint GetDrawObjectID();
	void SetDrawObjectID(GLuint drawObjectID);

	GLuint GetChunkObjectID();
	void SetChunkObjectID(GLuint drawObjectID);

	void SetBuffer(RasterBuffer buffer);

	unsigned int GetVertiSize();
	void SetVertiSize(unsigned int size);

	unsigned int GetIndiSize();
	void SetIndiSize(unsigned int size);

	glm::vec3 GetColor();
	void SetColor(glm::vec3);

	bool GetTextureSet();
	void SetTextureSet(bool set);

	glm::mat4 GetScale();
	void SetScale(glm::mat4);

	void TransformObject();

	glm::mat4 GetRotation();
	void SetRotation(glm::mat4);

	glm::mat4 GetTranslation();
	void SetTranslation(glm::mat4);

	void AddTexture(Texture texture);

	std::vector<GLuint> GetTextureIDs();
	void AddTextureID(GLuint texID);

	BoundBox& GetBoundingBox();
	void SetBoundingBox(BoundBox bbox);

	bool GetCollision();
	void SetCollision(bool val);

};

class RTGameObject
{

private:

	GLuint objectID,drawObjectID;

	unsigned int vertSize, indiSize;

	std::vector<GLuint> textureIDs;
	glm::vec3 position, color;
	bool collision;

protected:
	glm::mat4 translation, rotation, scale;
	BoundBox boundingbox;

public:

	RTGameObject();
	RTGameObject(glm::mat4 scale, glm::mat4 translation, glm::mat4 rotation);
	RTGameObject(GLuint drawObjectID, GLuint textureID);

	~RTGameObject();

	GLuint GetObjectID();
	void SetObjectID(GLuint objectID);

	void SetBuffer(ShaderStorageBuffer buffer);

	unsigned int GetVertiSize();
	void SetVertiSize(unsigned int size);

	unsigned int GetIndiSize();
	void SetIndiSize(unsigned int size);

	glm::mat4 GetScale();
	void SetScale(glm::mat4);

	glm::mat4 GetRotation();
	void SetRotation(glm::mat4);

	glm::mat4 GetTranslation();
	void SetTranslation(glm::mat4);

	void AddTexture(Texture texture);

	std::vector<GLuint> GetTextureIDs();
	void AddTextureID(GLuint texID);

	BoundBox& GetBoundingBox();
	void SetBoundingBox(BoundBox bbox);

	bool GetCollision();
	void SetCollision(bool val);

};

