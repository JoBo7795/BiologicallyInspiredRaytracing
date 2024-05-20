#pragma once

#include "GL/glew.h"

#include "Buffer.h"
#include "BufferCreationConfig.h"
#include "../Objects/VerticeMap.h"

#include "../../Dependencies/glm/glm/glm.hpp"
#include "../../Dependencies/glm/glm/gtc/matrix_transform.hpp"
#include "../../Dependencies/glm/glm/gtc/type_ptr.hpp"

namespace CreateBuffer
{

	RasterBuffer CreateIndexedBuffer(float* vertices, unsigned int* indices, unsigned int vertiSize, unsigned int indiSize, VerticeMap verticeMap);
	RasterBuffer CreateIndexedBufferCustom(float* vertices, unsigned int* indices, unsigned int vertiSize, unsigned int indiSize, BufferCreationConfig bufferCreationConfig);
	ShaderStorageBuffer GenerateShaderStorageBufferObject(float* data, int size);
	ShaderStorageBuffer GenerateShaderStorageBufferObject(int* data, int size);
	ShaderStorageBuffer GenerateShaderStorageBufferObject(GLuint* data, int size);
	ShaderStorageBuffer GenerateUniformBufferObject(GLuint* data, int size);
	ShaderStorageBuffer GenerateUniformBufferObject(float* data, int size);
	GLuint GenerateDataBuffer(int target,int size);
	void SetVertexAttribPointer(GLuint vao, int target, int bufferId, GLuint index, GLuint size, GLuint stride, GLuint pointer);
	void SetInstanceMatrixPointers(GLuint vao, int target, int bufferId, GLuint index, GLuint size, GLuint stride, GLuint pointer);
	void SetVertexAttribDivisor(GLuint vao, int target, int bufferId, int index, int divisor = 1);
	void SetBufferSubData(int target, GLuint bufferId, int offset, int size, /*glm::mat4*/ float& data);
	GLuint CreateTextureBuffer(unsigned char* data, GLuint width, GLuint height, GLenum format);
	GLuint CreateFrameBufferObject();
};

