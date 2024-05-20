#include "CreateBuffer.h"

int offset = 0;
RasterBuffer CreateBuffer::CreateIndexedBuffer(float* vertices, unsigned int* indices, unsigned int vertiSize, unsigned int indiSize, VerticeMap verticeMap) {

	RasterBuffer indexedBuffer;
	unsigned int vbo, vao, ibo;

	glGenVertexArrays(1, &vao);
	std::cout << "glGenVertexError: " << glGetError() << std::endl;
	glGenBuffers(1, &vbo);
	std::cout << "glGenBuffersError: " << glGetError() << std::endl;
	glGenBuffers(1, &ibo);
	std::cout << "glGenBuffersError: " << glGetError() << std::endl;
	glBindVertexArray(vao);
	std::cout << "glBindBufferError: " << glGetError() << std::endl;
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	std::cout << "glBindBufferError: " << glGetError() << std::endl;
	glBufferData(GL_ARRAY_BUFFER, vertiSize, vertices, GL_STATIC_DRAW);
	std::cout << "glBufferDataError: " << glGetError() << std::endl;
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	std::cout << "glBindBufferError: " << glGetError() << std::endl;
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indiSize, indices, GL_STATIC_DRAW);
	std::cout << "glBufferDataError: " << glGetError() << std::endl;
	//glBufferData(GL_ARRAY_BUFFER, vertiSize, vertices, GL_DYNAMIC_DRAW);
	//std::cout << "glBufferDataError: " << glGetError() << std::endl;
	glVertexAttribPointer(0, verticeMap.positionSize, GL_FLOAT, GL_FALSE, verticeMap.stride * sizeof(float), (void*)0);
	std::cout << "glVertexAttribPointerError0: " << glGetError() << std::endl;
	glVertexAttribPointer(1,/* verticeMap.colorSize*/3, GL_FLOAT, GL_FALSE, verticeMap.stride * sizeof(float), (void*)((verticeMap.positionSize + verticeMap.normalSize) * sizeof(float)));
	std::cout << "glVertexAttribPointerError1: " << glGetError() << std::endl;
	glVertexAttribPointer(2, verticeMap.textureSize, GL_FLOAT, GL_FALSE, verticeMap.stride * sizeof(float), (void*)((verticeMap.positionSize + verticeMap.normalSize + verticeMap.colorSize) * sizeof(float)));
	std::cout << "glVertexAttribPointerError2: " << glGetError() << std::endl;
	glVertexAttribPointer(3, /*verticeMap.normalSize*/3, GL_FLOAT, GL_FALSE, verticeMap.stride * sizeof(float), (void*)((verticeMap.positionSize) * sizeof(float)));
	std::cout << "glVertexAttribPointerError3: " << glGetError() << std::endl;
	offset = (verticeMap.positionSize + verticeMap.normalSize + verticeMap.colorSize) * sizeof(float);

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);
	glEnableVertexAttribArray(3);
	std::cout << "glEnableVertexAttribArrayError: " << glGetError() << std::endl;

	indexedBuffer.SetIBO(ibo);
	indexedBuffer.AddVBO(vbo);
	indexedBuffer.SetVAO(vao);

	return indexedBuffer;
}

RasterBuffer CreateBuffer::CreateIndexedBufferCustom(float* vertices, unsigned int* indices, unsigned int vertiSize, unsigned int indiSize, BufferCreationConfig bufferCreationConfig) {

	RasterBuffer indexedBuffer;
	unsigned int vbo, vao, ibo;

	glGenVertexArrays(1, &vao);
	std::cout << "glGenVertexError: " << glGetError() << std::endl;
	glGenBuffers(1, &vbo);
	std::cout << "glGenBuffersError: " << glGetError() << std::endl;
	glGenBuffers(1, &ibo);
	std::cout << "glGenBuffersError: " << glGetError() << std::endl;
	glBindVertexArray(vao);
	std::cout << "glBindBufferError: " << glGetError() << std::endl;
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	std::cout << "glBindBufferError: " << glGetError() << std::endl;
	glBufferData(GL_ARRAY_BUFFER, vertiSize, vertices, GL_STATIC_DRAW);
	std::cout << "glBufferDataError: " << glGetError() << std::endl;
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	std::cout << "glBindBufferError: " << glGetError() << std::endl;
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indiSize, indices, GL_STATIC_DRAW);
	std::cout << "glBufferDataError: " << glGetError() << std::endl;

	int pointerOffset = 0;

	for (VertObj vertObj : bufferCreationConfig.vertObj) {

		glVertexAttribPointer(vertObj.PointerID, vertObj.size, vertObj.dataType, GL_FALSE, vertObj.stride * sizeof(float), (void*)(pointerOffset+vertObj.size));
		std::cout << "glVertexAttribPointerError"<< vertObj.PointerID<<": " << glGetError() << std::endl;

		glEnableVertexAttribArray(vertObj.PointerID);
		std::cout << "glEnableVertexAttribArrayError"<< vertObj.PointerID << ": " << glGetError() << std::endl;

		pointerOffset += vertObj.size;
	}

	

	indexedBuffer.SetIBO(ibo);
	indexedBuffer.AddVBO(vbo);
	indexedBuffer.SetVAO(vao);

	return indexedBuffer;
}


ShaderStorageBuffer CreateBuffer::GenerateShaderStorageBufferObject(float* data, int size) {

	ShaderStorageBuffer ssbo;

	GLuint ssboID;
	glGenBuffers(1, &ssboID);
	std::cout << "before:" << glGetError() << std::endl;
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboID);
	std::cout << "bind:" << glGetError() << std::endl;
	glBufferData(GL_SHADER_STORAGE_BUFFER, size, data, GL_DYNAMIC_DRAW);
	std::cout << "bufferdata:" << glGetError() << std::endl;
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, ShaderStorageManager::GetCounter(), ssboID);
	std::cout << "bufferbase:" << glGetError() << std::endl;
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	std::cout << "after:" << glGetError() << std::endl;
	ssbo.SetBufferID(ssboID);
	ssbo.SetBufferIndex(ShaderStorageManager::GetCounter());

	ShaderStorageManager::IncreaseCounter();

	return ssbo;
}

ShaderStorageBuffer CreateBuffer::GenerateShaderStorageBufferObject(int* data, int size) {

	ShaderStorageBuffer ssbo;

	GLuint ssboID;
	glGenBuffers(1, &ssboID);
	std::cout <<"before:"<< glGetError() << std::endl;
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboID);
	glBufferData(GL_SHADER_STORAGE_BUFFER, size, data, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, ShaderStorageManager::GetCounter(), ssboID);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	std::cout << "after:" << glGetError() << std::endl;
	ssbo.SetBufferID(ssboID);
	ssbo.SetBufferIndex(ShaderStorageManager::GetCounter());

	ShaderStorageManager::IncreaseCounter();

	return ssbo;
}

ShaderStorageBuffer CreateBuffer::GenerateShaderStorageBufferObject(GLuint* data, int size) {

	ShaderStorageBuffer ssbo;

	GLuint ssboID;
	glGenBuffers(1, &ssboID);
	std::cout << "before:" << glGetError() << std::endl;
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboID);
	glBufferData(GL_SHADER_STORAGE_BUFFER, size, data, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, ShaderStorageManager::GetCounter(), ssboID);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	std::cout << "after:" << glGetError() << std::endl;
	ssbo.SetBufferID(ssboID);
	ssbo.SetBufferIndex(ShaderStorageManager::GetCounter());

	ShaderStorageManager::IncreaseCounter();

	return ssbo;
}

ShaderStorageBuffer CreateBuffer::GenerateUniformBufferObject(GLuint* data, int size) {

	ShaderStorageBuffer ubo;

	GLuint uboID;
	glGenBuffers(1, &uboID);
	std::cout << "before:" << glGetError() << std::endl;
	glBindBuffer(GL_UNIFORM_BUFFER, uboID);
	glBufferData(GL_UNIFORM_BUFFER, size, data, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_UNIFORM_BUFFER, UniformBufferManager::GetCounter(), uboID);
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	std::cout << "after:" << glGetError() << std::endl;
	ubo.SetBufferID(uboID);
	ubo.SetBufferIndex(UniformBufferManager::GetCounter());

	UniformBufferManager::IncreaseCounter();

	return ubo;
}

ShaderStorageBuffer CreateBuffer::GenerateUniformBufferObject(float* data, int size) {

	ShaderStorageBuffer ubo;

	GLuint uboID;
	glGenBuffers(1, &uboID);
	std::cout << "before:" << glGetError() << std::endl;
	glBindBuffer(GL_UNIFORM_BUFFER, uboID);
	glBufferData(GL_UNIFORM_BUFFER, size, data, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_UNIFORM_BUFFER, UniformBufferManager::GetCounter(), uboID);
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	std::cout << "after:" << glGetError() << std::endl;
	ubo.SetBufferID(uboID);
	ubo.SetBufferIndex(UniformBufferManager::GetCounter());

	UniformBufferManager::IncreaseCounter();

	return ubo;
}

GLuint CreateBuffer::GenerateDataBuffer(int target,int size) {

	GLuint bufferId;
	glGenBuffers(1,&bufferId);
	//std::cout << "glGenBuffers: " << glGetError() << std::endl;
	glBindBuffer(target, bufferId);
	//std::cout << "glBindBuffer: " << glGetError() << std::endl;
	glBufferData(target,size,0,GL_DYNAMIC_DRAW);
	//std::cout << "glBufferData: " << glGetError() << std::endl;

	return bufferId;
}

void CreateBuffer::SetBufferSubData(int target, GLuint bufferId, int offset, int size, /*glm::mat4&*/ float& data) {

	glBindBuffer(target, bufferId);

	void* ptr = glMapBuffer(target, GL_WRITE_ONLY);

	if (ptr == NULL) {
		return;
	}
	
	memcpy(ptr, &data, size);

	glUnmapBuffer(target);

}

void CreateBuffer::SetInstanceMatrixPointers(GLuint vao, int target, int bufferId, GLuint index, GLuint size, GLuint stride, GLuint pointer) {

	glBindVertexArray(vao);
	//std::cout << "glBindVertexArrayError: " << glGetError() << std::endl;
	glBindBuffer(target, bufferId);
	//std::cout << "glBindBufferError: " << glGetError() << std::endl;
	//glBufferSubData(target, 64, size, &data);

	// Attributes for the mat4 matrix
	glVertexAttribPointer(4, size, GL_FLOAT, GL_FALSE, stride, (void*)0);
	//std::cout << "glVertexAttribPointerError4: " << glGetError() << std::endl;
	glEnableVertexAttribArray(4);
	glVertexAttribPointer(5, size, GL_FLOAT, GL_FALSE, /*4 * vec4Size*/stride, (void*)((1) * pointer));
	//std::cout << "glVertexAttribPointerError5: " << glGetError() << std::endl;
	glEnableVertexAttribArray(5);
	glVertexAttribPointer(6, size, GL_FLOAT, GL_FALSE, /*4 * vec4Size*/stride, (void*)((2) * pointer));
	//std::cout << "glVertexAttribPointerError6: " << glGetError() << std::endl;
	glEnableVertexAttribArray(6);
	glVertexAttribPointer(7, size, GL_FLOAT, GL_FALSE, /*4 * vec4Size*/stride, (void*)((3) * pointer));
	//std::cout << "glVertexAttribPointerError7: " << glGetError() << std::endl;
	glEnableVertexAttribArray(7);

	glVertexAttribDivisor(4, 1);
	//std::cout << "glVertexAttribDivisorError4: " << glGetError() << std::endl;
	glVertexAttribDivisor(5, 1);
	//std::cout << "glVertexAttribDivisorError5: " << glGetError() << std::endl;
	glVertexAttribDivisor(6, 1);
	//std::cout << "glVertexAttribDivisorError6: " << glGetError() << std::endl;
	glVertexAttribDivisor(7, 1);
	//std::cout << "glVertexAttribDivisorError7: " << glGetError() << std::endl;

}

void CreateBuffer::SetVertexAttribPointer(GLuint vao, int target, int bufferId, GLuint index, GLuint size, GLuint stride, GLuint pointer) {

	glBindVertexArray(vao);

	glBindBuffer(target, bufferId);

	glVertexAttribPointer(index, size, GL_FLOAT, GL_FALSE, stride, (void*)pointer);

	glEnableVertexAttribArray(index);
	
}

void CreateBuffer::SetVertexAttribDivisor(GLuint vao, int target, int bufferId, int index, int divisor) {

	glBindVertexArray(vao);

	glBindBuffer(target, bufferId);

	glVertexAttribDivisor(index, divisor);

}

GLuint CreateBuffer::CreateTextureBuffer(unsigned char* data, GLuint width, GLuint height,GLenum format) {

	GLuint textureID;
	glGenTextures(1, &textureID);
	//std::cout << "glGenTexturesError: " << glGetError() << std::endl;
	glBindTexture(GL_TEXTURE_2D, textureID);
	//std::cout << "glBindTexturesError: " << glGetError() << std::endl;
	glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
	//std::cout << "glTexImage2DError: " << glGetError() << std::endl;
	glGenerateMipmap(GL_TEXTURE_2D);
	//std::cout << "glGenerateMipmap: " << glGetError() << std::endl;
	return textureID;
}

GLuint CreateBuffer::CreateFrameBufferObject() {

	// Erstelle und binde das FBO
	unsigned int fbo;
	glGenFramebuffers(1, &fbo);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);

	// Erstelle ein Texture-Attachment
	unsigned int texColorBuffer;
	glGenTextures(1, &texColorBuffer);
	glBindTexture(GL_TEXTURE_2D, texColorBuffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 800, 600, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texColorBuffer, 0);

	// Erstelle ein Renderbuffer-Attachment
	unsigned int rbo;
	glGenRenderbuffers(1, &rbo);
	glBindRenderbuffer(GL_RENDERBUFFER, rbo);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, 800, 600);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rbo);

	// Überprüfe die Vollständigkeit des FBO
	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
		std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	return fbo;

}