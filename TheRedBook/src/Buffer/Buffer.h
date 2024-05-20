#pragma once

#include "GL/glew.h"
#include <iostream>
#include <vector>


class RasterBuffer
{
private:

	GLuint VAO, IBO;
	std::vector<GLuint> VBOlist;

public:

	void SetVAO(GLuint vao);
	GLuint GetVAO();

	GLuint AddVBO(GLuint vbo);
	GLuint GetVBO(GLuint id);

	void SetIBO(GLuint ibo);
	GLuint GetIBO();

	void BindVertexBuffer();
	void DeleteBuffer();
};



class ShaderStorageBuffer
{

private:

	GLuint BufferID, BufferIndex, dataRangeL,dataRangeR;

public:

	void SetBufferID(GLuint ibo);
	GLuint GetBufferID();

	void SetBufferIndex(GLuint ibo);
	GLuint GetBufferIndex();

	void SetDataRange(GLuint pointLeft,GLuint pointRight);
	GLuint GetPointerLeft();
	GLuint GetPointerRight();

	void BindBuffer();

	template<typename T>
	inline void ModifyBuffer(std::vector<T>& data) {

			ShaderStorageBuffer ssbo;

			GLuint ssboID,oldID;
			glGenBuffers(1, &ssboID);

			oldID = this->BufferID;

			glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboID);
			std::cout << glGetError() << std::endl;
			glBufferData(GL_SHADER_STORAGE_BUFFER, this->dataRangeR + (data.size()*sizeof(T)), NULL/* & data[0]*/, GL_DYNAMIC_DRAW);
			std::cout << glGetError() << std::endl;
			glCopyNamedBufferSubData(this->BufferID,ssboID,NULL,NULL, this->dataRangeR);
			std::cout << glGetError() << std::endl;
			glNamedBufferSubData(ssboID, this->dataRangeR, (data.size()*sizeof(T)), &data[0]);
			std::cout << glGetError() << std::endl;
			glBindBufferBase(GL_SHADER_STORAGE_BUFFER, this->BufferIndex, ssboID);
			std::cout << glGetError() << std::endl;
			glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
			std::cout << glGetError() << std::endl;

			this->dataRangeR = this->dataRangeR + (data.size() * sizeof(T));

			this->BufferID = ssboID;

			ssbo.SetBufferIndex(this->BufferIndex);

			glDeleteBuffers(1,&oldID);
	}

	/*template<typename T>
	inline void ModifyBufferInRange(std::vector<T>& data,int From, int To) {

		glBindBuffer(GL_SHADER_STORAGE_BUFFER, this->BufferID);
		std::cout << glGetError() << std::endl;
		glBufferData(GL_SHADER_STORAGE_BUFFER, this->dataRangeR + (data.size() * sizeof(T)), NULL/* & data[0]*///, GL_DYNAMIC_DRAW);
		/*std::cout << glGetError() << std::endl;
		glCopyNamedBufferSubData(this->BufferID, ssboID, NULL, NULL, this->dataRangeR);
		std::cout << glGetError() << std::endl;
		glNamedBufferSubData(ssboID, this->dataRangeR, (data.size() * sizeof(T)), &data[0]);
		std::cout << glGetError() << std::endl;
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, this->BufferIndex, this->BufferID);
		std::cout << glGetError() << std::endl;
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
		std::cout << glGetError() << std::endl;

		this->dataRangeR = this->dataRangeR + (data.size() * sizeof(T));

		this->BufferID = ssboID;

		ssbo.SetBufferIndex(this->BufferIndex);

		glDeleteBuffers(1, &oldID);
	}*/


	void DeleteBuffer();
};

class ShaderStorageManager {

private:
	static GLuint ssboCounter;
	static std::vector<ShaderStorageBuffer> storageBuffers;

public:
	static void IncreaseCounter();
	static GLuint GetCounter();
};

class UniformBufferManager {

private:
	static GLuint uboCounter;
	static std::vector<ShaderStorageBuffer> uniformBuffers;

public:
	static void IncreaseCounter();
	static GLuint GetCounter();
};
