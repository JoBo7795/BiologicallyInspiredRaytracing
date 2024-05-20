#include "Buffer.h"

void RasterBuffer::SetVAO(unsigned int vao) {

	this->VAO = vao;
}

unsigned int RasterBuffer::GetVAO() {

	return this->VAO;
}

GLuint RasterBuffer::AddVBO(GLuint vbo) {
	this->VBOlist.push_back(vbo);
	return this->VBOlist.size() - 1;
}

GLuint RasterBuffer::GetVBO(GLuint id) {

	return this->VBOlist[id];
}

void RasterBuffer::SetIBO(unsigned int ibo) {

	this->IBO = ibo;
}

unsigned int RasterBuffer::GetIBO() {

	return this->IBO;
}

void RasterBuffer::BindVertexBuffer() {
	glBindVertexArray(this->VAO);
}

void RasterBuffer::DeleteBuffer() {

	glDeleteVertexArrays(1, &this->VAO);
	glDeleteBuffers(VBOlist.size(), &this->VBOlist[0]);
	glDeleteBuffers(1, &this->IBO);
}

GLuint ShaderStorageManager::ssboCounter = 0;

void ShaderStorageManager::IncreaseCounter() {
	ssboCounter++;
}

GLuint ShaderStorageManager::GetCounter() {
	return ssboCounter;
}

void ShaderStorageBuffer::SetBufferID(GLuint bufferID) {
	this->BufferID = bufferID;
}

GLuint ShaderStorageBuffer::GetBufferID() {
	return this->BufferID;
}

void ShaderStorageBuffer::SetBufferIndex(GLuint bufferIndex) {

	this->BufferIndex = bufferIndex;
}
GLuint ShaderStorageBuffer::GetBufferIndex() {
	return this->BufferIndex;
}

void ShaderStorageBuffer::SetDataRange(GLuint pointLeft, GLuint pointRight) {
	this->dataRangeL = pointLeft;
	this->dataRangeR = pointRight;
}

GLuint ShaderStorageBuffer::GetPointerLeft() {
	return this->dataRangeL;
}

GLuint ShaderStorageBuffer::GetPointerRight() {
	return this->dataRangeR;
}

void ShaderStorageBuffer::BindBuffer(){
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, this->BufferID);
}

/*template<typename T>
void ShaderStorageBuffer::ModifyBuffer(std::vector<T>& data) {
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, this->BufferID);
	glNamedBufferSubData(this->BufferID,this->dataRangeR,data.size()*sizeof(T),&data[0]);
}*/

void ShaderStorageBuffer::DeleteBuffer(){
	glDeleteBuffers(1,&this->BufferID);
}

GLuint UniformBufferManager::uboCounter = 0;
void UniformBufferManager::IncreaseCounter() {
	uboCounter++;	
}

GLuint UniformBufferManager::GetCounter() {
	return uboCounter;
}
