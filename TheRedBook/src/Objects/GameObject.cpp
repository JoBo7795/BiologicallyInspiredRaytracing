#include "GameObject.h"

GameObject::GameObject(){}
GameObject::~GameObject() {}

GameObject::GameObject(GLuint drawObjectID,GLuint textureID) {
	this->DrawObjectID = drawObjectID;	
	this->textureIDs.push_back(textureID);

	this->scale = glm::mat4(1.0f);
	this->translation = glm::mat4(1.0f);
	this->rotation = glm::mat4(1.0f);
	this->collision = true;
	this->color = glm::vec3(1.0f,1.0f,1.0f);
}

// Constructor to initialize Player
GameObject::GameObject(glm::mat4 scale, glm::mat4 translation, glm::mat4 rotation) {	

	this->scale = scale;
	this->translation = translation;
	this->rotation = rotation;
	this->color = glm::vec3(0.0f, 0.0, 1.0f);
	this->collision = true;
}

GLuint GameObject::GetObjectID() {
	return this->objectID;
}

void GameObject::SetObjectID(GLuint objectID) {
	this->objectID = objectID;
}

GLuint GameObject::GetDrawObjectID() {
	return this->DrawObjectID;
}

void GameObject::SetDrawObjectID(GLuint drawObjectID) {
	this->DrawObjectID = drawObjectID;
}

GLuint GameObject::GetChunkObjectID() {
	return this->chunkObjectID;
}

void GameObject::SetChunkObjectID(GLuint objectID) {
	this->chunkObjectID = objectID;
}

void GameObject::SetBuffer(RasterBuffer buffer) {
	//this->buffer = buffer;
}

unsigned int GameObject::GetVertiSize() {
	return this->vertSize * sizeof(float);
}

glm::vec3 GameObject::GetColor(){
	return this->color;
}

void GameObject::SetColor(glm::vec3 color) {
	this->color = color;
}

bool GameObject::GetTextureSet() {
	return this->textureSet;
}

void GameObject::SetTextureSet(bool set) {
	this->textureSet = set;
}

void GameObject::SetVertiSize(unsigned int size) {
	this->vertSize = size;
}

unsigned int GameObject::GetIndiSize(){
	return this->indiSize * sizeof(unsigned int);
}

void GameObject::SetIndiSize(unsigned int size){
	this->indiSize = size;
}

glm::mat4 GameObject::GetScale() {
	return this->scale;
}

void GameObject::SetScale(glm::mat4 scale) {
	this->scale = scale;
}

glm::mat4 GameObject::GetRotation() {
	return this->rotation;
}

void GameObject::SetRotation(glm::mat4 rotation) {
	this->rotation = rotation;
}

glm::mat4 GameObject::GetTranslation() {
	return this->translation;
}

void GameObject::SetTranslation(glm::mat4 translation) {
	this->translation = translation;
}

/*calculate the resulting transformation position of the object from translation, rotation and it's scale*/
void GameObject::TransformObject() {
		
	ShapeData sd = DrawObjectManager::GetShapeDataByID(DrawObjectID);
	std::vector<glm::vec4> transformVertList;

	for (int i = 0; i < sd.vertexCount;i++) {		
		transformVertList.push_back( this->translation * this->rotation * this->scale * glm::vec4(sd.vertices[i], sd.vertices[i + 1], sd.vertices[i + 2], 1.0f));
	}

	this->boundingbox = Mathematics::GenerateBoundingBox(transformVertList);
}

void GameObject::AddTexture(Texture texture) {
	//this->textures.push_back(texture);
}

std::vector<GLuint> GameObject::GetTextureIDs() {
	return this->textureIDs;
}

void GameObject::AddTextureID(GLuint texID) {
	this->textureIDs.push_back(texID);
}

BoundBox& GameObject::GetBoundingBox() {
	return this->boundingbox;
}

void GameObject::SetBoundingBox(BoundBox bbox) {
	this->boundingbox = bbox;
}

bool GameObject::GetCollision() {
	return this->collision;
}

void GameObject::SetCollision(bool val) {
	this->collision = val;
}

RTGameObject::RTGameObject(GLuint drawObjectID, GLuint textureID) {
	this->drawObjectID = drawObjectID;
	this->textureIDs.push_back(textureID);
}

RTGameObject::~RTGameObject() {};

void RTGameObject::SetBoundingBox(BoundBox bbox) {
	this->boundingbox = bbox;
}