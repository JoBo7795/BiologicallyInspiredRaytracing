#include "ObjectCreator.h"


GameObject ObjectCreator::GenerateSquare() {

	std::vector<float> vertices = {
		// vordere Fläche
		-1.0, -1.0,  1.0, 0.0f,1.0f,
		 1.0, -1.0,  1.0, 1.0f,1.0f,
		 1.0,  1.0,  1.0, 1.0f,0.0f,
		-1.0,  1.0,  1.0, 0.0f,0.0f,

		// hintere Fläche
		-1.0, -1.0, -1.0, 1.0f,1.0f,
		-1.0,  1.0, -1.0, 1.0f,0.0f,
		 1.0,  1.0, -1.0, 0.0f,0.0f,
		 1.0, -1.0, -1.0, 0.0f,1.0f,

		 // obere Fläche
		 -1.0,  1.0, -1.0, 0.0f,0.0f,
		 -1.0,  1.0,  1.0, 1.0f,0.0f,
		  1.0,  1.0,  1.0, 1.0f,1.0f,
		  1.0,  1.0, -1.0, 0.0f,1.0f,

		  // untere Fläche
		  -1.0, -1.0, -1.0, 0.0f,0.0f,
		   1.0, -1.0, -1.0, 1.0f,0.0f,
		   1.0, -1.0,  1.0, 1.0f,1.0f,
		  -1.0, -1.0,  1.0, 0.0f,1.0f,

		  // rechte Fläche
		   1.0, -1.0, -1.0, 0.0f,1.0f,
		   1.0,  1.0, -1.0, 0.0f,0.0f,
		   1.0,  1.0,  1.0, 1.0f,0.0f,
		   1.0, -1.0,  1.0, 1.0f,1.0f,

		   // linke Fläche
		   -1.0, -1.0, -1.0, 0.0f,1.0f,
		   -1.0, -1.0,  1.0, 1.0f,1.0f,
		   -1.0,  1.0,  1.0, 1.0f,0.0f,
		   -1.0,  1.0, -1.0, 0.0f,0.0f,
	};

	std::vector<unsigned int> indices = {
	0,  1,  2,      0,  2,  3,    // vorne
	4,  5,  6,      4,  6,  7,    // hinten
	8,  9,  10,     8,  10, 11,   // oben
	12, 13, 14,     12, 14, 15,   // unten
	16, 17, 18,     16, 18, 19,   // rechts
	20, 21, 22,     20, 22, 23    // links
	};

	RasterBuffer buffer = CreateBuffer::CreateIndexedBuffer(&vertices[0], &indices[0], vertices.size() * sizeof(float), indices.size() * sizeof(unsigned int), VerticeMap(3, 0, 0, 2));

	ShapeData shapeData;
	shapeData.vertices = vertices;
	shapeData.indices = indices;
	shapeData.vertexSize = vertices.size() * sizeof(float);
	shapeData.indexSize = indices.size() * sizeof(unsigned int);

	DrawObject drawObject = DrawObject(buffer, shapeData);

	GameObject Square = GameObject(DrawObjectManager::RegisterDrawObject(drawObject),0);//GameObject(buffer, vertices, indices, vertices.size(), indices.size());

	//Square.SetDrawObjectID();

	Square.SetBoundingBox(Mathematics::GenerateBoundingBox(vertices));

	std::vector<glm::vec4> tmpPos;

	int offset = 5;

	for (int i = 0; i < vertices.size() / offset; i++) {
		tmpPos.push_back(glm::vec4(vertices[i * offset], vertices[i * offset + 1], vertices[i * offset + 2], 1.0f));
	}

	//Square.SetTransformPosition(tmpPos);
	
	return Square;
}

GameObject ObjectCreator::CustomObject(GLuint drawObjectID, GLuint textureID) {
	
	GameObject Obj = GameObject(drawObjectID,textureID);

	ShapeData shapeData = DrawObjectManager::GetShapeDataByID(drawObjectID);

	Obj.SetBoundingBox(shapeData.boundingBox);

	std::vector<glm::vec4> tmpPos;

	int offset = 5;

	for (int i = 0; i < shapeData.vertices.size() / offset; i++) {
		tmpPos.push_back(glm::vec4(shapeData.vertices[i * offset], shapeData.vertices[i * offset + 1], shapeData.vertices[i * offset + 2], 1.0f));
	}

	return Obj;
}

RTGameObject ObjectCreator::CustomRTObject(GLuint drawObjectID, GLuint textureID) {

	RTGameObject Obj = RTGameObject(drawObjectID, textureID);

	ShapeData shapeData = DrawObjectManager::GetShapeDataByID(drawObjectID);

	Obj.SetBoundingBox(shapeData.boundingBox);

	std::vector<glm::vec4> tmpPos;

	int offset = 5;

	for (int i = 0; i < shapeData.vertices.size() / offset; i++) {
		tmpPos.push_back(glm::vec4(shapeData.vertices[i * offset], shapeData.vertices[i * offset + 1], shapeData.vertices[i * offset + 2], 1.0f));
	}

	return Obj;
}

BoundBox ObjectCreator::GenerateClickPoint(glm::vec3 position) {

	std::vector<float> vertices = {
		position.x,position.y,position.z
	};

	return Mathematics::GenerateBoundingBox(vertices);
}

BoundBox ObjectCreator::GeneratePlayerBB(glm::vec3 position) {

	float size = 0.1f, height = .2f;

	std::vector<glm::vec4> vertices = {
		glm::vec4(position.x + size,position.y + height ,position.z + size,1.0f),
		glm::vec4(position.x - size,position.y + height ,position.z + size,1.0f),
		glm::vec4(position.x + size,position.y + height ,position.z - size,1.0f),
		glm::vec4(position.x - size,position.y + height ,position.z - size,1.0f),

		glm::vec4(position.x + size,position.y - height ,position.z + size,1.0f),
		glm::vec4(position.x - size,position.y - height ,position.z + size,1.0f),
		glm::vec4(position.x + size,position.y - height ,position.z - size,1.0f),
		glm::vec4(position.x - size,position.y - height ,position.z - size,1.0f),
	};

	return Mathematics::GenerateBoundingBox(vertices);
}

GameObject ObjectCreator::Generate2DSquare() {

	std::vector<float> vertices = {
		-1.0f, -1.0f, 0.0f, 1.0f,
		 1.0f, -1.0f, 1.0f, 1.0f,
		 1.0f,  1.0f, 1.0f, 0.0f,
		-1.0f,  1.0f, 0.0f, 0.0f,
	};

	std::vector<unsigned int> indices = {
		0,  1,  2,      0,  2,  3,    // vorne
	};	

	RasterBuffer buffer = CreateBuffer::CreateIndexedBuffer(&vertices[0], &indices[0], vertices.size() * sizeof(float), indices.size() * sizeof(unsigned int), VerticeMap(2, 0, 0, 2));
		
	ShapeData shapeData;
	shapeData.vertices = vertices;
	shapeData.indices = indices;
	shapeData.vertexSize = vertices.size() * sizeof(float);
	shapeData.indexSize = indices.size() * sizeof(unsigned int);

	DrawObject drawObject = DrawObject(buffer, shapeData);	

	GameObject Square = GameObject(DrawObjectManager::RegisterDrawObject(drawObject),0);	

	Square.SetIndiSize(shapeData.indexSize);

	return Square;
}

void ObjectCreator::SerializeToSSBO(GameObject& go, int objectSize) {
	auto drawO = DrawObjectManager::GetDrawObjectByID(go.GetDrawObjectID());

	std::vector<float> serializedObjects;
	serializedObjects.push_back(1.0f);
	serializedObjects.push_back(go.GetObjectID());
	serializedObjects.push_back(drawO.VertBufferRangeL / objectSize);
	serializedObjects.push_back(0.0f);

	serializedObjects.push_back(drawO.IndiBufferRangeL);
	serializedObjects.push_back(drawO.IndiBufferRangeR);
	serializedObjects.push_back(drawO.VertBufferRangeL);
	serializedObjects.push_back(drawO.VertBufferRangeR);


	auto translation = go.GetTranslation();

	// translation matrix
	serializedObjects.push_back(translation[0][0]);
	serializedObjects.push_back(translation[0][1]);
	serializedObjects.push_back(translation[0][2]);
	serializedObjects.push_back(translation[0][3]);

	serializedObjects.push_back(translation[1][0]);
	serializedObjects.push_back(translation[1][1]);
	serializedObjects.push_back(translation[1][2]);
	serializedObjects.push_back(translation[1][3]);

	serializedObjects.push_back(translation[2][0]);
	serializedObjects.push_back(translation[2][1]);
	serializedObjects.push_back(translation[2][2]);
	serializedObjects.push_back(translation[2][3]);

	serializedObjects.push_back(translation[3][0]);
	serializedObjects.push_back(translation[3][1]);
	serializedObjects.push_back(translation[3][2]);
	serializedObjects.push_back(translation[3][3]);


	if (!DrawObjectManager::drawObjectBufferInit) {
		std::cout << "before drawbuffer:" << glGetError() << std::endl;
		DrawObjectManager::drawObjectBuffer = CreateBuffer::GenerateShaderStorageBufferObject(&serializedObjects[0], sizeof(float) * serializedObjects.size());
		DrawObjectManager::drawObjectBuffer.SetDataRange(0, sizeof(float) * (serializedObjects.size()));
		DrawObjectManager::drawObjectBufferInit = true;
		std::cout << "after drawbuffer:" << glGetError() << std::endl;

	}
	else {
		DrawObjectManager::drawObjectBuffer.ModifyBuffer(serializedObjects);
		std::cout << "mod drawbuffer:" << glGetError() << std::endl;
	}
}

void ObjectCreator::SerializeToSSBO(std::vector<float>& sphereArr) {
	CreateBuffer::GenerateShaderStorageBufferObject(&sphereArr[0], sizeof(float) * sphereArr.size());
}

void ObjectCreator::SerializeToSSBO(std::vector<Sphere>& sphereArr) {

	std::vector<float> dataArr;
	int size = 0;
	for (Sphere& sphere : sphereArr) {
		size += sphere.floatSize;
		for (float val : sphere.serialize()) {			
			dataArr.push_back(val);
		}
	}

	CreateBuffer::GenerateShaderStorageBufferObject(&dataArr[0], sizeof(float) * size);
}

ShaderStorageBuffer ObjectCreator::SerializeToUBO(std::vector<float>& data) {
	return CreateBuffer::GenerateUniformBufferObject(&data[0], sizeof(float) * data.size());
}