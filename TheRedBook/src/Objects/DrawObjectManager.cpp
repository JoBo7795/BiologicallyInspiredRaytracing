#include "DrawObjectManager.h"

std::vector<DrawObject> DrawObjectManager::ObjectSlotQueue;
std::vector<Sphere> DrawObjectManager::SphereQueue;
std::vector<GLuint> DrawObjectManager::BoundObjects;
GLuint DrawObjectManager::boundID = -1;
GLuint DrawObjectManager::currentBound = 0;

bool DrawObjectManager::SphereStorageBufferInit = false,
DrawObjectManager::geoBufferInit = false,
DrawObjectManager::indiceBufferInit = false,
DrawObjectManager::floatListBufferInit = false,
DrawObjectManager::primArrBufferInit = false,
DrawObjectManager::drawObjectBufferInit = false;

ShaderStorageBuffer DrawObjectManager::SphereStorageBuffer,
DrawObjectManager::geoBuffer,
DrawObjectManager::indiceBuffer,
DrawObjectManager::floatListBuffer,
DrawObjectManager::primArrBuffer,
DrawObjectManager::drawObjectBuffer;

int DrawObjectManager::RegisterDrawObject(ShapeData shapeData) {

	VerticeMap vertMap = VerticeMap(shapeData.vertMapArr[0][0], shapeData.vertMapArr[0][1], shapeData.vertMapArr[0][2], shapeData.vertMapArr[0][3]);

	RasterBuffer buffer = CreateBuffer::CreateIndexedBuffer(&shapeData.vertices[0], &shapeData.indices[0], shapeData.vertices.size() * sizeof(float), shapeData.indices.size() * sizeof(unsigned int), vertMap);

	DrawObject drawObject(buffer, shapeData);
	drawObject.id = ObjectSlotQueue.size();
	ObjectSlotQueue.push_back(drawObject);
	return ObjectSlotQueue.size() - 1;
}


int DrawObjectManager::RegisterRTDrawObject(ShapeData shapeData) {
	std::cout << "begin func:" << glGetError() << std::endl;
	VerticeMap vertMap = VerticeMap(shapeData.vertMapArr[0][0], shapeData.vertMapArr[0][1], shapeData.vertMapArr[0][2], shapeData.vertMapArr[0][3]);

	std::vector<float> plyGeo;
	float center = 0.0f, maxRad = 0.0f;

	shapeData.CalcRadius();
	maxRad = shapeData.maxRad;
	
	shapeData.GenerateMeshVerts();

	for (int i = 0; i < shapeData.geometryPoints.size(); i++) {

		plyGeo.push_back(shapeData.geometryPoints[i].x);
		plyGeo.push_back(shapeData.geometryPoints[i].y);
		plyGeo.push_back(shapeData.geometryPoints[i].z);
		plyGeo.push_back(1.0f);

		plyGeo.push_back(shapeData.normals[i].x);
		plyGeo.push_back(shapeData.normals[i].y);
		plyGeo.push_back(shapeData.normals[i].z);
		plyGeo.push_back(1.0f);

		//plyGeo.push_back(2500.0f);
		//plyGeo.push_back(0.2);
		//plyGeo.push_back(0.5);
		//plyGeo.push_back(1.0f);

		// plyGeo.push_back(shapeData.albedo.x);
		// plyGeo.push_back(shapeData.albedo.y);
		// plyGeo.push_back(shapeData.albedo.z);
		// plyGeo.push_back(shapeData.albedo.w);

		plyGeo.push_back(shapeData.vertices[i*vertMap.stride + (vertMap.positionSize+vertMap.normalSize)]);
		plyGeo.push_back(shapeData.vertices[i*vertMap.stride + (vertMap.positionSize+vertMap.normalSize)+1]);
		plyGeo.push_back(shapeData.vertices[i*vertMap.stride + (vertMap.positionSize+vertMap.normalSize)+2]);
		plyGeo.push_back(1.0f);

		plyGeo.push_back(0.0);
		//plyGeo.push_back(2.0);
		plyGeo.push_back(shapeData.matPointer);
		plyGeo.push_back(center); // center
		plyGeo.push_back(maxRad); // bbCircleRadius

	}

	int root, depth = 5;
	auto primitiveList = shapeData.meshVerts;
	auto list = Mathematics::CreateKDTree(primitiveList, 5, root);
	Mathematics::SortPrimitivesToNodes(list, primitiveList, root, depth);
	auto primArr = Mathematics::CreatePrimitiveArrayList(list);
	auto floatList = Mathematics::NodeToFloatList(list);
	int vL, vR, iL, iR;
	std::cout << "before geo:" << glGetError() << std::endl;

	if (!geoBufferInit) {
		geoBuffer = CreateBuffer::GenerateShaderStorageBufferObject(&plyGeo[0], sizeof(float) * plyGeo.size());
		geoBuffer.SetDataRange(0, sizeof(float) * (plyGeo.size()));
		vL = 0;
		vR = plyGeo.size()-1;
		geoBufferInit = true;
		std::cout << "before geo:" << glGetError() << std::endl;
	}
	else {
		vL = ((geoBuffer.GetPointerRight())  /sizeof(float));
		vR = ((geoBuffer.GetPointerRight()) / sizeof(float)) + plyGeo.size()-1;
		geoBuffer.ModifyBuffer(plyGeo);
		std::cout << "mod geo:" << glGetError() << std::endl;
	}

	if (!indiceBufferInit) {
		indiceBuffer = CreateBuffer::GenerateShaderStorageBufferObject(&shapeData.indices[0], sizeof(GLuint) * shapeData.indices.size());
		indiceBuffer.SetDataRange(0, sizeof(GLuint) * (shapeData.indices.size()));
		iL = 0;
		iR = shapeData.indices.size()-1;
		indiceBufferInit = true;
		std::cout << "before geo:" << glGetError() << std::endl;

	}
	else {
		iL = (indiceBuffer.GetPointerRight() /sizeof(GLuint) );
		iR = (indiceBuffer.GetPointerRight() / sizeof(GLuint)) + shapeData.indices.size()-1;
		indiceBuffer.ModifyBuffer(shapeData.indices);
		std::cout << "indice mod:" << glGetError() << std::endl;
	}

	if (!floatListBufferInit) {
		floatListBuffer = CreateBuffer::GenerateShaderStorageBufferObject(&floatList[0], sizeof(float) * floatList.size());
		floatListBuffer.SetDataRange(0, sizeof(float) * floatList.size());
		floatListBufferInit = true;
		std::cout << "before geo:" << glGetError() << std::endl;
	}
	else {
		floatListBuffer.ModifyBuffer(floatList);
		std::cout << "float mod:" << glGetError() << std::endl;
	}

	if (!primArrBufferInit) {
		primArrBuffer = CreateBuffer::GenerateShaderStorageBufferObject(&primArr[0], sizeof(int) * primArr.size());
		primArrBuffer.SetDataRange(0, sizeof(int) * primArr.size());
		primArrBufferInit = true;
		std::cout << "before geo:" << glGetError() << std::endl;
	}
	else {
		primArrBuffer.ModifyBuffer(primArr);
		std::cout << "mod prim:" << glGetError() << std::endl;
	}

	DrawObject drawObject(geoBuffer, shapeData);
	drawObject.VertBufferRangeL = vL;
	drawObject.VertBufferRangeR = vR;
	drawObject.IndiBufferRangeL = iL;
	drawObject.IndiBufferRangeR = iR;
	drawObject.id = ObjectSlotQueue.size();
	ObjectSlotQueue.push_back(drawObject);
	// type = TRIANGLE = 1
	
		
	return ObjectSlotQueue.size() - 1;
}




int DrawObjectManager::RegisterRTSphere(Sphere sphere, std::vector<float> sphereArr) {
	sphere.position = glm::vec3(sphereArr[0], sphereArr[1], sphereArr[2]);
	sphere.radius = sphereArr[3];
	sphere.color = glm::vec4(sphereArr[4], sphereArr[5], sphereArr[6], sphereArr[7]);
	sphere.albedo = glm::vec4(sphereArr[8], sphereArr[9], sphereArr[10], sphereArr[11]);
	sphere.fuzz = sphereArr[12];
	sphere.mat_ptr = sphereArr[13];

	SphereQueue.push_back(sphere);
	SphereQueue.back().ID = SphereQueue.size() - 1;
	return SphereQueue.back().ID;
}

int DrawObjectManager::RegisterRTSphere(Sphere sphere) {

	SphereQueue.push_back(sphere);
	SphereQueue.back().ID = SphereQueue.size()-1;
	return SphereQueue.back().ID;
}

int DrawObjectManager::RegisterDrawObject(DrawObject drawObject) {	

	ObjectSlotQueue.push_back(drawObject);
	return ObjectSlotQueue.size() - 1;
}

void DrawObjectManager::UseDrawObject(GLuint drawObjectID) {

	for (int i = 0; i < ObjectSlotQueue.size(); i++) {
		if (i == drawObjectID && drawObjectID != boundID) {
			ObjectSlotQueue[drawObjectID].buffer.BindVertexBuffer();
			boundID = drawObjectID;
			currentBound = drawObjectID;
		}
	}
}

ShapeData DrawObjectManager::GetShapeDataByID(GLuint drawObjectID) {

	for (int i = 0; i < ObjectSlotQueue.size(); i++) {
		if (i == drawObjectID ) {
			return ObjectSlotQueue[drawObjectID].shapeData;
		}
	}
}

DrawObject& DrawObjectManager::GetDrawObjectByID(GLuint drawObjectID) {

	for (int i = 0; i < ObjectSlotQueue.size(); i++) {
		if (i == drawObjectID) {
			return ObjectSlotQueue[drawObjectID];
		}
	}
}

DrawObject& DrawObjectManager::GetCurrentBound() {
	return ObjectSlotQueue[currentBound];
}

std::vector<Sphere>& DrawObjectManager::GetSphereQueue() {
	return SphereQueue;
}

std::vector<DrawObject>& DrawObjectManager::GetObjectSlotQueue() {
	return ObjectSlotQueue;
}