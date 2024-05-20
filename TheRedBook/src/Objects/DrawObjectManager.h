#pragma once

#include <vector>
#include "DrawObjects.h"
#include "../Buffer/CreateBuffer.h"
#include "Sphere.h"


class DrawObjectManager
{
private:

	static std::vector<DrawObject> ObjectSlotQueue;
	static std::vector<Sphere> SphereQueue;
	static std::vector<GLuint> BoundObjects;
	static GLuint boundID;
	static GLuint currentBound;
	static bool SphereStorageBufferInit, geoBufferInit, indiceBufferInit, floatListBufferInit, primArrBufferInit;
	static ShaderStorageBuffer SphereStorageBuffer, geoBuffer, indiceBuffer, floatListBuffer, primArrBuffer;

	

public:

	static int RegisterDrawObject(ShapeData shapeData);
	static int RegisterRTDrawObject(ShapeData shapeData);
	static int RegisterRTSphere(Sphere sphere);
	static int RegisterRTSphere(Sphere sphere, std::vector<float> sphereArr);
	static int RegisterDrawObject(DrawObject drawObject);
	static void UseDrawObject(GLuint drawObjectID);
	static ShapeData GetShapeDataByID(GLuint drawObjectID);
	static DrawObject& GetDrawObjectByID(GLuint drawObjectID);
	static DrawObject& GetCurrentBound();
	static std::vector<Sphere>& GetSphereQueue();
	static std::vector<DrawObject>& GetObjectSlotQueue();
	

	static ShaderStorageBuffer drawObjectBuffer;
	static bool drawObjectBufferInit;

};

