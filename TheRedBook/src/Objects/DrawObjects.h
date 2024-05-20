#pragma once

#include "../ShapeData/ShapeData.h"
#include "../Buffer/Buffer.h"

struct DrawObject
{
	RasterBuffer buffer;

	ShaderStorageBuffer RTBuffer;
	// Range for SSBO
	GLuint VertBufferRangeL, VertBufferRangeR;
	GLuint IndiBufferRangeL, IndiBufferRangeR;
	// type eg. TRIANGLE = 1, SPHERE = 0
	GLuint type;

	ShapeData shapeData;
	int id;
	bool isRT = false;


	DrawObject(){}

	DrawObject(ShaderStorageBuffer buffer, ShapeData shapeData) {
		this->RTBuffer = buffer;
		this->shapeData = shapeData;
		this->isRT = true;
	}

	DrawObject(RasterBuffer buffer, ShapeData shapeData) {
		this->buffer = buffer;
		this->shapeData = shapeData;		
	}


};

