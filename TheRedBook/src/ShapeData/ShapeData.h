#pragma once

#include "GL/glew.h"
#include "../../Dependencies/glm/glm/glm.hpp"
#include "../Mathematics/Mathematics.h"
#include <iostream>
#include <vector>

struct ShapeData
{	
	std::vector<GLuint> indices;	
	std::vector<GLfloat> vertices;
	std::vector<glm::vec3> normals;
	std::vector<glm::vec3> geometryPoints;
	std::vector<glm::vec3> edges;
	std::vector<std::array<glm::vec3, 3>> meshVerts;

	float maxRad;

	GLuint indexSize;
	GLuint vertexSize;
	GLuint normalSize;
	GLuint edgeSize;

	GLuint indexCount;
	GLuint vertexCount;
	GLuint normalCount;
	GLuint edgeCount;

	// array to map positions, normals, colors and textures
	GLuint vertMapArr[4][1] = { {0},{0},{0},{0} };
	GLuint mapSize;
	
	BoundBox boundingBox;

	GLuint matPointer = 2;
	glm::vec4 albedo = glm::vec4(1.0,0.0,0.0,1.0);

	inline void GenerateMeshVerts() {
		for (int i = 0; i < indices.size(); i += 3) {

			std::array<glm::vec3, 3> tri;

			tri[0] = geometryPoints[indices[i]];

			tri[1] = geometryPoints[indices[i + 1]];

			tri[2] = geometryPoints[indices[i + 2]];

			meshVerts.push_back(tri);

		}

	}

	void CalcRadius() {
		float center = 0.0f;
		for (int i = 0; i < geometryPoints.size(); i++) {
			if (glm::distance(geometryPoints[i], glm::vec3(center)) > maxRad)
				maxRad = glm::distance(geometryPoints[i], glm::vec3(center));
		}
	}
};