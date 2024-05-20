#pragma once

#include <vector>
#include "../../Dependencies/glm/glm/glm.hpp"

struct BoundBox {

	float xmax, xmin, ymax, ymin, zmax, zmin;
	std::vector<float> xList;
	std::vector<float> yList;
	std::vector<float> zList;
	glm::vec3 normals[6] = {
		glm::vec3(0.0f,0.0f,1.0f),
		glm::vec3(0.0f,0.0f,-1.0f),
		glm::vec3(0.0f,1.0f,0.0f),
		glm::vec3(0.0f,-1.0f,0.0f),
		glm::vec3(1.0f,0.0f,0.0f),
		glm::vec3(-1.0f,0.0f,0.0f)
	};
};