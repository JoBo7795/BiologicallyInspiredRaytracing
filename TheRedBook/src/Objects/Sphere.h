#pragma once

#include "../../Dependencies/glm/glm/glm.hpp"
#include "../../Dependencies/glm/glm/gtc/matrix_transform.hpp"
#include "../../Dependencies/glm/glm/gtc/type_ptr.hpp"
#include "../../Dependencies/glm/glm/gtx/string_cast.hpp""

#include <memory>
#include <vector>

struct Sphere
{
	int ID;
	glm::vec3 position;
	float radius, fuzz, mat_ptr;
    glm::vec4 color;
    glm::vec4 albedo;

	int floatSize = 20;


	Sphere() {};
	Sphere(glm::vec3 in_position, float in_radius, float in_fuzz, float in_mat_ptr, glm::vec4 in_color, glm::vec4 in_albedo) {
		ID = -1;
		position = in_position;
		radius = in_radius;
		fuzz = in_fuzz;
		mat_ptr = in_mat_ptr;
		color = in_color;
		albedo = in_albedo;
	};

	std::vector<float> serialize() {
		std::vector<float> valArr;
		valArr.resize(20);

		valArr[0] = position[0];
		valArr[1] = position[1];
		valArr[2] = position[2];

		valArr[3] = radius;

		valArr[4] = color[0];
		valArr[5] = color[1];
		valArr[6] = color[2];
		valArr[7] = color[3];

		valArr[8] = albedo[0];
		valArr[9] = albedo[1];
		valArr[10] = albedo[2];
		valArr[11] = albedo[3];

		valArr[12] = fuzz;

		valArr[13] = mat_ptr;

		return valArr;
		/*sphere.position = glm::vec3(sphereArr[i], sphereArr[i + 1], sphereArr[i + 2]);
		sphere.radius = sphereArr[i + 3];
		sphere.color = glm::vec4(sphereArr[i + 4], sphereArr[i + 5], sphereArr[i + 6], sphereArr[i + 7]);
		sphere.albedo = glm::vec4(sphereArr[i + 8], sphereArr[i + 9], sphereArr[i + 10], sphereArr[i + 11]);
		sphere.fuzz = sphereArr[i + 12];
		sphere.mat_ptr = sphereArr[i + 13];*/
		
	}
};