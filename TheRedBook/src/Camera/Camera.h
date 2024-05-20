#pragma once

#include "../../Dependencies/glm/glm/glm.hpp"
#include "../../Dependencies/glm/glm/gtc/matrix_transform.hpp"
#include "../../Dependencies/glm/glm/gtc/type_ptr.hpp"
#include "../Buffer/CreateBuffer.h"
#include "../Objects/ObjectCreator.h"

#include <vector>

class Camera
{
private:
	glm::vec3 position, target, direction, right, front, up;
	glm::mat4 view;
	

public:

	float zoom, fov = 20, aspect_ratio = 16.0 / 9.0, aperture = 0.1f, focus_dist = 1.0f, imageWidth = 1920, imageHeight = 1080,
		theta,h, viewport_height, viewport_width, focal_length, lens_radius;
	glm::vec3 w, u, v, horizontal, vertical, lower_left_corner;
	ShaderStorageBuffer ubo;
	bool uboInit = false, updateData = true;
	std::vector<float> serData;
	
	Camera() {};
	Camera(glm::vec3 position);
	~Camera();

	glm::mat4 GetViewMat();
	void SetViewMat(glm::mat4 viewMat);

	glm::vec3 GetPosition();
	void SetPosition(glm::vec3 position);

	glm::vec3 GetTarget();
	void SetTarget(glm::vec3 target);

	glm::vec3 GetDirection();
	void SetDirection(glm::vec3 direction);

	glm::vec3 GetUp();
	void SetUp(glm::vec3 up);

	float GetZoom();
	void SetZoom(float zoom);

	glm::mat4 LookAt();

	double rawXposMouse, rawYposMouse;

	void Serialize(std::vector<float>& serialized);
	void Update();
	void InitRTCamera();

};

