#include "Camera.h"

Camera::Camera(glm::vec3 position) {
	this->position = position;
	this->target = glm::vec3(0.0f, 0.0f, 0.0f);
	this->direction = glm::vec3(0.0f, 0.0f, -1.0f);
	this->right = glm::vec3(0.0f, 0.0f, 0.0f);
	this->up = glm::vec3(0.0f, 1.0f, 0.0f);
	this->view = glm::mat4(1.0f);
	this->zoom = 45.0f;
}

Camera::~Camera() {


}

glm::mat4 Camera::GetViewMat() {
	return this->view;
}

void Camera::SetViewMat(glm::mat4 viewMat) {
	this->view = viewMat;
}

glm::vec3 Camera::GetPosition() {
	return this->position;
}

void Camera::SetPosition(glm::vec3 position) {
	this->position = position;
}

glm::vec3 Camera::GetTarget() {
	return this->target;
}

void Camera::SetTarget(glm::vec3 target) {
	this->target = target;
}

glm::vec3 Camera::GetDirection() {
	return this->direction;
}

void Camera::SetDirection(glm::vec3 direction) {
	this->direction = direction;
}

glm::vec3 Camera::GetUp() {
	return this->up;
}

void Camera::SetUp(glm::vec3 up) {
	this->up = up;
}

float Camera::GetZoom() {
	return this->zoom;
}

void Camera::SetZoom(float zoom) {
	this->zoom = zoom;
}

glm::mat4 Camera::LookAt() {

	return glm::lookAt(this->position,this->position+this->direction,this->up);
}

void Camera::Serialize(std::vector<float>& serialized) {


	serialized.push_back(this->theta);
	serialized.push_back(this->h);
	serialized.push_back(this->viewport_height);
	serialized.push_back(this->viewport_width);

	serialized.push_back(this->fov);
	serialized.push_back(this->aspect_ratio);
	serialized.push_back(this->aperture);
	serialized.push_back(this->focus_dist);

	serialized.push_back(this->focal_length);
	serialized.push_back(this->lens_radius);
	serialized.push_back(0.0f);
	serialized.push_back(0.0f);

	serialized.push_back(this->position[0]);
	serialized.push_back(this->position[1]);
	serialized.push_back(this->position[2]);
	serialized.push_back(0.0f);

	serialized.push_back(this->horizontal[0]);
	serialized.push_back(this->horizontal[1]);
	serialized.push_back(this->horizontal[2]);
	serialized.push_back(0.0f);

	serialized.push_back(this->vertical[0]);
	serialized.push_back(this->vertical[1]);
	serialized.push_back(this->vertical[2]);
	serialized.push_back(0.0f);

	serialized.push_back(this->lower_left_corner[0]);
	serialized.push_back(this->lower_left_corner[1]);
	serialized.push_back(this->lower_left_corner[2]);
	serialized.push_back(0.0f);

	glm::vec3 view = this->GetPosition() + this->GetDirection();

	serialized.push_back(this->position[0]);
	serialized.push_back(this->position[1]);
	serialized.push_back(this->position[2]);
	serialized.push_back(1.0f);

	serialized.push_back(view[0]);
	serialized.push_back(view[1]);
	serialized.push_back(view[2]);
	serialized.push_back(1.0f);

	serialized.push_back(this->up[0]);
	serialized.push_back(this->up[1]);
	serialized.push_back(this->up[2]);
	serialized.push_back(1.0f);

	serialized.push_back(this->w[0]);
	serialized.push_back(this->w[1]);
	serialized.push_back(this->w[2]);
	serialized.push_back(1.0f);


	serialized.push_back(this->u[0]);
	serialized.push_back(this->u[1]);
	serialized.push_back(this->u[2]);
	serialized.push_back(0.0f);

	serialized.push_back(this->v[0]);
	serialized.push_back(this->v[1]);
	serialized.push_back(this->v[2]);
	serialized.push_back(0.0f);

}

// Repeatedly update for position and directional data
void Camera::Update() {	

	glm::vec3 view = GetPosition() + GetDirection();// camera->LookAt();
	//view = GetPosition() + GetDirection();// camera->LookAt();
	this->view = glm::translate(glm::mat4(1.0),GetPosition() + GetDirection());
	
	w = glm::normalize(GetPosition() - view);
	u = glm::normalize(glm::cross(GetUp(), w));
	v = glm::cross(w, u);

	horizontal = glm::vec3(focus_dist * viewport_width * u);
	vertical = glm::vec3(focus_dist * viewport_height * v);
	lower_left_corner = GetPosition() - horizontal / 2.0f - vertical / 2.0f - focus_dist * w;

	serData.clear();
	Serialize(serData);

	if (!uboInit) {
		this->ubo = ObjectCreator::SerializeToUBO(serData);
		uboInit = true;
	}
	else {
		CreateBuffer::SetBufferSubData(GL_UNIFORM_BUFFER, this->ubo.GetBufferID(), 0, sizeof(float) * serData.size(), serData[0]);
	}

}

void Camera::InitRTCamera() {

	glm::vec3 view = GetPosition() + GetDirection();
	theta = glm::radians(fov);
	h = glm::tan(float(theta / 2));
	viewport_height = 2.0 * h;
	viewport_width = aspect_ratio * viewport_height;
	focal_length = 1.0f;
	lens_radius = aperture / 2;

}