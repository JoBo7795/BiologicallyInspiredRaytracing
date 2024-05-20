#pragma once

#ifndef __COLLISIONDETECTION_HEADER__
#define __COLLISIONDETECTION_HEADER__

#include <vector>
#include "../Objects/BBStruct.h"
#include "../../Dependencies/glm/glm/glm.hpp"
#include "../Renderer/Renderer.h"

#endif

#define EPSILON 0.0000001

namespace CollisionDetection
{

    struct hit_record {
        glm::vec3 p;
        glm::vec3 normal;
        double t;
        bool front_face;
        int mat_ptr;
        int arrIndex;
        int geoType;
        //Triangle tri;
    };

	bool AABBCollisionDetection(BoundBox* boxA, BoundBox* boxB);
	bool RayPlaneIntersection(glm::vec3& n, glm::vec3& l, glm::vec3& p0, glm::vec3& l0);
	bool RaySphereIntersection(glm::vec3& spherePosition, double sphereRadius, glm::vec3& rayPosition, glm::vec3& rayDirection, double t_min, double& t_max);
    bool RaySphereIntersectionSimple(glm::vec3& spherePosition, float sphereRadius, glm::vec3& rayPosition, glm::vec3& rayDirection);
    bool RayTriangleIntersect(
        glm::vec3& rayPosition, glm::vec3& rayDirection,
        double t_min,
        double& t_max,
        std::array<glm::vec3, 3>& position,
        glm::vec3& outIntersectionPoint);
};

