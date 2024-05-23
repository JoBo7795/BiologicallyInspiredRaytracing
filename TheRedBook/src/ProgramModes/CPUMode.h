#pragma once

#include <iostream>

#include "../Definitions.h"
#include "../Objects/GameObjectManager.h"
#include "../DataReader/DataReader.h"
#include "../Camera/Camera.h"
#include "ProgramParams.h"
#include "../Objects/LensData.h"
#include "../DataReader/DataWriter.h"
#include <fstream>
#include <map>

#define EPSILON 0.0000001
#define CHANNEL_NUM 3




struct Ray {
	glm::vec3 origin, direction;
	glm::vec4 color;
	int currDepth;
	bool inLense = false;
    int id = -1;
    glm::vec3 hit_normal_dir;
    glm::vec3 hit_normal_orig;
    GLuint rgb;
    float brightness;
};

struct Triangle {

    glm::vec4 position[3];
    glm::vec4 normal;
    glm::vec4 albedo;
    glm::vec4 fuzzAndmat_ptr;

};

struct hit_record {

    glm::vec3 p;
    glm::vec3 normal;
    double t;
    bool front_face;
    int mat_ptr;
    int arrIndex;
    int geoType;
    Triangle tri;

};

struct BeginEnd {
    glm::vec3 p1, p2;
    bool inLense;
    int numHits;
    glm::vec3 hit_normal_orig;
    glm::vec3 hit_normal_dir;
    int rgb;
};


class CPUMode
{
public:

    static LensData lens;
    static Camera debugCam;
    static std::ofstream dataText;
	static int depth;
    static const float infinity;
    static const float M_PI;
    static float global_refract_index;
    static uint8_t* pixels;
    static uint8_t* fBuffer;
    static glm::vec3 lowerLeftCorner, horizontal, vertical, lookFrom;
    static glm::vec3 initlowerLeftCorner, initHorizontal, initVertical, initCamPos;
    static std::vector<glm::vec3> circlePoints;
    static int numPoints;
    static std::map<float,float> lambdaVal, lambdaPrimeVal;

    static void Render();
    static void InitCPUMode(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir);
    static void SetDebugParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir);

    static glm::vec4 ray_color(Ray& r, int depth);
    static bool CameraGetRay(float s, float t, Ray& ray, int o, bool lp);
    static bool camera_get_debug_ray(float s, float t, Ray& ray, int o, bool lp);
    static glm::vec4 ray_color2(Ray& r, bool& hit, bool& no_reflect, bool& no_shadow);
    static bool WorldHit(float t_min, Ray& r, float t_max, hit_record& rec);
    static glm::vec3 ray_at(glm::vec3 origin, glm::vec3 direction, float t);
    static glm::vec3 refract_s(glm::vec3& uv, glm::vec3 n, float etai_over_etat);
    static bool SphereHit(float t_min, Ray& r, float t_max, hit_record& rec, Sphere& sphere);
    static bool check_front_face(const Ray r, const glm::vec3 outward_normal);
    static glm::vec3 set_face_normal(const Ray r, const glm::vec3 outward_normal);
    static bool RayTriangleIntersect(Ray r, hit_record& rec, float t_min, float t_max, Triangle tri, glm::vec3& outIntersectionPoint);
    static glm::vec3 triInterpolNormal(glm::vec3 tri0, glm::vec3 tri1, glm::vec3 tri2, glm::vec3 normal0, glm::vec3 normal1, glm::vec3 normal2, glm::vec3 vertex);
    static bool lambertian_scatter(Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Sphere& sphere);
    static bool lambertian_scatter(Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Triangle& tri);
    static bool plain_color(Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Sphere& sphere);
    static bool plain_color(Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Triangle& tri);
    static bool metal_scatter(Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Sphere& sphere);
    static bool metal_scatter(Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Triangle& tri);
    static bool DielectricScatter(Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Sphere& sphere);
    static bool DielectricScatter(Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Triangle& tri);

    static bool IntersectRaySphere(glm::vec3 rayOrigin, glm::vec3 rayDirection, glm::vec3 sphereOrigin, double sphereRadius, glm::vec3& iPoint, double& t);
    static bool IntersectRayLens(glm::vec3 rayOrigin, glm::vec3 rayDirection, LensData lens, glm::vec3& normal, double& t, bool inLense);

    static glm::vec3 random_unit_vector();
    static glm::vec3 random_in_unit_sphere();
    static glm::vec3 randomVec3(float min, float max);
    static glm::vec3 random_pcg3d(glm::uvec3 v);
    static void addCirclePoints(glm::vec3 center, float radius, int numCircles, int numPointsPerCircle, glm::vec3 normal);
    static glm::vec3 rotatePoint(glm::vec3 point, glm::vec3 center, glm::vec3 axis, float angle);
    static bool rayIntersectsPlane(glm::vec3 rayOrigin, glm::vec3 rayDirection, glm::vec3 planeOrigin, glm::vec3 planeNormal, glm::vec3 planeU, glm::vec3 planeV);
    static bool drawColoredRays(Ray& r);
    static bool intersectionRayCylinder(glm::vec3 rayOrigin, glm::vec3 rayDir, glm::vec3 cylinderStart, glm::vec3 cylinderEnd, float cylinderRadius, float& t);
    static bool intersectRayPlane(glm::vec3 rayOrigin, glm::vec3 rayDir, glm::vec3 planePoint, glm::vec3 planeNormal, float& t);
    static std::vector<BeginEnd> ComputeDebugRays(glm::vec3 imagePlanePos);
    static bool CheckDebugRayHit(std::vector<BeginEnd>& rayArr, Ray& r);

    static double logInterpolation(double x1, double y1, double x2, double y2, double x);

    static std::ofstream OpenData();
    static float lerp(float x, float y, float t);

};

