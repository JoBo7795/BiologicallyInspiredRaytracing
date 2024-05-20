#include "CollisionDetection.h"

bool CollisionDetection::AABBCollisionDetection(BoundBox* boxA, BoundBox* boxB) {

    return (((boxA->xmax <= boxB->xmax && boxA->xmax >= boxB->xmin) || (boxA->xmin <= boxB->xmax && boxA->xmin >= boxB->xmin)) &&
        ((boxA->ymax <= boxB->ymax && boxA->ymax >= boxB->ymin) || (boxA->ymin <= boxB->ymax && boxA->ymin >= boxB->ymin)) &&
        ((boxA->zmax <= boxB->zmax && boxA->zmax >= boxB->zmin) || (boxA->zmin <= boxB->zmax && boxA->zmin >= boxB->zmin)));

}

//Checks Collision based on the worldposition of the Objects, ignores object with collision turned off
/*GLuint CollisionDetection::PlayerWorldObjectCollisionDetection(BoundBox* boxA) {

    glm::vec3 playerPos = Renderer::GetCamera()->GetPosition();

    //std::vector<GameObject*> collBlocks = WorldGenerator::GetSurroundingBlocks(playerPos);
    std::vector<ChunkObject*> collBlocks = WorldGenerator::GetSurroundingBlocks(playerPos);
    int collBlockSize = collBlocks.size();    

    for (int i = 0; i < collBlockSize; i++) {
        //boxB = &collBlocks[i]->GetBoundingBox();
        BoundBox boxB = Mathematics::GenerateBoundingBox(glm::vec3(glm::vec4(1.0f)*collBlocks[0]->transform), 0.1f);
        if (((boxA->xmax <= boxB.xmax && boxA->xmax >= boxB.xmin) || (boxA->xmin <= boxB.xmax && boxA->xmin >= boxB.xmin)) &&
            ((boxA->ymax <= boxB.ymax && boxA->ymax >= boxB.ymin) || (boxA->ymin <= boxB.ymax && boxA->ymin >= boxB.ymin)) &&
            ((boxA->zmax <= boxB.zmax && boxA->zmax >= boxB.zmin) || (boxA->zmin <= boxB.zmax && boxA->zmin >= boxB.zmin))) {
            if (collBlocks[i]->collision) {
                //std::cout << collBlocks[i]->GetCollision() << std::endl;
                return collBlocks[i]->ID;
            }
        }
    }    
    return -1;
}

//Checks Collision based on the world position of the objects, checks all objects
GLuint CollisionDetection::PlayerWorldObjectCollisionDetectionNoCheck(BoundBox* boxA) {

    glm::vec3 playerPos = Renderer::GetCamera()->GetPosition();

    //std::vector<GameObject*> collBlocks = WorldGenerator::GetSurroundingBlocks(playerPos, 7, 7, 7);
    std::vector<ChunkObject*> collBlocks = WorldGenerator::GetSurroundingBlocks(playerPos, 7, 7, 7);
    int collBlockSize = collBlocks.size();
    BoundBox* boxB;

    GLuint id;
    //BoundBox baseBB = Mathematics::GenerateBoundingBox(glm::mat4(1.0f));
    for (int i = 0; i < collBlockSize; i++) {
        BoundBox boxB = Mathematics::GenerateBoundingBox(glm::vec3(glm::vec4(1.0f) * collBlocks[0]->transform), 0.1f);
        //boxB = &collBlocks[i]->GetBoundingBox();
        if (((boxA->xmax <= boxB.xmax && boxA->xmax >= boxB.xmin) || (boxA->xmin <= boxB.xmax && boxA->xmin >= boxB.xmin)) &&
            ((boxA->ymax <= boxB.ymax && boxA->ymax >= boxB.ymin) || (boxA->ymin <= boxB.ymax && boxA->ymin >= boxB.ymin)) &&
            ((boxA->zmax <= boxB.zmax && boxA->zmax >= boxB.zmin) || (boxA->zmin <= boxB.zmax && boxA->zmin >= boxB.zmin))) {            
            return collBlocks[i]->ID;
        }
    }
    return -1;
}

//Checks Collision based on the world position of the objects, checks all objects
std::vector<GLuint>* CollisionDetection::PlayerWorldObjectCollisionDetectionList(BoundBox* boxA) {

    glm::vec3 playerPos = Renderer::GetCamera()->GetPosition();

    //std::vector<GameObject*> collBlocks = WorldGenerator::GetSurroundingBlocks(playerPos,7,7,7);
    std::vector<ChunkObject*> collBlocks = WorldGenerator::GetSurroundingBlocks(playerPos,7,7,7);
    int collBlockSize = collBlocks.size();
    //BoundBox* boxB;

    std::vector<GLuint>* idList = new std::vector<GLuint>;

    for (int i = 0; i < collBlockSize; i++) {
        //boxB = &collBlocks[i]->GetBoundingBox();
        BoundBox boxB = Mathematics::GenerateBoundingBox(glm::vec3(glm::vec4(1.0f) * collBlocks[0]->transform), 0.1f);
        if (((boxA->xmax <= boxB.xmax && boxA->xmax >= boxB.xmin) || (boxA->xmin <= boxB.xmax && boxA->xmin >= boxB.xmin)) &&
            ((boxA->ymax <= boxB.ymax && boxA->ymax >= boxB.ymin) || (boxA->ymin <= boxB.ymax && boxA->ymin >= boxB.ymin)) &&
            ((boxA->zmax <= boxB.zmax && boxA->zmax >= boxB.zmin) || (boxA->zmin <= boxB.zmax && boxA->zmin >= boxB.zmin))) {
                std::cout << "list ids: " << collBlocks[i]->ID << std::endl;
                idList->push_back(collBlocks[i]->ID);    
        }
    }
    return idList;
}

bool CollisionDetection::PlayerWorldCollisionDetection(BoundBox* boxA) {

    glm::vec3 playerPos = Renderer::GetCamera()->GetPosition();

    std::vector<ChunkObject*> collBlocks = WorldGenerator::GetSurroundingBlocks(playerPos);
    int collBlockSize = collBlocks.size();

    glm::vec4 boxCoords;

    for (int i = 0; i < collBlockSize; i++) {

        if (!collBlocks[i]->collision)
            continue;

        boxCoords = collBlocks[i]->transform * glm::vec4(1.0f);
 
        BoundBox boxB = Mathematics::GenerateBoundingBox(glm::vec3(boxCoords) - glm::vec3(0.1f,0.0f,0.1f), 0.2f);

        if (((boxA->xmax <= boxB.xmax && boxA->xmax >= boxB.xmin) || (boxA->xmin <= boxB.xmax && boxA->xmin >= boxB.xmin)) &&
            ((boxA->ymax <= boxB.ymax && boxA->ymax >= boxB.ymin) || (boxA->ymin <= boxB.ymax && boxA->ymin >= boxB.ymin)) &&
            ((boxA->zmax <= boxB.zmax && boxA->zmax >= boxB.zmin) || (boxA->zmin <= boxB.zmax && boxA->zmin >= boxB.zmin))) {
                return true;
        }
    }
    return false;
}*/

bool CollisionDetection::RayPlaneIntersection(glm::vec3& n, glm::vec3& l, glm::vec3& p0, glm::vec3& l0) {

    // assuming vectors are all normalized
    float denom = glm::dot(n, l);
    if (denom > 1e-6) {
        glm::vec3 p0l0 = p0 - l0;
        float t = glm::dot(p0l0, n) / denom;
        
        return (t >= 0);
    }

    return false;
}

bool CollisionDetection::RaySphereIntersection(glm::vec3& spherePosition, double sphereRadius, glm::vec3& rayPosition, glm::vec3& rayDirection, double t_min, double& t_max) {

        glm::vec3 oc = rayPosition - spherePosition;
        double a = glm::length(rayDirection) * glm::length(rayDirection);
        double half_b = glm::dot(oc, rayDirection);
        double c = glm::length(oc) * glm::length(oc) - sphereRadius * sphereRadius;

        double discriminant = half_b * half_b - a * c;
        if (discriminant < 0) return false;
        double sqrtd = glm::sqrt(discriminant);

        // Find the nearest root that lies in the acceptable range.
        double root = (-half_b - sqrtd) / a;
        if (root < t_min || t_max < root) {
            root = (-half_b + sqrtd) / a;
            if (root < t_min || t_max < root)
                return false;
        }

        t_min = root;

        return true;
    
}

bool CollisionDetection::RaySphereIntersectionSimple(glm::vec3& spherePosition, float sphereRadius, glm::vec3& rayPosition, glm::vec3& rayDirection) {

    glm::vec3 oc = rayPosition - spherePosition;
    double a = glm::length(rayDirection) * glm::length(rayDirection);
    double half_b = glm::dot(oc, rayDirection);
    double c = glm::length(oc) * glm::length(oc) - sphereRadius * sphereRadius;

    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return false;


    return true;
}

bool CollisionDetection::RayTriangleIntersect(
    glm::vec3& rayPosition, glm::vec3& rayDirection,
    double t_min,
    double& t_max,
    std::array<glm::vec3,3>& position,
    //glm::vec3 normal,
    glm::vec3& outIntersectionPoint)
{

    glm::vec3 vertex0 = position[0];
    glm::vec3 vertex1 = position[1];
    glm::vec3 vertex2 = position[2];

    glm::vec3 edge1, edge2, h, s, q;

    double a, f, u, v;

    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;

    h = glm::cross(rayDirection, edge2);
    a = glm::dot(edge1, h);

    if (a > -EPSILON && a < EPSILON)
        return false;    // This ray is parallel to this triangle.

    f = 1.0 / a;
    s = rayPosition - vertex0;
    u = f * glm::dot(s, h);
    if (u < 0.0 || u > 1.0)
        return false;
    q = glm::cross(s, edge1);
    v = f * glm::dot(rayDirection, q);
    if (v < 0.0 || u + v > 1.0)
        return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    float t = f * glm::dot(edge2, q);

    if (t < t_min || t_max < t)
        return false;


    //rec.t = t;
    //rec.p = ray_at(r.origin, r.direction, rec.t);

    //vec3 outward_normal = tri.normal.xyz;//normalize(cross(edge1, edge2));//(rec.p - sphere.positionAndradius.xyz) / sphere.positionAndradius.w;
    //glm::vec3 outward_normal = normalize(cross(edge1, edge2));//(rec.p - sphere.positionAndradius.xyz) / sphere.positionAndradius.w;
    //rec.front_face = check_front_face(r, outward_normal);
    //rec.normal = set_face_normal(r, outward_normal);
    //rec.mat_ptr = int(tri.fuzzAndmat_ptr.x);

    if (t > EPSILON) // ray intersection
    {
        outIntersectionPoint = rayPosition + rayDirection * t;
        t_max = t;
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return false;
}