#version 460


uniform vec3 view;
uniform vec3 viewPos;
uniform vec3 up;
uniform vec3 lightPos;
uniform vec2 factXY;
uniform mat4 lenseTransform;
uniform vec2 widthHeight;
uniform vec3 initCamPos;

uniform vec3 initlowerLeftCorner;
uniform vec3 initHorizontal;
uniform vec3 initVertical;
//float a[5] = float[](, 4.2, 5.0, 5.2, 1.1);

//uniform float 

/*struct Camera {

    double theta;
    double h;
    double viewport_height;
    double viewport_width;

    double fov;
    double aspect_ratio;
    double aperture;
    double focus_dist;

    //        auto aspect_ratio = 16.0 / 9.0;

    double focal_length;

    vec3 origin;
    vec3 horizontal;
    vec3 vertical;
    vec3 lower_left_corner;

    double lens_radius;

    vec3 lookfrom;
    vec3 lookat;
    vec3 vup;

    vec3 w;
    vec3 u;
    vec3 v;
};*/

layout(std140, binding = 0) uniform cam
{
    vec4 theta_H_viewHeight_viewWidth;
    vec4 fov_aspRatio_aperture_focus_dist;
    vec4 focLength_lensRadius;

    vec4 origin;
    vec4 horizontal;
    vec4 vertical;
    vec4 lowerLeftCorner;

    vec4 lookFrom;
    vec4 lookAt;
    vec4 vup;

    vec4 w;
    vec4 u;
    vec4 v;
};// block_cam;

uniform bool shadows_active;
uniform bool renderImageOnly;
uniform int  renderDepth;
uniform float global_refract_index;

in vec3 fPos;
in vec3 camPos;

out vec4 FragColor;

const float infinity = 1. / 0.;
const float M_PI = 3.14159265359;

#define METAL 1
#define DIELECTRIC 2
#define LAMBERTIAN 3
#define PLAIN 4

#define SPHERE 0
#define TRIANGLE 1

#define EPSILON 0.0000001

vec4 globCol = vec4(0.0f);



uint wang_hash(inout uint seed)
{
    seed = uint(seed ^ uint(61)) ^ uint(seed >> uint(16));
    seed *= uint(9);
    seed = seed ^ (seed >> 4);
    seed *= uint(0x27d4eb2d);
    seed = seed ^ (seed >> 15);
    return seed;
}

float RandomFloat01(inout uint state)
{
    return float(wang_hash(state)) / 4294967296.0;
}

vec3 RandomUnitVector(inout uint state)
{
    float z = RandomFloat01(state) * 2.0f - 1.0f;
    float a = RandomFloat01(state);// * c_twopi;
    float r = sqrt(1.0f - z * z);
    float x = r * cos(a);
    float y = r * sin(a);
    return vec3(x, y, z);
}

vec3 random_pcg3d(uvec3 v) {
    v = v * 1664525u + 1013904223u;
    v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
    v ^= v >> 16u;
    v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
    return vec3(v) * (1.0 / float(0xffffffffu));
}

float rand()
{
    vec2 co = gl_FragCoord.xy;

    co *= 10.0; // Scale the coordinate system by 10
    vec2 ipos = floor(co);  // get the integer coords
    //vec2 fpos = fract(co);  // get the fractional coords

    return fract(sin(dot(ipos.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

float rand(float min, float max)
{
    vec2 co = gl_FragCoord.xy;

    co *= 10.0; // Scale the coordinate system by 10
    vec2 ipos = floor(co);  // get the integer coords

    return min + (max - min) * fract(sin(dot(ipos.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 random() {

    //return vec3(rand(), rand(), rand());
    uvec2 i = uvec2(gl_FragCoord.xy);
    return random_pcg3d(uvec3(i.x, i.y, 0));
}

vec3 randomVec3(float min, float max) {
    //return vec3(rand(min, max), rand(min, max), rand(min, max));
    uvec2 i = uvec2(gl_FragCoord.xy);
    return random_pcg3d(uvec3(i.x, i.y, 0));
}

vec3 random_in_unit_sphere() {

    //while (true) {
    for (int i = 0; i < 100;i++) {
        vec3 p = randomVec3(-1, 1);
        if ((length(p) * length(p)) >= 1) continue;
        return p;
    }
}

vec3 random_unit_vector() {
    return normalize(random_in_unit_sphere());
}

vec3 random_in_hemisphere(const vec3 normal) {
    vec3 in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

vec3 random_in_unit_disk() {
    for (int i = 0; i < 100; i++) {
        vec3 p = vec3(rand(-1, 1), rand(-1, 1), 0);
        if ((length(p)* length(p)) >= 1) continue;
        return p;
    }
}


/////////////---RAY---/////////////
struct Ray {
    vec3 origin, direction;
    vec4 color;
    int currDepth;
    bool inLense;
};

vec3 ray_at(vec3 origin, vec3 direction, double t) {
    return vec3(origin + t * direction);
}

bool ray_near_zero(vec3 ray_dir) {
    // Return true if the vector is close to zero in all dimensions.
    double s = 1e-8;
    return (abs(ray_dir.x) < s) && (abs(ray_dir.y) < s) && (abs(ray_dir.z) < s);
}

/////////////---SCENE-OBJECTS---/////////////
//SIZE 16 FLOAT
struct Sphere {

    vec4 positionAndradius;
    vec4 color;
    vec4 albedo;
    vec4 fuzzAndmat_ptr;
    vec4 spacer;
          
};

struct Triangle {

    vec4 position[3];
    vec4 normal;
    vec4 albedo;
    vec4 fuzzAndmat_ptr;

};

struct Vertex {

    vec4 position;
    vec4 normal;
    vec4 albedo;
    vec4 fuzzAndmat_ptr;

};

struct Node {

    vec4 data;
    vec4 indexLeftRight;
    vec4 min;
    vec4 max;
    vec4 PrimRangeLH;
};

struct drawObjectData {

    vec4 IDType;
    vec4 IndiLeftRightVertiLeftRight;
    mat4 translation;

};

/////////////---BUFFER-SCENE-DATA---/////////////
//layout(std430, binding = 0) buffer sphereSSBO
layout(std430, binding = 4) buffer sphereSSBO
{        
    Sphere sphereArr[];
    //material
};


//layout(std430, binding = 1) buffer vertSSBO
layout(std430, binding = 0) buffer vertSSBO
{
    Vertex vertArr[];
    //material
};

//layout(std430, binding = 2) buffer indiSSBO
layout(std430, binding = 1) buffer indiSSBO
{
    uint indiArr[];
    //material
};

//layout(std430, binding = 3) buffer nodeListSSBO
layout(std430, binding = 2) buffer nodeListSSBO
{
    float root;
    Node nodeList[];
    //material
};

//layout(std430, binding = 4) buffer PrimIndexSSBO
layout(std430, binding = 3) buffer PrimIndexSSBO
{
    float primList[];
    //material
};

layout(std430, binding = 5) buffer GameObjectSSBO
{
    drawObjectData drawObjectDataList[];
    //material
};

vec3 refract_s(in vec3 uv, in vec3 n, float etai_over_etat) {
    float cos_theta = min(dot(-uv, n), 1.0f);
    vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
    vec3 r_out_parallel = -sqrt(abs(1.0f - r_out_perp.length() * r_out_perp.length())) * n;
    return r_out_perp + r_out_parallel;

    //return refract(uv,-n,etai_over_etat);
}

struct Leaf { // a leaf node in the kd-tree
    // triangleCount > 0
    float triangleCount;
    int firstTriangle;
    int parentPointer;
};

struct Split { // an internal node in the kd-tree
    float splitValue; // data
    int leftChildPointer;
    int rightChildPointer;
    // sign bits store child type (leaf/split) and split axis
};

struct ExtendedSplit { // further data for kd-backtrack
    vec3 boundingBoxMin;
    vec3 boundingBoxMax;
    vec2 parentPointer;
};

struct TraversalState {
    vec2 nodePointer;
    float tmin;
    float tmax;
    // sign bits store node type (leaf, split)
    // and ray state (traverse, intersect, done)
};
struct IntersectState {
    vec2 triangleIndex;
    float triangleCount;
    float tmax;
};
struct HitState {
    vec2 bestTriangleIndex;
    float thit;
    float globalTmax;
};


/////////////---HIT-RECORD---/////////////

bool check_front_face(const Ray r, const vec3 outward_normal) {
    return dot(r.direction, outward_normal) < 0;
}

vec3 set_face_normal(const Ray r, const vec3 outward_normal) {
    return (check_front_face(r, outward_normal) ? outward_normal : -outward_normal);
}

struct hit_record {

    vec3 p;
    vec3 normal;
    double t;
    bool front_face;
    int mat_ptr;
    int arrIndex;
    int geoType;
    Triangle tri;

};


float hit_sphere(vec3 center, float radius, Ray r) {

    vec3 oc = r.origin - center;
    float a = length(r.direction)*length(r.direction);
    float half_b = dot(oc, r.direction);
    float c = length(oc)*length(oc) - radius * radius;
    float discriminant = half_b * half_b - a * c;

    if (discriminant < 0) {
        return -1.0f;
    }
    else {
        return (-half_b - sqrt(discriminant)) / (2.0 * a);
    }

}

double reflectance(double cosine, double ref_idx) {
    // Use Schlick's approximation for reflectance.
    double r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1 - r0) * pow(float(1 - cosine), 5);
 
}

bool sphere_hit( double t_min, in Ray r, double t_max, out hit_record rec,in Sphere sphere) {

    vec3 oc = r.origin - sphere.positionAndradius.xyz;
    double a = length(r.direction) * length(r.direction);
    double half_b = dot(oc, r.direction);
    double c = length(oc) * length(oc) - sphere.positionAndradius.w * sphere.positionAndradius.w;

    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return false;
    double sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    double root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = ray_at(r.origin, r.direction, rec.t);

    vec3 outward_normal = (rec.p - sphere.positionAndradius.xyz) / sphere.positionAndradius.w;
    rec.front_face = check_front_face(r,outward_normal);
    rec.normal = set_face_normal(r, outward_normal);
    rec.mat_ptr = int(sphere.fuzzAndmat_ptr.x);

    return true;
}

bool sphere_hit_simple(in Ray r, in Sphere sphere) {

    vec3 oc = r.origin - sphere.positionAndradius.xyz;
    double a = length(r.direction) * length(r.direction);
    double half_b = dot(oc, r.direction);
    double c = length(oc) * length(oc) - sphere.positionAndradius.w * sphere.positionAndradius.w;

    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return false;


    return true;
}

bool RayTriangleIntersect(
    Ray r,
    out hit_record rec,
    double t_min,
    double t_max,
    Triangle tri,
    out vec3 outIntersectionPoint)
{
    

    vec3 vertex0 = tri.position[0].xyz;
    vec3 vertex1 = tri.position[1].xyz;
    vec3 vertex2 = tri.position[2].xyz;

    vec4 normal = tri.normal;

    vec3 edge1, edge2, h, s, q;

    float a, f, u, v;

    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;

    h = cross(r.direction,edge2);
    a = dot(edge1,h);

    if (a > -EPSILON && a < EPSILON)
        return false;    // This ray is parallel to this triangle.

    f = 1.0 / a;
    s = r.origin - vertex0;
    u = f * dot(s,h);
    if (u < 0.0 || u > 1.0)
        return false;
    q = cross(s,edge1);
    v = f * dot(r.direction,q);
    if (v < 0.0 || u + v > 1.0)
        return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    float t = f * dot(edge2,q);

    if (t < t_min || t_max < t)
        return false;

    rec.t = t;
    rec.p = ray_at(r.origin, r.direction, rec.t);

    //vec3 outward_normal = tri.normal.xyz;//normalize(cross(edge1, edge2));//(rec.p - sphere.positionAndradius.xyz) / sphere.positionAndradius.w;
    vec3 outward_normal = normalize(cross(edge1, edge2));//(rec.p - sphere.positionAndradius.xyz) / sphere.positionAndradius.w;
    rec.front_face = check_front_face(r, outward_normal);
    rec.normal = set_face_normal(r, outward_normal);
    rec.mat_ptr = int(tri.fuzzAndmat_ptr.x);

    if (t > EPSILON) // ray intersection
    {
        outIntersectionPoint = r.origin + r.direction * t;
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
        return false;
}

vec3 triInterpolNormal(vec3 tri0,vec3 tri1, vec3 tri2,vec3 normal0, vec3 normal1, vec3 normal2, vec3 vertex) {   
    
    vec3 C;
    float u_n, v_n, w_n;

    vec3 v0v1 = tri1 - tri0;
    vec3 v0v2 = tri2 - tri0;

    vec3 N = cross(v0v1, v0v2);

    vec3 edge1 = tri2 - tri1;
    vec3 vp1 = vertex - tri1;

    C = cross(edge1, vp1);

    u_n = dot(N, C);

    vec3 edge2 = tri0 - tri2;
    vec3 vp2 = vertex - tri2;

    C = cross(edge2, vp2);

    v_n = dot(N, C);

    float denom = dot(N, N);

    u_n /= denom;
    v_n /= denom;
    w_n = 1.0f - u_n - v_n;


    return u_n * normal0 + v_n * normal1 + w_n * normal2;

}

bool metal_scatter(inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Sphere sphere){
    vec3 reflected = reflect(normalize(r_in.direction), rec.normal);
    scattered.origin = rec.p;
    scattered.direction = reflected;// +(sphere.fuzzAndmat_ptr.y * random_in_unit_sphere());
    attenuation = sphere.albedo;
    return (dot(scattered.direction, rec.normal) > 0);
}

bool metal_scatter(inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Triangle tri) {
    vec3 reflected = reflect(normalize(r_in.direction), rec.normal);
    scattered.origin = rec.p;
    scattered.direction = reflected;// +(sphere.fuzzAndmat_ptr.y * random_in_unit_sphere());
    attenuation = tri.albedo;
    return (dot(scattered.direction, rec.normal) > 0);
}


bool dielectric_scatter(
    inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Sphere sphere
) {
    attenuation = vec4(1.0, 1.0, 1.0, 1.0);
    double refraction_ratio = rec.front_face ? (1.0 / sphere.fuzzAndmat_ptr.x) : sphere.fuzzAndmat_ptr.x;


    vec3 unit_direction = normalize(r_in.direction);

    double cos_theta = min(dot(-unit_direction, rec.normal), 1.0f);
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    vec3 direction;
    //direction = refract(unit_direction, rec.normal, float(refraction_ratio));
    // checks wether the material can refrect or not
    if (refraction_ratio * sin_theta > 1.0 /* || reflectance(cos_theta, refraction_ratio) > random().y*/) {
        direction = reflect(unit_direction, rec.normal);    
    }
    else {
        //direction = refract_s(normalize(unit_direction), normalize(rec.normal), float(refraction_ratio));
        //direction = refract(unit_direction, rec.normal, 1.52f);
        direction = refract_s(normalize(unit_direction), normalize(rec.normal), global_refract_index);
    }

    scattered.origin = rec.p;
    scattered.direction = direction;
    
    return true;
}

bool dielectric_scatter(
     in Ray r_in, hit_record rec, inout vec4 attenuation, inout Ray scattered, Triangle tri
) {
    attenuation = vec4(1.0f, 1.0f, 1.0f, 1.0f);
    //double refraction_ratio = rec.front_face ? (1.0 / tri.fuzzAndmat_ptr.x) : tri.fuzzAndmat_ptr.x;
    double refraction_ratio = rec.front_face ? (1.0 / global_refract_index) : global_refract_index;

    vec3 unit_direction = normalize(r_in.direction);

    double cos_theta = min(dot(-unit_direction, rec.normal), 1.0f);
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    
    //direction = refract(normalize(unit_direction), normalize(rec.normal), 1.52f);
    vec3 direction = refract_s(unit_direction, rec.normal, global_refract_index);

    //direction = refract_s(normalize(r_in.direction), normalize(rec.normal), global_refract_index);
    //direction = refract_s(unit_direction, rec.normal, global_refract_index);
    // checks wether the material can refrect or not
    if (refraction_ratio * sin_theta > 1.0 /* || reflectance(cos_theta, refraction_ratio) > random().y*/) {
       direction = reflect(unit_direction, normalize(rec.normal));
    }
    else {
      direction = refract(unit_direction, normalize(rec.normal), float(refraction_ratio));
       // direction = refract(normalize(unit_direction), normalize(rec.normal), float(refraction_ratio));
    }

   // direction = refract(unit_direction, normalize(rec.normal), float(refraction_ratio));

    scattered.origin = rec.p;
    scattered.direction = direction;

    return true;
}



bool lambertian_scatter(
    inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Sphere sphere)
{

    vec3 scatter_direction = rec.normal +random_unit_vector();

    // Catch degenerate scatter direction
    //if (ray_near_zero(scattered.direction))
      //  scatter_direction = rec.normal;

    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = sphere.albedo;
    return true;
}

bool lambertian_scatter(
    inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Triangle tri)
{

    vec3 scatter_direction = rec.normal + random_unit_vector();

    // Catch degenerate scatter direction
    //if (ray_near_zero(scattered.direction))
      //  scatter_direction = rec.normal;

    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = tri.albedo;
    return true;
}

bool plain_color(
    inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Sphere sphere)
{

    vec3 scatter_direction = rec.normal;

    // Catch degenerate scatter direction
    //if (ray_near_zero(scattered.direction))
      //  scatter_direction = rec.normal;

    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = sphere.albedo;
    return true;
}

bool plain_color(
    inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Triangle tri)
{

    vec3 scatter_direction = rec.normal;

    // Catch degenerate scatter direction
    //if (ray_near_zero(scattered.direction))
      //  scatter_direction = rec.normal;

    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = tri.albedo;
    return true;
}

// box [0] = min, box [1] = max
bool ray_box_intersect(Ray r,vec4 min, vec4 max)
{
    float tmin = (min.x - r.origin.x) / r.direction.x;
    float tmax = (max.x - r.origin.x) / r.direction.x;

    if (tmin > tmax) {
        float tmp = tmin;
        tmin = tmax;
        tmax = tmp;
    }


    float tymin = (min.y - r.origin.y) / r.direction.y;
    float tymax = (max.y - r.origin.y) / r.direction.y;

    if (tymin > tymax) {
        float tmp = tymin;
        tymin = tymax;
        tymax = tmp;
    }

    if ((tmin > tymax) || (tymin > tmax))
        return false;

    if (tymin > tmin)
        tmin = tymin;

    if (tymax < tmax)
        tmax = tymax;

    float tzmin = (min.z - r.origin.z) / r.direction.z;
    float tzmax = (max.z - r.origin.z) / r.direction.z;

    if (tzmin > tzmax) {
        float tmp = tzmin;
        tzmin = tzmax;
        tzmax = tmp;
    }

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    if (tzmin > tmin)
        tmin = tzmin;

    if (tzmax < tmax)
        tmax = tzmax;

    return true;
}

/*
struct Node {

    vec4 data;
    vec4 indexLeftRight; + z = parent
    vec4 min;
    vec4 max;
    vec4 PrimRangeLH;
};
*/

bool kd_tree_search(double t_min, inout Ray ray, double t_max, inout hit_record rec) {
 

    hit_record temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;
    vec3 p;

    int currentNode = 217; int(root);
    bool upLeft = false,upRight = false;
    Triangle tri;
    int j = 0;
    
    while (j < 100){//primList.length()) {

        //TODO erster BUG hier, kommt nicht über fkt raus
       // globCol = vec4(currentNode + .5f, 0.0f, 0.0f, 1.0f);
        if ((!upLeft && !upRight) && ((nodeList[currentNode].indexLeftRight[0] != -1) || (nodeList[currentNode].indexLeftRight[1] != -1))) {
            
            // left child
            vec4 min = nodeList[int(nodeList[currentNode].indexLeftRight.x)].min,
                  max = nodeList[int(nodeList[currentNode].indexLeftRight.x)].max;
            if (ray_box_intersect(ray, min, max)) {
                globCol = vec4(currentNode + .5f, 0.0f, 0.0f, 1.0f);
                currentNode = int(nodeList[currentNode].indexLeftRight.x);
            }

            // right child
            min = nodeList[int(nodeList[currentNode].indexLeftRight.y)].min,
            max = nodeList[int(nodeList[currentNode].indexLeftRight.y)].max;
            if (ray_box_intersect(ray, min, max)) {
                globCol = vec4(currentNode + .5f, 0.0f, 0.0f, 1.0f);
                currentNode = int(nodeList[currentNode].indexLeftRight.y);
            }
            
            continue;
        }
        

        if ((!upLeft && !upRight) && ((nodeList[currentNode].indexLeftRight.x == -1) && (nodeList[currentNode].indexLeftRight.y == -1))) {
           // if ((int(nodeList[currentNode].PrimRangeLH.x) - int(nodeList[currentNode].PrimRangeLH.y)) <= 0)
             //   continue;
            
            for (int i = int(nodeList[currentNode].PrimRangeLH.x); i < int(nodeList[currentNode].PrimRangeLH.y); i++) {
               /* if (i >= primList.length() - 1 || i < 0)
                    return false;*/
                // TODO bug wohl hier

               /* tri.position[0] = vertArr[indiArr[i]].position;
                tri.position[1] = vertArr[indiArr[i]].position;
                tri.position[2] = vertArr[indiArr[i]].position;

                tri.albedo = vertArr[indiArr[i]].albedo;
                tri.fuzzAndmat_ptr = vertArr[indiArr[i]].fuzzAndmat_ptr;*/
                
                tri.position[0] = vec4(0.0f, 100.0f, 0.0f, 1.0f);
                tri.position[1] = vec4(0.0f, -100.0f, 0.0f, 1.0f);
                tri.position[2] = vec4(0.0f, 50.0f, 50.0f, 1.0f);



                tri.albedo = vec4(1.0f);
                tri.fuzzAndmat_ptr = vertArr[indiArr[0]].fuzzAndmat_ptr;

                globCol = vec4(0.0f, 1.0f, 0.0f, 1.0f);
                if (RayTriangleIntersect(ray, temp_rec, t_min, closest_so_far, tri, p)) {
                    hit_anything = true;
                    closest_so_far = temp_rec.t;
                    rec = temp_rec;
                    rec.arrIndex = i;
                    rec.geoType = TRIANGLE;
                    rec.tri = tri;
                }
               
                
                
            }
            if (!hit_anything) {
                
                if (!upLeft) {
                    upLeft = true;
                    continue;
                }
                if (!upRight) {
                    upRight = true;
                    continue;
                }
            }
            
            return hit_anything;
        }

        vec4 min = nodeList[int(nodeList[currentNode].indexLeftRight.y)].min,
            max = nodeList[int(nodeList[currentNode].indexLeftRight.y)].max;

        if (upLeft && ray_box_intersect(ray, min, max)) {
            
            // right node traversal
            currentNode = int(nodeList[currentNode].indexLeftRight.y);
            upLeft = false;
            
            continue;

        }
        
        if (upRight) {

            // go up to parent
            currentNode = int(nodeList[currentNode].indexLeftRight.z);
            upRight = false;
            
            continue;
        }
        
        j++;
    }
    
    return hit_anything;

}

/*struct drawObjectData {

    vec4 IDType;
    vec4 IndiLeftRightVertiLeftRight;

};*/

bool world_hit( double t_min, inout Ray r, double t_max, inout hit_record rec) {
    
    

        hit_record temp_rec;
        bool hit_anything = false;
        double closest_so_far = t_max;
        vec3 p;

        for (int i = 0; i < sphereArr.length(); i++) {
            if (sphere_hit(t_min, r, closest_so_far, temp_rec, sphereArr[i])) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
                rec.arrIndex = i;
                rec.geoType = SPHERE;

            }
        }

        // if (hit_anything)
          //  return hit_anything;



        Triangle tri;
        Sphere boundingSphere;
        vec3 boundingSpherePos;
        int offsetI = 0;

        for (int o = 0; o < drawObjectDataList.length(); o++) {
            //break;
            vec3 boundingSpherePos;
            if (o <= 0) {
                boundingSpherePos = vec3(lenseTransform * drawObjectDataList[o].translation * vec4(0.0f, 0.0f, 0.0f, 1.0f));
            }
            else
            {
                boundingSpherePos = vec3( drawObjectDataList[o].translation * vec4(0.0f, 0.0f, 0.0f, 1.0f));
            }
            boundingSphere.positionAndradius =   vec4(boundingSpherePos, drawObjectDataList[o].translation[0][0] * vertArr[0].fuzzAndmat_ptr.a);
            //for (int i = 0; i < indiArr.length(); i += 3) {
            for (int i = int(drawObjectDataList[o].IndiLeftRightVertiLeftRight.x);
                i < int(drawObjectDataList[o].IndiLeftRightVertiLeftRight.y); i += 3) {

                if (!sphere_hit_simple(r, boundingSphere)) break;
                
                

                if (o > 0) {
                    offsetI = int(drawObjectDataList[o].IDType.z);
                }

                
                //tri.position[0] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i]].position;
                //tri.position[1] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i + 1]].position;
                //tri.position[2] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i + 2]].position;
                if (o <= 0) {
                    tri.position[0] = lenseTransform * drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i]].position;
                    tri.position[1] = lenseTransform * drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i + 1]].position;
                    tri.position[2] = lenseTransform * drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i + 2]].position;

                    r.inLense = !r.inLense;
                }
                else {
                    tri.position[0] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i]].position;
                    tri.position[1] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i + 1]].position;
                    tri.position[2] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i + 2]].position;
                }
                //tri.normal = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i]].normal;
                

                tri.albedo = vertArr[indiArr[i]].albedo;
                //tri.albedo = tri.normal;
                tri.fuzzAndmat_ptr = vertArr[indiArr[i]].fuzzAndmat_ptr;
               // tri.fuzzAndmat_ptr.y = 4.0f;


                //  tri.fuzzAndmat_ptr = vertArr[indiArr[i]].fuzzAndmat_ptr;

                if (RayTriangleIntersect(r, temp_rec, t_min, closest_so_far, tri, p)) {
                    //mat4 inverseMat = (drawObjectDataList[o].translation);
                    /*
                    tri.normal =vec4(triInterpolNormal(tri.position[0].xyz, tri.position[1].xyz, tri.position[2].xyz,
                        ( vertArr[offsetI + indiArr[i]].normal).xyz, ( vertArr[offsetI + indiArr[i + 1]].normal).xyz, ( vertArr[offsetI + indiArr[i + 2]].normal).xyz, p.xyz), 1.0f);
                   */
                    
                    // tri.normal =vec4(triInterpolNormal(tri.position[0].xyz, tri.position[1].xyz, tri.position[2].xyz,  
                    //    (inverseMat * vertArr[offsetI + indiArr[i]].normal).xyz, (inverseMat * vertArr[offsetI + indiArr[i + 1]].normal).xyz, (inverseMat * vertArr[offsetI + indiArr[i + 2]].normal).xyz, p.xyz),1.0f);
                    
                    vec3 tmpNorm1 = (vertArr[offsetI + indiArr[i]].normal).xyz,
                         tmpNorm2 = (vertArr[offsetI + indiArr[i + 1]].normal).xyz,
                         tmpNorm3 = (vertArr[offsetI + indiArr[i + 2]].normal).xyz;



                   /* if (acos(abs(dot((vertArr[offsetI + indiArr[i]].normal).xyz, temp_rec.normal))) > M_PI * 0.75) {
                        tmpNorm1 = temp_rec.normal;
                    }

                    if (acos(abs(dot((vertArr[offsetI + indiArr[i + 1]].normal).xyz, temp_rec.normal))) > M_PI * 0.75) {
                        tmpNorm2 = temp_rec.normal;
                    }

                    if (acos(abs(dot((vertArr[offsetI + indiArr[i + 2]].normal).xyz, temp_rec.normal))) > M_PI * 0.75) {
                        tmpNorm3 = temp_rec.normal;
                    }*/

                    tri.normal = vec4(triInterpolNormal(tri.position[0].xyz, tri.position[1].xyz, tri.position[2].xyz,
                        tmpNorm1, tmpNorm2, tmpNorm3, p.xyz), 1.0f);
                    
                    //tri.normal = vec4(temp_rec.normal,1.0f);
                    hit_anything = true;
                    closest_so_far = temp_rec.t;
                    rec = temp_rec;
                    rec.arrIndex = i;
                    rec.geoType = TRIANGLE;
                    //tri.albedo = vec4(float(check_front_face(r, tri.normal.xyz)),0.0f,0.0f,1.0f);// * vec4(set_face_normal(r, tri.normal.xyz),1.0f);
                    rec.tri = tri;
                    //rec.normal = tri.normal.xyz;
                   
                    rec.front_face = check_front_face(r, tri.normal.xyz);
                    rec.normal = set_face_normal(r, tri.normal.xyz);
                    //rec.mat_ptr = int(tri.fuzzAndmat_ptr.x);


                }
            }
        }
    //return kd_tree_search(t_min, r, t_max, rec);
    return hit_anything;
}

vec4 ray_color2(inout Ray r, out bool hit, inout bool no_reflect, inout bool no_shadow) {
    
    hit_record rec;
    Ray scattered;
    vec4 attenuation;

    if (world_hit(0.01f, r, infinity, rec)) {

        if (rec.geoType == SPHERE) {

            switch (int(sphereArr[rec.arrIndex].fuzzAndmat_ptr.y)) {
            case LAMBERTIAN:

                if (lambertian_scatter(r, rec, attenuation, scattered, sphereArr[rec.arrIndex])) {

                    hit = true;

                    r = scattered;

                    return attenuation;


                }


                break;

            case DIELECTRIC:
                if (dielectric_scatter(r, rec, attenuation, scattered, sphereArr[rec.arrIndex])) {

                    hit = true;

                    r = scattered;

                    no_shadow = true;

                    return attenuation;


                }


                break;

            case METAL:
                if (metal_scatter(r, rec, attenuation, scattered, sphereArr[rec.arrIndex])) {

                    hit = true;

                    r = scattered;

                    return attenuation;

                }
                break;

            case PLAIN:
                if (plain_color(r, rec, attenuation, scattered, sphereArr[rec.arrIndex])) {

                    hit = true;
                    no_reflect = true;

                    r = scattered;

                    return attenuation;

                }
                break;
            default:
                vec4(0.0f, 0.0f, 0.0f, 1.0f);
                //hit = false;
                break;
            }
        }
        else {
            switch (int(rec.tri.fuzzAndmat_ptr.y)) {
            case LAMBERTIAN:

                if (lambertian_scatter(r, rec, attenuation, scattered, rec.tri)) {

                    hit = true;

                    r = scattered;

                    return attenuation;

                }


                break;

            case DIELECTRIC:
                if (dielectric_scatter(r, rec, attenuation, scattered,rec.tri)) {

                    hit = true;

                    r = scattered;

                    no_shadow = true;

                    return attenuation;


                }


                break;

            case METAL:
                if (metal_scatter(r, rec, attenuation, scattered, rec.tri)) {

                    hit = true;

                    r = scattered;
                    
                    return attenuation;

                }
                break;
            case PLAIN:
                if (plain_color(r, rec, attenuation, scattered, rec.tri)) {

                    hit = true;

                    r = scattered;

                    return attenuation;

                }
                break;
            default:
                vec4(0.0f, 0.0f, 0.0f, 1.0f);
                //hit = false;
                break;
            }
        }


    }
    
    no_shadow = true;
    hit = false;
    vec3 unit_direction = normalize(r.direction);
    double t = 0.5 * (unit_direction.y + 1.0);
   // return vec4(0.0f);//vec4((1.0 - t) * vec4(1.0, 1.0, 1.0, 1.0) + t * vec4(0.5, 0.7, 1.0, 1.0));
    return vec4((1.0 - t) * vec4(1.0, 1.0, 1.0, 1.0) + t * vec4(0.5, 0.7, 1.0, 1.0));
    
}

bool shadow_ray(vec3 origin, vec3 direction) {

    Ray ray;
    hit_record rec;

    ray.origin = origin;
    ray.direction = direction;

    if(world_hit(0.01f, ray, infinity, rec))
        return true;

    return false;
}

vec4 ray_color(in Ray r, int depth) {

    vec4 fColor  = vec4(1.0f);
    bool hit = false;
    vec4 color;

    
    for (int i = depth - 1; i >= 0; i--) {
        r.currDepth = depth;
        bool no_reflect = false;
        bool no_shadow =  false;
        color = ray_color2(r, hit, no_reflect, no_shadow);

        fColor *= float(!shadow_ray(r.origin, -normalize((r.origin - lightPos)))  || no_shadow || !shadows_active) * color;


        if (no_reflect) {
            break;
        }
    }

    if (hit) {
        return fColor;
    }

    vec3 unit_direction = normalize(r.direction);
    double t = 0.5 * (unit_direction.y + 1.0);
    return (vec4((1.0 - t) * vec4(1.0, 1.0, 1.0, 1.0) + t * vec4(0.5, 0.7, 1.0, 1.0)));


    //return vec4(1.0f, 0.0f, 0.0f, 1.0f);
}

float RandomVal(inout uint state, float min, float max) {

    state = state * 747796405 + 2891336453;
    uint result = ((state >> ((state >> 28) + 4)) ^ state) * 277803737;
    result = (result >> 22) ^ result;
    return (min + (max - min)) *  result / 4294967295.0f;
}

bool camera_get_ray(double s, double t, inout Ray ray, int o, bool lp) {
  
    if (lp) {
        // lense rays
        ray.origin = vec3((lowerLeftCorner.xyz + s * horizontal.xyz + t * vertical.xyz) + lookFrom.xyz).xyz;
        ray.direction = vec3((lenseTransform * drawObjectDataList[0].translation * vertArr[indiArr[o]].position) + 0.00005f) - ray.origin.xyz;
    }
    else {


        //Perspective Rays
        ray.origin = lookFrom.xyz;
        ray.direction = vec3(lowerLeftCorner.xyz + s * horizontal.xyz + t * vertical.xyz - vec3(lookFrom).xyz).xyz;
    }


    return true;
}

bool camera_get_debug_ray(double s, double t, inout Ray ray, int o, bool lp) {

    if (lp) {
        // lense rays

       
        //vec3 rayDir = initCamPos.xyz - initRayPos;

        ray.origin = vec3((initlowerLeftCorner.xyz + s * initHorizontal.xyz + t * initVertical.xyz) + initCamPos.xyz).xyz;
        ray.direction = vec3((lenseTransform * drawObjectDataList[0].translation * vertArr[indiArr[o]].position) + 0.00005f) - ray.origin.xyz;
    }
    else {


        //Perspective Rays
        ray.origin = initCamPos.xyz;// vec3((initlowerLeftCorner.xyz + 0.5 * initHorizontal.xyz + 0.5 * initVertical.xyz) - initCamPos.xyz).xyz;
        ray.direction = vec3((initlowerLeftCorner.xyz + s * initHorizontal.xyz + t * initVertical.xyz) - initCamPos.xyz).xyz;
    }


    return true;
}

bool camera_get_ray_retinal(double s, double t, out Ray ray) {
    // vec3 rd = vec3(camera.lens_radius * random_in_unit_disk());
    // vec3 offset = camera.u * rd.x + camera.v * rd.y;
    int offsetX = (1920 / 2);
    int offsetY = (1080 / 2);
    float phi = float(s/2),// +offsetX,
        lambda = asin(tanh(float(t/2)));





   // if (pow((gl_FragCoord.x* factXY.x) - offsetX ,2)+ pow((gl_FragCoord.y * factXY.y) - offsetY, 2) > 500000.0f)
   // if (pow((gl_FragCoord.x) - offsetX ,2)+ pow((gl_FragCoord.y) - offsetY, 2) > 500000.0f)
   //     return false;

    vec3 circleOrigin = vec3(viewPos);
    float r = factXY.x;

    //vec3 cp = normalize(vec3(lowerLeftCorner.xyz + s * horizontal.xyz + t * vertical.xyz - vec3(lookFrom).xyz) - circleOrigin);
    //vec3 vp = vec3(lowerLeftCorner.xyz + s * horizontal.xyz + t * vertical.xyz - vec3(lookFrom).xyz).xyz;
    float   x = r * sin(lambda) * cos(phi),
        y = r * sin(lambda) * sin(phi),
        z = r * cos(lambda);
    vec3 vp = vec3(x,y,z).xyz;
    ray.origin = vp;//normalize(circleOrigin - vp);

    //vec3 vp = normalize(vec3(lowerLeftCorner.xyz + s * horizontal.xyz + t * vertical.xyz - vec3(lookFrom).xyz).xyz - vec3(factXY.x, 0.0f, 0.0f));
    //vec3 vp = (viewPos- vec3(0.0f)) * factXY.x;
    ray.direction = circleOrigin - ray.origin ;//vec3(lowerLeftCorner.xyz + s * horizontal.xyz + t * vertical.xyz - vec3(lookFrom).xyz).xyz;// -offset);

    return true;
}

bool rayIntersectsPlane(vec3 rayOrigin, vec3 rayDirection, vec3 planeOrigin, vec3 planeNormal, vec3 planeU, vec3 planeV) {
    float d = dot(planeNormal, rayDirection);
    if (abs(d) < 1e-6) return false; // No intersection
    float t = dot(planeNormal, planeOrigin - rayOrigin) / d;
    if (t < 0.0) return false; // Intersection behind ray origin
    vec3 intersectionPoint = rayOrigin + rayDirection * t;
    vec3 localIntersectionPoint = intersectionPoint - planeOrigin;
    float u = dot(localIntersectionPoint, normalize(planeU));
    float v = dot(localIntersectionPoint, normalize(planeV));
    float planeWidth = length(planeU);
    float planeHeight = length(planeV);
    return u >= 0.0 && u <= planeWidth && v >= 0.0 && v <= planeHeight;
}



bool intersectRayPlane(vec3 rayOrigin, vec3 rayDir, vec3 planePoint, vec3 planeNormal, out float t) {
    float denom = dot(planeNormal, rayDir);

    if (abs(denom) > 1e-6) {
        vec3 planeToRayOrigin = planePoint - rayOrigin;
        t = dot(planeToRayOrigin, planeNormal) / denom;

        return (t >= 0.0);
    }

    return false;
}

bool intersectionRayCylinder(vec3 rayOrigin, vec3 rayDir, vec3 cylinderStart, vec3 cylinderEnd, float cylinderRadius, out float t) {
    vec3 AB = cylinderEnd - cylinderStart;
    vec3 AO = rayOrigin - cylinderStart;
    vec3 AOxAB = cross(AO, AB);
    vec3 VxAB = cross(rayDir, AB);
    float ab2 = dot(AB, AB);
    float a = dot(VxAB, VxAB);
    float b = 2.0 * dot(VxAB, AOxAB);
    float c = dot(AOxAB, AOxAB) - (cylinderRadius * cylinderRadius * ab2);
    float d = b * b - 4.0 * a * c;
    if (d < 0.0) {
        return false;
    }
    d = sqrt(d);
    t = (-b - d) / (2.0 * a);
    float t1 = (-b + d) / (2.0 * a);
    if (t < 0.0 && t1 < 0.0) {
        return false;
    }
    if (t < 0.0) {
        t = t1;
    }
    vec3 P = rayOrigin + rayDir * t;
    float y = dot(P - cylinderStart, AB) / ab2;
    if (y < 0.0 || y > 1.0) {
        return false;
    }

    // Test end caps
    float tCap1, tCap2;
    bool hitCap1 = intersectRayPlane(rayOrigin, rayDir, cylinderStart, normalize(AB), tCap1);
    bool hitCap2 = intersectRayPlane(rayOrigin, rayDir, cylinderEnd, normalize(AB), tCap2);

    if (hitCap1 && distance(rayOrigin + rayDir * tCap1, cylinderStart) <= cylinderRadius) {
        if (tCap1 > 0.0 && tCap1 < t) {
            t = tCap1;
            return true;
        }
    }

    if (hitCap2 && distance(rayOrigin + rayDir * tCap2, cylinderEnd) <= cylinderRadius) {
        if (tCap2 > 0.0 && tCap2 < t) {
            t = tCap2;
            return true;
        }
    }

    return true;
}





bool drawColoredRays(Ray r) {
    
    int raysPerPixel =  vertArr.length(); // minimum teilen durch 4 oder mehr
    bool shadow=false,reflect=false;
    float imageWidth = widthHeight.x, imageHeight = widthHeight.y;   
    int div = 500;
    //for (int i = 0; i < (imageWidth / div); i++) {
         //for (int j = 0; j < (imageHeight / div); j++) {
             for (int k = 0; k < raysPerPixel; k++) {
               // float sp = ( i * (imageWidth/ div)),
                 //   tp = ( j * (imageHeight/ div));
                Ray debugRay;
                camera_get_debug_ray(0.5f, 0.5f, debugRay, k * (vertArr.length() / raysPerPixel), true);
                //vec3 rayDir = vec3((initlowerLeftCorner.xyz + (0.45+0.05*i * (100 / imageWidth)) * initHorizontal.xyz + ((0.45 + 0.05 * j * (100 / imageHeight)) * initVertical.xyz) - initCamPos.xyz).xyz);
                vec3 rayDir = debugRay.direction;

                Ray camRay;
                camRay.origin = debugRay.origin;
                camRay.direction = rayDir;
                int depth = 5;

                struct BeginEnd {
                    vec3 p1, p2;
                };

                int hitCount = 0;
                BeginEnd rayArr[40];
                rayArr[0].p1 = camRay.origin;
                rayArr[0].p2 = camRay.origin + normalize(camRay.direction) * 100;
                for (int i = 1; i <= depth; i++) {
                    bool hitPoint = false;

                    ray_color2(camRay, hitPoint, shadow, reflect);

                    rayArr[i].p1 = camRay.origin;
                    rayArr[i].p2 = camRay.origin + normalize(camRay.direction) * 100;

                    if (hitPoint) {
                        hitCount++;
                        rayArr[i - 1].p2 = camRay.origin;
                    }
                    else {
                        rayArr[i - 1].p2 = camRay.origin + normalize(camRay.direction) * 100;
                        break;
                    }



                }

                float t;
                for (int i = 0; i <= hitCount; i++) {
                    if (intersectionRayCylinder(r.origin, r.direction, rayArr[i].p1, rayArr[i].p2, .05f, t)) {
                        return true;
                    }
                }

    }
    return false;
}



void main()
{   

    int depth = 3;
    //int depth = renderDepth;
    //+float imageWidth = 1920, imageHeight = 1080;
    float imageWidth = widthHeight.x, imageHeight = widthHeight.y;
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;   
    
    //float imageAspectRatio =  imageHeight / imageWidth;//imageWidth / imageHeight;
    //float Px = (2 * ((x + 0.5) / imageWidth) - 1) * tan(fov / 2 * M_PI / 180) * aspect_ratio;//imageAspectRatio;
    //float Py = (1 - 2 * ((y + 0.5) / imageHeight) * tan(fov / 2 * M_PI / 180));

    double s = (x /* + rand()*/) / ((imageWidth - 1));
    double t = (y /* + rand()*/) / ((imageHeight - 1));

    Ray r;

    int raysPerPixel = int(factXY.y);
    vec4 col = vec4(0.0f);

    if (renderImageOnly) {
        raysPerPixel = 1;// vertArr.length();

        uint rand = 100;

        int loops = 0;
        r.origin = vec3((lowerLeftCorner.xyz + s * horizontal.xyz + t * vertical.xyz) + lookFrom.xyz).xyz;
        for (int i = 0; i < raysPerPixel; i++) {
            camera_get_ray(s, t, r, i,false);
            col += ray_color(r, depth);

            //for (int i = 0; i < (imageWidth/300);i++) {
              //   for (int j = 0; j < (imageHeight/300); j++) {

          


                 //vec3 initRayPos = vec3((initlowerLeftCorner.xyz + 0.5 * initHorizontal.xyz + 0.5 * initVertical.xyz) - initCamPos.xyz).xyz;
                 //vec3 rayDir = initCamPos.xyz - initRayPos;
                 //vec3 rayDir = vec3((initlowerLeftCorner.xyz + (i * 300 / imageWidth) * initHorizontal.xyz + (i * 300 / imageHeight) * initVertical.xyz) - initCamPos.xyz).xyz;
                 //vec3 rayDir = vec3((initlowerLeftCorner.xyz + 0.5f * initHorizontal.xyz + 0.5f * initVertical.xyz) - initCamPos.xyz).xyz;
                     Ray debugRay;
                 //camera_get_debug_ray(s, t, debugRay, i, false);
                 if (rayIntersectsPlane(r.origin, r.direction, initlowerLeftCorner.xyz, cross(initHorizontal.xyz, initVertical.xyz), initHorizontal.xyz, initVertical.xyz))
                     col = vec4(1.0f, 0.0f, 0.0f, 1.0f);
				 
				 
                 if (drawColoredRays(r)) {
                     col = vec4(0.0f, 1.0f, 0.0f, 1.0f);
                 }



            loops++;
        }

        col /= loops;
    }
    else {
        raysPerPixel = 1;

        uint rand = 100;

        int loops = 0;
        //r.origin = vec3((lowerLeftCorner.xyz + s * horizontal.xyz + t * vertical.xyz) + lookFrom.xyz).xyz;

        for (int i = 0; i < raysPerPixel; i++) {

            camera_get_ray(s, t, r, i * (vertArr.length() / raysPerPixel), false);

            col += ray_color(r, depth);

            //vec3 initRayPos = vec3((initlowerLeftCorner.xyz + 0.5 * initHorizontal.xyz + 0.5 * initVertical.xyz) - initCamPos.xyz).xyz;
            //vec3 rayDir = initCamPos.xyz - initRayPos;
            //vec3 rayDir = vec3((initlowerLeftCorner.xyz + 0.5 * initHorizontal.xyz + 0.5 * initVertical.xyz) - initCamPos.xyz).xyz;

            if (rayIntersectsPlane(r.origin, r.direction, initlowerLeftCorner.xyz, cross(initHorizontal.xyz, initVertical.xyz), initHorizontal.xyz, initVertical.xyz))
                col = vec4(1.0f, 0.0f, 0.0f, 1.0f);
            
            //for (int i = 0; i < (imageWidth/100);i++) {
              //  for (int j = 0; j < (imageHeight/100); j++) {
                    //rayDir = vec3((initlowerLeftCorner.xyz + (i*100/ imageWidth) * initHorizontal.xyz + (j*100/ imageHeight) * initVertical.xyz) - initCamPos.xyz).xyz;
                    //rayDir = vec3((initlowerLeftCorner.xyz +initHorizontal.xyz + ( initVertical.xyz) - initCamPos.xyz).xyz);
                    //if (intersectCylinder(r, initCamPos.xyz, initCamPos.xyz + rayDir*100, .05f /*, out float t*/)) {

                    /*if (intersectRayCylinder(r, initCamPos.xyz, initCamPos.xyz + rayDir * 100, .05f )) {
                        col = vec4(0.0f, 1.0f, 0.0f, 1.0f);
                    }*/
                    /*Ray colRay;
                    colRay.origin = initCamPos.xyz;
                    colRay.direction = rayDir;*/

             //   }
            //}
            if (drawColoredRays(r)) {
                col = vec4(0.0f, 1.0f, 0.0f, 1.0f);
            }
            loops++;
        }

       col /= loops;

    }

    FragColor = col;
};