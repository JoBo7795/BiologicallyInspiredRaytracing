#version 460


uniform vec3 view;
uniform vec3 viewPos;
uniform vec3 up;
uniform vec3 lensPos;
uniform vec3 lightPos;
uniform vec2 factXY;
uniform mat4 lenseTransform;
uniform vec2 widthHeight;
uniform vec3 initCamPos;

uniform vec3 initlowerLeftCorner;
uniform vec3 initHorizontal;
uniform vec3 initVertical;


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
};

uniform bool shadows_active;
uniform bool debugMode;
uniform bool renderImageOnly;
uniform int  renderDepth;
uniform float global_refract_index;
uniform float sceneBrightness;

in vec3 fPos;
in vec3 camPos;

out vec4 FragColor;

const float infinity = 1. / 0.;
const float M_PI = 3.14159265359;


vec3 circlePoints[1000];
int numPoints = 0, notLitPoints = 0;

#define METAL 1
#define DIELECTRIC 2
#define LAMBERTIAN 3
#define PLAIN 4

#define SPHERE 0
#define TRIANGLE 1
#define LENSE 2
#define  IRIS 3

#define EPSILON 0.0000001

vec4 globCol = vec4(0.0f);


/////////////---RANDOM-FUNCTIONS---/////////////

 vec3 random_pcg3d(uvec3 v) {
     v = v * 1664525u + 1013904223u;
     v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
     v ^= v >> 16u;
     v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
     return vec3(v) * (1.0 / float(0xffffffffu));
 }

 vec3 randomVec3(float min, float max) {
    
     uvec2 i = uvec2(gl_FragCoord.xy);
     return random_pcg3d(uvec3(i.x, i.y, 0));
 }

 vec3 random_in_unit_sphere() {
 
     for (int i = 0; i < 100;i++) {
         vec3 p = randomVec3(-1, 1);
         if ((length(p) * length(p)) >= 1) continue;
         return p;
     }
 }

vec3 random_unit_vector() {
    return normalize(random_in_unit_sphere());
}



/////////////---RAY---/////////////
struct Ray {
    vec3 origin, direction;
    vec4 color;
    int currDepth;
    bool inLense;
    vec3 hit_normal_dir;
    vec3 hit_normal_orig;
    int rgb;
    float brightness;
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

struct drawObjectData {

    vec4 IDType;
    vec4 IndiLeftRightVertiLeftRight;
    mat4 translation;

};

struct LensData {

    vec3 lensOrigin;

    vec3 lensOriginLeft;
    vec3 lensOriginRight;

    double lensRadiusLeft, lensRadiusRight;
    double d;
    double n;

};

// lens parameters
double dpt, b, g, f;
LensData lens;
LensData lensArr[2];
int numLens = 1;

// construct lens
void lens_data_construct(inout LensData lens, vec3 lensOrigin_in = vec3(0.0), float n_in = 1.59, vec3 lensOriginLeft_in = vec3(0.9 * 2, 0.0, 0.0), vec3 lensOriginRight_in = vec3(0.9 * 2, 0.0, 0.0), double lensRadiusLeft_in = 2.0, double lensRadiusRight_in = 2.0) {
    lens.lensOrigin = lensOrigin_in;

    lens.lensOriginLeft = lensOrigin_in - lensOriginLeft_in;
    lens.lensOriginRight = lensOrigin_in + lensOriginRight_in;

    lens.lensRadiusLeft = lensRadiusLeft_in;
    lens.lensRadiusRight = lensRadiusRight_in;

    lens.n = n_in;

    lens.d = abs(lens.lensRadiusLeft) + abs(lens.lensRadiusRight) - length(lens.lensOriginLeft - lens.lensOriginRight);
}

// calculate parameter 1/f for lens
double calc_one_by_f(inout LensData lens) {
    return (lens.n - 1.0) * ((1.0 / lens.lensRadiusLeft) - (1.0 / (-lens.lensRadiusRight)) + ((lens.n - 1.0) * lens.d) / (lens.n * lens.lensRadiusLeft * (-lens.lensRadiusRight)));
}

// calculate parameter f for lens
double calc_f(inout LensData lens) {
    return 1.0 / ((lens.n - 1.0) * ((1.0 / lens.lensRadiusLeft) - (1.0 / (-lens.lensRadiusRight)) + ((lens.n - 1.0) * lens.d) / (lens.n * lens.lensRadiusLeft * (-lens.lensRadiusRight))));
}

// calculate picture width b
void calc_image_dist(inout LensData lens,float offsetVal = 0.0f) {

    dpt = calc_one_by_f(lens) + offsetVal + 1.0f;
    b = 1.0 / (dpt - (1.0 / g));

}



/////////////---BUFFER-SCENE-DATA---/////////////

layout(std430, binding = 4) buffer sphereSSBO
{        
    Sphere sphereArr[];

};



layout(std430, binding = 0) buffer vertSSBO
{
    Vertex vertArr[];

};


layout(std430, binding = 1) buffer indiSSBO
{
    uint indiArr[];

};


layout(std430, binding = 3) buffer PrimIndexSSBO
{
    float primList[];

};

layout(std430, binding = 5) buffer GameObjectSSBO
{
    drawObjectData drawObjectDataList[];

};



/////////////---GEOMETRY/RAY-COLLISION-FUNCTIONS---/////////////

// checks if a ray hits the front face of a sphere
bool check_front_face(const Ray r, const vec3 outward_normal) {
    return dot(r.direction, outward_normal) < 0;
}

// set face normal of a sphere, if not front face invert it
vec3 set_face_normal(const Ray r, const vec3 outward_normal) {
    return (check_front_face(r, outward_normal) ? outward_normal : -outward_normal);
}

// struct to record hit information data
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

// Disc struct for iris
struct Disc {
    vec3 origin, normal;
    float radius;
};

// function to check for sphere ray intersection, stores information in hit_record struct
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


// reduced sphere hit check, doesnt store information in hit_record struct
bool sphere_hit_simple(in Ray r, in Sphere sphere) {

    vec3 oc = r.origin - sphere.positionAndradius.xyz;
    double a = length(r.direction) * length(r.direction);
    double half_b = dot(oc, r.direction);
    double c = length(oc) * length(oc) - sphere.positionAndradius.w * sphere.positionAndradius.w;

    double discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return false;


    return true;
}

// ray sphere intersection function, special for lens calculation
bool intersect_ray_sphere(vec3 rayOrigin, vec3 rayDirection, vec3 sphereOrigin, double sphereRadius, inout vec3 iPoint, inout double t) {
    vec3 oc = vec3(rayOrigin - sphereOrigin);
    double a = dot(rayDirection, rayDirection);
    double b = 2.0 * dot(oc, rayDirection);
    double c = dot(oc, oc) - sphereRadius * sphereRadius;
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return false;
    }
    else {
        double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
        double t2 = (-b + sqrt(discriminant)) / (2.0 * a);
        if (t1 < 0 && t2 < 0) {
            // both intersection points in negative direction of the ray
            return false;
        }
        else if (t1 > 0 && t2 > 0) {
            // both intersection points in positive direction of the ray
            // take smaller value of t
            iPoint = rayOrigin + float(min(t1, t2)) * rayDirection;
            t = min(t1, t2);
        }
        else {
            // one intersection points in positive direction of the ray
            // take bigger value of t
            iPoint = rayOrigin + float(max(t1, t2)) * rayDirection;
            t = max(t1, t2);
        }
        return true;
    }
}

// function to calculate a intersection with a ray and a lens
bool intersect_ray_lens(vec3 rayOrigin, vec3 rayDirection, LensData lens, out vec3 normal, inout float t, bool inLense) {

    vec3 intersection1;
    vec3 intersection2;
    double t1 = 0.0, t2 = 0.0;

    bool inter1 = intersect_ray_sphere(rayOrigin, normalize(rayDirection), lens.lensOriginLeft, lens.lensRadiusLeft, intersection1, t1);
    bool inter2 = intersect_ray_sphere(rayOrigin, normalize(rayDirection), lens.lensOriginRight, lens.lensRadiusRight, intersection2, t2);

    float angle1 = degrees(acos(dot(normalize((lens.lensOrigin - lens.lensOriginRight)), normalize(lens.lensOrigin - intersection1))));
    float angle2 = degrees(acos(dot(normalize((lens.lensOrigin - lens.lensOriginLeft)), normalize(lens.lensOrigin - intersection2))));

    if (inter2 && (angle1 <= 90) && inLense) {

        t = float(t1);
        normal = normalize(vec3(intersection1 - lens.lensOriginLeft));

        return true;
    }

    if (inter2 && (angle2 <= 90) && !inLense) {

        t = float(t2);
        normal = normalize(vec3(intersection2 - lens.lensOriginRight));

        return true;
    }

    return false;
}

// Function to check if ray intersects with triangle, stores hit data in hit_record struct
bool ray_triangle_intersect(
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


    vec3 outward_normal = normalize(cross(edge1, edge2));
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


// function to check for end caps of ray cylinder
bool check_end_caps(vec3 rayOrigin, vec3 rayDir, vec3 planePoint, vec3 planeNormal, out float t) {
    float denom = dot(planeNormal, rayDir);

    if (abs(denom) > 1e-6) {
        vec3 planeToRayOrigin = planePoint - rayOrigin;
        t = dot(planeToRayOrigin, planeNormal) / denom;

        return (t >= 0.0);
    }

    return false;
}

bool intersection_ray_cylinder(vec3 rayOrigin, vec3 rayDir, vec3 cylinderStart, vec3 cylinderEnd, float cylinderRadius, out float t) {
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
    bool hitCap1 = check_end_caps(rayOrigin, rayDir, cylinderStart, normalize(AB), tCap1);
    bool hitCap2 = check_end_caps(rayOrigin, rayDir, cylinderEnd, normalize(AB), tCap2);

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

// function to interpolate surface normal by vertice normal and position data
// tri 0 - 2 vertice positions, normal 0 - 2 vertice normals, vertex is position to interpolate
vec3 triInterpolNormal(vec3 tri0, vec3 tri1, vec3 tri2, vec3 normal0, vec3 normal1, vec3 normal2, vec3 vertex) {
    
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



// Ray disc intersection for iris calculation
bool does_ray_intersect_disc(inout Ray ray, Disc disc) {

    float dotProduct = ray.direction.x * disc.normal.x + ray.direction.y * disc.normal.y + ray.direction.z * disc.normal.z;

    // if dot product is 0 ray and disc are parallel and don't intersect
    if (dotProduct == 0) {
        return false;
    }

    // vector from ray origin to disc origin
    vec3 originToDisc = { disc.origin.x - ray.origin.x, disc.origin.y - ray.origin.y, disc.origin.z - ray.origin.z };

    // distance of ray and disc origin
    float t = (originToDisc.x * disc.normal.x + originToDisc.y * disc.normal.y + originToDisc.z * disc.normal.z) / dotProduct;

    // if t negative, disc is behind camera
    if (t < 0) {
        return false;
    }

    // calculate intersection point
    vec3 intersection = { ray.origin.x + t * ray.direction.x, ray.origin.y + t * ray.direction.y, ray.origin.z + t * ray.direction.z };


    float distanceSquared = (intersection.x - disc.origin.x) * (intersection.x - disc.origin.x) +
        (intersection.y - disc.origin.y) * (intersection.y - disc.origin.y) +
        (intersection.z - disc.origin.z) * (intersection.z - disc.origin.z);

    // if distance squared bigger than quadratic radius of disc, intersection located outside of disc
    if (distanceSquared > (disc.radius * disc.radius)) {
        return true;
    }

    return false;
}



/////////////---SCATTER-AND-COLORING-FUNCTIONS---/////////////

// function to calculate metallic light scatter for sphere
bool metal_scatter(inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Sphere sphere){
    vec3 reflected = reflect(normalize(r_in.direction), rec.normal);
    scattered.origin = rec.p;
    scattered.direction = reflected;
    attenuation = sphere.albedo;
    return (dot(scattered.direction, rec.normal) > 0);
}

// function to calculate metallic light scatter for triangle
bool metal_scatter(inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Triangle tri) {
    vec3 reflected = reflect(normalize(r_in.direction), rec.normal);
    scattered.origin = rec.p;
    scattered.direction = reflected;
    attenuation = tri.albedo;
    return (dot(scattered.direction, rec.normal) > 0);
}


// function to calculate dielectric light scatter for sphere
bool dielectric_scatter(
    inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Sphere sphere
) {
    attenuation = vec4(1.0, 1.0, 1.0, 1.0);

    double refVal;

    if (r_in.rgb == 0) {
        refVal = 1.59;

    }
    else if (r_in.rgb == 1) {

        refVal = 1.59;
    }
    else {
        refVal = 1.59;

    }

    double refraction_ratio = rec.front_face ? (1.0f / refVal) : (refVal);


    vec3 unit_direction = normalize(r_in.direction);

    double cos_theta = min(dot(-unit_direction, rec.normal), 1.0f);
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);   

    vec3 direction = refract(unit_direction, normalize(rec.normal), float(refraction_ratio));

    scattered = r_in;    

    scattered.origin = rec.p;
    scattered.direction = normalize(direction);
    
    return true;
}

// function to calculate dielectric light scatter for triangle
bool dielectric_scatter(
    in Ray r_in, hit_record rec, inout vec4 attenuation, inout Ray scattered, Triangle tri
) {
    attenuation = vec4(1.0f, 1.0f, 1.0f, 1.0f);
    double refraction_ratio = rec.front_face ? (1.0f / global_refract_index) : ( global_refract_index);
    vec3 normal = rec.front_face ? normalize(-rec.normal) : normalize(rec.normal);
    vec3 unit_direction = normalize(r_in.direction);

    double cos_theta = min(dot(-unit_direction, rec.normal), 1.0f);
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    vec3 direction;

    direction = refract(unit_direction, normalize(rec.normal), float(refraction_ratio));

    scattered = r_in;

    scattered.origin = rec.p;
    scattered.direction = direction;

    return true;
}

// function to calculate lambertian light scatter for sphere
bool lambertian_scatter(
    inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Sphere sphere)
{

    vec3 scatter_direction = rec.normal +random_unit_vector();

    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = sphere.albedo;
    return true;
}

// function to calculate lambertian light scatter for triangle
bool lambertian_scatter(
    inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Triangle tri)
{

    vec3 scatter_direction = rec.normal + random_unit_vector();

    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = tri.albedo;
    return true;
}

// function to calculate plain non-reflecting light scatter for sphere
bool plain_color(
    inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Sphere sphere)
{

    vec3 scatter_direction = rec.normal;


    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = sphere.albedo;
    return true;
}

// function to calculate plain non-reflecting light scatter for triangle
bool plain_color(
    inout Ray r_in, inout hit_record rec, inout vec4 attenuation, inout Ray scattered, in Triangle tri)
{

    vec3 scatter_direction = rec.normal;

    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = tri.albedo;
    return true;
}



// function to check if different scene geometry gets hit by ray
bool world_hit( double t_min, inout Ray r, double t_max, inout hit_record rec) {
    
    
     Disc iris;
     iris.origin = vec3(0.0);
     iris.normal = vec3(1.0, 0.0, 0.0);
     iris.radius = .01;

     hit_record temp_rec;
     bool hit_anything = false;
     double closest_so_far = t_max;
     vec3 p;

     // check for spheres
     for (int i = 0; i < sphereArr.length(); i++) {
         if (sphere_hit(t_min, r, closest_so_far, temp_rec, sphereArr[i])) {
             hit_anything = true;
             closest_so_far = temp_rec.t;
             rec = temp_rec;
             rec.arrIndex = i;
             rec.geoType = SPHERE;

         }
     }


     float t;
     vec3 normal_out;
     // check for ray intersection of ray with lens(es)
     for (int i = 0; i < numLens; i++) {            

         if (intersect_ray_lens(r.origin, r.direction, lensArr[i], normal_out, t, r.inLense)) {
             if (closest_so_far > t) {
                 hit_anything = true;
                 rec = temp_rec;

                 rec.t = t;
                 closest_so_far = t;

                 rec.p = ray_at(r.origin, r.direction, t);

                 vec3 outward_normal = normal_out;
                 rec.front_face = check_front_face(r, outward_normal);
                 rec.normal = set_face_normal(r, outward_normal);
                 rec.mat_ptr = DIELECTRIC;

                 rec.arrIndex = 3;
                 rec.geoType = SPHERE;
                 r.inLense = !r.inLense;


                 r.hit_normal_dir = normalize(normal_out);
                 r.hit_normal_orig = rec.p;

             }
         }
     }

     Triangle tri;
     Sphere boundingSphere;
     vec3 boundingSpherePos;
     int offsetI = 0;


     // check for intersection of ray with triangle objects
     for (int o = 0; o < drawObjectDataList.length(); o++) {

         // use bounding sphere to avoid unneccessary calculations
         vec3 boundingSpherePos;
         if (o <= 0) {
             boundingSpherePos = vec3(lenseTransform * drawObjectDataList[o].translation * vec4(0.0f, 0.0f, 0.0f, 1.0f));
         }
         else
         {
             boundingSpherePos = vec3( drawObjectDataList[o].translation * vec4(0.0f, 0.0f, 0.0f, 1.0f));
         }
         boundingSphere.positionAndradius =   vec4(boundingSpherePos, drawObjectDataList[o].translation[0][0] * vertArr[0].fuzzAndmat_ptr.a);

         for (int i = int(drawObjectDataList[o].IndiLeftRightVertiLeftRight.x);
             i < int(drawObjectDataList[o].IndiLeftRightVertiLeftRight.y); i += 3) {              
             

             if (o > 0) {
                 offsetI = int(drawObjectDataList[o].IDType.z);
             }

             tri.position[0] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i]].position;
             tri.position[1] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i + 1]].position;
             tri.position[2] = drawObjectDataList[o].translation * vertArr[offsetI + indiArr[i + 2]].position;
             
             // albedo and material are stored in all vertices identically
             tri.albedo = vertArr[indiArr[i]].albedo;
             tri.fuzzAndmat_ptr = vertArr[indiArr[i]].fuzzAndmat_ptr;


             if (ray_triangle_intersect(r, temp_rec, t_min, closest_so_far, tri, p)) {

                 vec3 tmpNorm1 = (vertArr[offsetI + indiArr[i]].normal).xyz,
                      tmpNorm2 = (vertArr[offsetI + indiArr[i + 1]].normal).xyz,
                      tmpNorm3 = (vertArr[offsetI + indiArr[i + 2]].normal).xyz;

                 vec3 tmpCol0 = (vertArr[offsetI + indiArr[i]].albedo).xyz,
                      tmpCol1 = (vertArr[offsetI + indiArr[i + 1]].albedo).xyz,
                      tmpCol2 = (vertArr[offsetI + indiArr[i + 2]].albedo).xyz;

                 // interpolate normal of triangle
                 tri.normal = vec4(triInterpolNormal(
                     tri.position[0].xyz, tri.position[1].xyz, tri.position[2].xyz, tmpNorm1, tmpNorm2, tmpNorm3, p.xyz), 1.0f);
                 
                 // interpolate color of triangle
                 tri.albedo = vec4(triInterpolNormal(
                     tri.position[0].xyz, tri.position[1].xyz, tri.position[2].xyz, tmpCol0.xyz, tmpCol1.xyz, tmpCol2.xyz, p.xyz),1.0f);

                 // fill hit struct
                 hit_anything = true;
                 closest_so_far = temp_rec.t;
                 rec = temp_rec;
                 rec.arrIndex = i;
                 rec.geoType = TRIANGLE;
                 rec.tri = tri;
                 rec.front_face = check_front_face(r, tri.normal.xyz);
                 rec.normal = set_face_normal(r, tri.normal.xyz);


                 r.hit_normal_dir = tri.normal.xyz;
                 r.hit_normal_orig = p;

             }
         }
     }

    return hit_anything;
}

// check if ray hits objects in the world and determine resulting color
vec4 world_hit_check(inout Ray r, inout bool hit, inout bool no_reflect, inout bool no_shadow) {
    
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
                    r.brightness = 1.0f;// sceneBrightness;

                    return attenuation;

                }
                break;
            default:
                vec4(0.0f, 0.0f, 0.0f, 1.0f);
                break;
            }
        }
        else if (rec.geoType == LENSE) {

            // Define dummy sphere
            Sphere sphere;
            sphere.fuzzAndmat_ptr.x = 1.0f;
            sphere.fuzzAndmat_ptr.y = DIELECTRIC;
            sphere.color = vec4(1.0f, 0.0f, 0.0f, 1.0f);
            sphere.albedo = vec4(0.0f,0.0f,1.0f,1.0f);

            if (dielectric_scatter(r, rec, attenuation, scattered, sphere)) {

                hit = true;

                r =  scattered;

                no_shadow = true;
                
                return attenuation;
            }          


        }
        else if (rec.geoType == IRIS) {

            hit = true;

            notLitPoints++;

            return vec4(0.0);
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
                    no_reflect = true;

                    r = scattered;

                    r.brightness = 1.0f;// sceneBrightness;
                    return attenuation;

                    // apply phong model on plain color
                    vec3 materialAmbient = vec3(0.05, 0.05, 0.05);
                    vec3 materialDiffuse = vec3(0.5, 0.5, 0.5);
                    vec3 materialSpecular = vec3(0.01, 0.01, 0.01);
                    
                    
                    vec3 fpos = r.origin;
                    
                    vec3 lightAmbient = vec3(0.2f, 0.2f, 0.2f);
                    vec3 lightDiffuse = vec3(0.5f,  0.5f, 0.5f);
                    vec3 lightSpecular = vec3(-10.0f, 100.0f, -10.0f);
                    
                    float materialShininess = 10;
                    
                    // Diffuse 
                    vec3 norm = normalize(rec.tri.normal.xyz);
                    vec3 lightDir = normalize(lightPos - fpos);
                    float diff = max(dot(norm, lightDir), 0.0);
                    vec3 diffuse = lightDiffuse * (diff * materialDiffuse);
                    
                    // Specular
                    vec3 viewDir = normalize(-fpos);
                    vec3 reflectDir = reflect(-lightDir, norm);
                    float spec = pow(max(dot(viewDir, reflectDir), 0.0), materialShininess);
                    vec3 specular = lightSpecular * (spec * materialShininess);

                    r.brightness = 1.0f;// sceneBrightness;
                                       

                    return attenuation +vec4((diffuse + specular), 1.0);

                }
                break;
            default:
                r.brightness = 1.0f;
                vec4(0.0f, 0.0f, 0.0f, 1.0f);

                break;
            }
        }


    }
    
    // if no hit detected background color returned
    no_shadow = true;
    hit = false;
    vec3 unit_direction = normalize(r.direction);
    double t = 0.5 * (unit_direction.y + 1.0);
    r.brightness = 1.0;
    return vec4((1.0 - t) * vec4(1.0, 1.0, 1.0, 1.0) + t * vec4(0.5, 0.7, 1.0, 1.0));
    
}

// check for shadow
bool shadow_ray(vec3 origin, vec3 direction) {

    Ray ray;
    hit_record rec;

    ray.origin = origin;
    ray.direction = direction;

    if(world_hit(0.01f, ray, infinity, rec))
        return true;

    return false;
}

// Determine color of ray
vec4 ray_color(inout Ray r, int depth) {

    vec4 fColor  = vec4(1.0f);
    bool hit = false;
    vec4 color;

    
    for (int i = depth - 1; i >= 0; i--) {
        r.currDepth = depth;
        bool no_reflect = false;
        bool no_shadow =  false;
        color = world_hit_check(r, hit, no_reflect, no_shadow);

        fColor *= float(!shadow_ray(r.origin, -normalize((r.origin - lightPos)))  || no_shadow || !shadows_active) * color;

        if (no_reflect) {
            break;
        }
    }

    if (hit) {
        return fColor;
    }

    // if no hit, return background color
    vec3 unit_direction = normalize(r.direction);
    double t = 0.5 * (unit_direction.y + 1.0);
    r.brightness = 1.0;
    return (vec4((1.0 - t) * vec4(1.0, 1.0, 1.0, 1.0) + t * vec4(0.5, 0.7, 1.0, 1.0)));

}

// helper function for addCirclePoints
vec3 rotate_point(vec3 point, vec3 center, vec3 axis, float angle) {
    vec3 v = point - center;
    vec3 k = axis / length(axis);
    return center + v * cos(angle) + cross(k, v) * sin(angle) + k * dot(k, v) * (1.0 - cos(angle));
}

// function to add circular sample points
// numcircles is the numbers of cirlces where numPointsPerCircle are located on
// numPointsPerCircle is the number of points per circle
void add_circle_points(vec3 center, float radius, int numCircles, int numPointsPerCircle, vec3 normal) {
    float radiusStep = radius / float(numCircles);
    for (int j = 0; j < numCircles; j++) {
        float currentRadius = radiusStep * float(j + 1);
        for (int i = 0; i < numPointsPerCircle; i++) {
            float angle = float(i) / float(numPointsPerCircle) * 6.28;
            vec3 point = center + vec3(cos(angle), sin(angle), 0.0) * currentRadius;
            point = rotate_point(point, center, normal, angle);
            circlePoints[numPoints] = point;
            numPoints++;
        }
    }
}

// create ray based on position s and t 
bool camera_get_ray(double s, double t, inout Ray ray, int o, bool lp) {
  
    if (lp) {
        // lense rays

        ray.origin = vec3((lowerLeftCorner.xyz + s * horizontal.xyz + t * vertical.xyz)).xyz;
        ray.origin.x -= 1.0f;
        ray.origin.x -= 0.5f;
        ray.direction = normalize(vec3((circlePoints[o])) - ray.origin.xyz);
    }
    else {


        //Perspective Rays
        ray.origin = lookFrom.xyz;
        ray.direction = normalize(vec3(lowerLeftCorner.xyz + s * horizontal.xyz + t * vertical.xyz - vec3(lookFrom).xyz).xyz);

    }


    return true;
}

// create debug view ray, for external debug view
bool camera_get_debug_ray(double s, double t, inout Ray ray, int o, bool lp) {

    if (lp) {
        // lense rays
        vec3 origin;
        ray.origin = vec3((initlowerLeftCorner + s * initHorizontal + t * initVertical));
       
        ray.origin.x -= 1; 

        b = 2.136;

        if (o == 0)
            ray.direction = normalize(vec3(vec3(0.0f, 0.15f, 0.0f)) - ray.origin);
        if (o == 1)
            ray.direction = normalize(vec3(vec3(0.0f, 0.0f, 0.0f)) - ray.origin);
        if (o == 2)
            ray.direction = normalize(vec3(vec3(0.0f, -0.15f, 0.0f)) - ray.origin);
        if (o == 3)
            ray.direction = normalize(vec3(vec3(0.0f, 0.15f, 0.0f)) - ray.origin);
        if (o == 4)
            ray.direction = normalize(vec3(vec3(0.0f, 0.0f, 0.0f)) - ray.origin);
        if (o == 5)
            ray.direction = normalize(vec3(vec3(0.0f, -0.15f, 0.0f)) - ray.origin);



    }
    else {


        //Perspective Rays
        ray.origin = initCamPos.xyz;
        ray.direction = vec3((initlowerLeftCorner.xyz + s * initHorizontal.xyz + t * initVertical.xyz) - initCamPos.xyz).xyz;
    }


    return true;
}

// function to draw bouncing rays through scene
bool draw_colored_rays(inout Ray r) {
    
    int raysPerPixel = 6; // number of rays originating in one pixel
    bool shadow=false;
    float imageWidth = widthHeight.x, imageHeight = widthHeight.y;   
    int div = 5;
    float valX = 0.0f;
    float valY = 0.0f;
    float partX = (imageWidth / raysPerPixel) / imageWidth;
    float partY = (imageHeight / raysPerPixel) / imageHeight;
    bool inLense = false;

    // traverse rays through scene
    int step = numPoints / 3;
    for (int col = 0; col < 1; col++) {
        valX = 0;
        for (int l = 0; l < div; l++) {

            valY = 0.0f;
            for (int j = 0; j < div; j++) {

                for (int k = 0; k < raysPerPixel; k++) {
                    Ray debugRay;
                    debugRay.inLense = false;
                    camera_get_debug_ray(0.5f, 0.5f, debugRay, k, true);

                    vec3 rayDir = debugRay.direction;

                    Ray camRay;
                    camRay.origin = debugRay.origin;
                    camRay.direction = rayDir;
                    camRay.inLense = false;
                    camRay.rgb = col;
                    int depth = 5;

                    struct BeginEnd {
                        vec3 p1, p2;
                        bool inLense;
                        int numHits;
                        vec3 hit_normal_orig;
                        vec3 hit_normal_dir;
                        int rgb;
                    };

                    int hitCount = 0;
                    BeginEnd rayArr[40];
                    rayArr[0].p1 = camRay.origin;
                    rayArr[0].p2 = camRay.origin + normalize(camRay.direction) * 200;
                    rayArr[0].inLense = false;
                    rayArr[0].numHits = 0;
                    bool no_reflect = false;
                    for (int i = 1; i <= depth; i++) {

                        bool hitPoint = false;

                        world_hit_check(camRay, hitPoint, shadow, no_reflect);

                        rayArr[i - 1].hit_normal_orig = camRay.hit_normal_orig;
                        rayArr[i - 1].hit_normal_dir = camRay.hit_normal_dir;

                        rayArr[i].p1 = camRay.origin;
                        rayArr[i].p2 = camRay.origin + normalize(camRay.direction) * 200;
                        rayArr[i].inLense = camRay.inLense;
                        rayArr[i].rgb = col;


                        if (hitPoint) {

                            rayArr[i - 1].p2 = camRay.origin;
                            hitCount++;

                        }
                        else {
                            rayArr[i - 1].p2 = camRay.origin + normalize(camRay.direction) * 200;

                            break;
                        }
                        rayArr[i].numHits = i;

                    }

                    float t;
                    // color rays
                    for (int i = 0; i <= hitCount; i++) {
                        if (intersection_ray_cylinder(r.origin, r.direction, rayArr[i].hit_normal_orig, (rayArr[i].hit_normal_orig + (normalize(rayArr[i].hit_normal_dir) * 0.1f)), .002f, t)) {
                            if (rayArr[i].hit_normal_orig != vec3(0.0f)) {
                                if (rayArr[i].hit_normal_dir.x > 0.0f)
                                    r.color = vec4(0.0f, 200.0f, 200.0f, 1.0f);
                                else
                                    r.color = vec4(200.0f, 0.0f, 0.0f, 1.0f);

                                return true;
                            }
                        }
                    }

                    // color rays
                    for (int i = 0; i <= hitCount; i++) {
                        if (intersection_ray_cylinder(r.origin, r.direction, rayArr[i].p1, rayArr[i].p2, .02f, t)) {

                            r.color = vec4(0.0f, 200.0f, 200.0f, 1.0f);

                            if (rayArr[i].rgb == 0) {
                                r.color = vec4(255.0f, 0.0f, 0.0f, 1.0f);
                            }
                            if (rayArr[i].rgb == 1) {
                                r.color = vec4(0.0f, 255.0f, 0.0f, 1.0f);
                            }
                            if (rayArr[i].rgb == 2) {
                                r.color = vec4(0.0f, 0.0f, 255.0f, 1.0f);
                            }


                            return true;
                        }
                    }

                }
                valY += partY;
            }
            valX += partX;
        }

    }
    return false;
}

void main()
{   
    int depth = 5; // number of ray bounces
    float imageWidth = widthHeight.x, imageHeight = widthHeight.y;
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;       


    double s = x  / ((imageWidth - 1));
    double t = y  / ((imageHeight - 1));    

    lens_data_construct(lens, lensPos);


    f = calc_f(lens);

    lensArr[0] = lens;

    Ray r;

    r.inLense = false;

    int raysPerPixel = int(factXY.y);
    vec4 col = vec4(0.0f);

    // define sampling parameters
    vec3 circleCenter = vec3(0.0f,0.0f,0.0f), normal = vec3(1.0f,0.0f,0.0f);
    float circleRad = .15f;
    int numCircles = 1, numPointsPerCircle = 1;
    // add sampling parameters
    add_circle_points(circleCenter, circleRad, numCircles, numPointsPerCircle, normal);

    if (renderImageOnly) {

        raysPerPixel = numPoints;

        int loops = 0;
        double dist = 0.0;
        hit_record rec;
        double closest_so_far = infinity;

        Ray ray;
        ray.origin = vec3((initlowerLeftCorner + 0.5f * initHorizontal + 0.5f * initVertical));
        ray.direction = vec3(1.0,0.0,0.0);

        for (int i = 0; i < sphereArr.length(); i++)
        {
            Sphere sphere = sphereArr[i];
            vec3 rayDir = vec3(1.0, 0.0, 0.0);

            if (sphere_hit(0.0, ray, closest_so_far, rec, sphereArr[i])) {

                dist = length(vec3(rec.p - ray.origin));
            }
        }

        for (int i = 0; i < raysPerPixel; i++) {
            int step = i;

            // Debug Parameters
            if (debugMode) {

                camera_get_ray(s, t, r, step, false);

                col += ray_color(r, depth);
                vec3 planePos = initlowerLeftCorner.xyz;
                planePos.x -= 1.0f;
              
                if (rayIntersectsPlane(r.origin, r.direction,planePos , cross(initHorizontal.xyz, initVertical.xyz), initHorizontal.xyz, initVertical.xyz))
                     col = vec4(0.0f, 0.0f, 0.0f, 1.0f);
                

                if (draw_colored_rays(r)) {
                    col = (r.color / numPoints);
                }
                

            }
            else {
                camera_get_ray(s, t, r, step, true);
                Ray initRay;
                float accBrightness = 0.0f;
                for (int color = 0; color < 3; color++) {
                    initRay = r;
                    initRay.rgb = color; 

                    vec4 colorVal = ray_color(initRay, depth);

                    // 3 render passes for chromatic aberration r/g/b
                    if (color == 0) {
                        colorVal.y = 0.0;
                        colorVal.z = 0.0;
                    
                        colorVal.x *= initRay.brightness;
                    }
                    else
                    if (color == 1) {
                        colorVal.x = 0.0;
                        colorVal.z = 0.0;
                    
                        colorVal.y *= initRay.brightness;
                    
                    }
                    else
                    if (color == 2) {
                        colorVal.x = 0.0;
                        colorVal.y = 0.0;
                    
                        colorVal.z *= initRay.brightness;
                    }
                    else {
                        colorVal = vec4(0.0);
                    }                        

                    col += colorVal;
                    notLitPoints = 0;
                }

                col *= 1.0;

                col.w = 1.0;
            }

            loops++;
        }

        col /= loops;
    }
    else {
        raysPerPixel = 1;

        int loops = 0;

        for (int i = 0; i < raysPerPixel; i++) {

            camera_get_ray(s, t, r, i , false);

            col += ray_color(r, depth);

            loops++;
        }

       col /= loops;

    }

    FragColor = col;
};