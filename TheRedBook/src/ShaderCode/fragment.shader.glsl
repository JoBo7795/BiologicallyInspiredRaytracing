#version 460 core

struct Sphere {
    vec3 position;
    float radius;
    vec4 color;
};

layout(std430, binding = 3) buffer dataSSBO
{
    Sphere sphereArr[];
};

in vec3 fPos;
in vec3 camPos;
//in vec4 fVecColor;
//in vec2 fTexCoords;
//in vec3 fNormal;
//in mat4 fMvp;

//uniform vec3 customObjectColor;
//uniform bool textureSet;
//uniform vec3 dir;

out vec4 FragColor;

const float infinity = 1. / 0.;
const float M_PI = 3.14159265359;

struct Ray {
    vec3 origin,direction;
    vec4 color;

};



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

    return vec3(rand(), rand(), rand());
}

vec3 randomVec3(float min, float max) {
    return vec3(rand(min, max), rand(min, max), rand(min, max));
}

vec3 random_in_unit_sphere() {
    while (true) {
        vec3 p = randomVec3(-1, 1);
        if ((length(p) * length(p)) >= 1) continue;
        return p;
    }
}

bool hitSphere_2(Ray ray, Sphere sphere, inout float t0, inout float t1, vec3 t_out)
{

    float radius2 = sphere.radius * sphere.radius;

    vec3 l = sphere.position - ray.origin;
    float tca = dot(l,ray.direction);
    if (tca < 0) return false;
    float d2 = dot(l,l) - tca * tca;
    if (d2 > radius2) return false;
    float thc = sqrt(radius2 - d2);
    t0 = tca - thc;
    t1 = tca + thc;

    t_out = ray.direction * t0;

    if (t0 < 0.0f)
        t_out = ray.direction * t1;

    return true;
}

vec3 ReflectionRay(Sphere sphere, Ray ray, vec3 collPoint) {

    vec3 normal = normalize((sphere.position - collPoint) / sphere.radius);

    return reflect(ray.direction,normal);
}

vec3 Refract(vec3 uv, vec3 n, float eta) {
    return refract(uv,n,eta);
}

vec4 RecursiveTrace(Ray ray,Sphere SphereArr[4],int SphereArrSize, int depth) {

    vec4 resColor;
    bool hit = false;
    vec3 t_out;
    float t0, t1;

    for (int j = 0; j <= depth; j++) {
        for (int i = 0; i < SphereArrSize; i++) {
            if (hitSphere_2(ray, SphereArr[i], t0, t1, t_out)) {

                vec3 reflection;

                if (i == 1) {
                    vec3 normal = normalize((SphereArr[i].position - ray.direction * t0) / SphereArr[i].radius);
                    reflection = Refract(normalize(ray.direction * t0), normalize(normal), .5f);
                }
                else {
                    reflection = ReflectionRay(SphereArr[i], ray, ray.direction * t0);
                }  

                ray.origin = ray.direction * t0;
                ray.direction = reflection;
                resColor = SphereArr[i].color;

                hit = true;
            }
        }
    }

    if (!hit) {
        resColor = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    }

    return resColor;

}



void main()
{	

    Sphere SphereArr[4];
    int SphereArrSize = 4;
    int depth = 1;
    float imageWidth = 1920, imageHeight = 1080;
    float x = gl_FragCoord.x;
    float y = gl_FragCoord.y;
    float fov = 90;
    
    
    vec3 rayOrigin = camPos;//vec3(0.0f,0.0f,0.0f);

    Sphere sphere;
    for (int i = 0; i < 3; i++) {

        SphereArr[i] = sphereArr[i];
    }

    /*Sphere sphere1, sphere2, sphere3, sphere4;
    sphere1.position = (vec4(0.0f, 0.0f, 30.0f, 1.0f)).xyz;       
    sphere1.radius = 1.0f;
    sphere1.color = vec4(1.0f, 1.0f, 0.0f, 1.0f);

    sphere2.position = (vec4(20.0f, 0.0f, 30.0f, 1.0f)).xyz;
    sphere2.radius = 1.0f;
    sphere2.color = vec4(0.5f, 0.7f, 0.3f, 1.0f);

    sphere3.position = (vec4(-10.0f, 0.0f, 30.0f, 1.0f)).xyz;
    sphere3.radius = 1.0f;
    sphere3.color = vec4(1.0f, 0.0f, 1.0f, 1.0f);

    sphere4.position = (vec4(10.0f, 10.0f, -30.0f, 1.0f)).xyz;
    sphere4.radius = 5.0f;
    sphere4.color = vec4(0.0f, 0.0f, 1.0f, 1.0f);

    SphereArr[0] = sphere1;
    SphereArr[1] = sphere2;
    SphereArr[2] = sphere3;
    SphereArr[3] = sphere4;*/

    Ray ray; 

    float imageAspectRatio = imageWidth / imageHeight;
    float Px = (2 * ((x + 0.5) / imageWidth) - 1) * tan(fov / 2 * M_PI / 180) * imageAspectRatio;
    float Py = (1 - 2 * ((y + 0.5) / imageHeight) * tan(fov / 2 * M_PI / 180));


    ray.origin = rayOrigin;
    ray.direction = normalize(vec3(vec4((vec3(Px, Py, 1) - ray.origin),1.0f)));

    FragColor = vec4((vec3(0.0f, 0.0f, 0.0f) + RecursiveTrace(ray, SphereArr, 3, depth).xyz),1.0f);
};
