#include "CPUMode.h"

//int CPUMode::image_width, CPUMode::image_height;
int CPUMode::depth = 5;
const float CPUMode::infinity = std::numeric_limits<float>::infinity();
const float CPUMode::M_PI = 3.14159265359;
float CPUMode::global_refract_index = 1.59;
int CPUMode::numPoints = 0;
uint8_t* CPUMode::pixels = new uint8_t[ProgramParams::windowWidth * ProgramParams::windowHeight * CHANNEL_NUM];
uint8_t* CPUMode::fBuffer = new uint8_t[ProgramParams::windowWidth * ProgramParams::windowHeight * CHANNEL_NUM];

glm::vec3 CPUMode::lowerLeftCorner, CPUMode::horizontal, CPUMode::vertical, CPUMode::lookFrom;
glm::vec3 CPUMode::initlowerLeftCorner, CPUMode::initHorizontal, CPUMode::initVertical, CPUMode::initCamPos;
std::vector<glm::vec3> CPUMode::circlePoints;
Camera CPUMode::debugCam;
std::ofstream CPUMode::dataText = std::ofstream("C:\\Users\\Johannes\\Desktop\\Renderings\\datei.txt");
LensData CPUMode::lens;
std::map<float, float> CPUMode::lambdaVal =
{
{400,	0.000396},
{410,	0.00121},
{420,	0.00400},
{430,	0.0116 },
{440,	0.0230 },
{450,	0.0380 },
{460,	0.0600 },
{470,	0.0910 },
{480,	0.139  },
{490,	0.208  },
{500,	0.323  },
{510,	0.503  },
{520,	0.710  },
{530,	0.862  },
{540,	0.954  },
{550,	0.995  },
{560,	0.995  },
{570,	0.952  },
{580,	0.870  },
{590,	0.757  },
{600,	0.631  },
{610,	0.503  },
{620,	0.381  },
{630,	0.265  },
{640,	0.175  },
{650,	0.107  },
{660,	0.061  },
{670,	0.032  },
{680,	0.017  },
{690,	0.0082 },
},




CPUMode::lambdaPrimeVal =
{
{400,	0.00929   },
{410,	0.03484   },
{420,	0.0966    },
{430,	0.1998    },
{440,	0.3281    },
{450,	0.4550    },
{460,	0.5670    },
{470,	0.6760    },
{480,	0.7930    },
{490,	0.9040    },
{500,	0.9820    },
{510,	0.9970    },
{520,	0.9350    },
{530,	0.8110    },
{540,	0.6500    },
{550,	0.4810    },
{560,	0.3288    },
{570,	0.2076    },
{580,	0.1212    },
{590,	0.0655    },
{600,	0.03315   },
{610,	0.01593   },
{620,	0.00737   },
{630,	0.003335  },
{640,	0.001497  },
{650,	0.000677  },
{660,	0.0003129 },
{670,	0.0001480 },
{680,	0.0000715 },
{690,	0.00003533},
};

int notLitPoints = 0;
float iris_radius = 2.05f;

void CPUMode::InitCPUMode(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir)
{

    Camera* cam = Renderer::GetInstance()->GetCamera();  

    lowerLeftCorner = cam->lower_left_corner;
    horizontal = cam->horizontal;
    vertical = cam->vertical;
    lookFrom = cam->GetPosition();   
}


void CPUMode::SetDebugParams(glm::vec3 in_imagePlanePos, glm::vec3 in_imagePlaneDir)
{

    debugCam = Camera(in_imagePlanePos);

    debugCam.SetPosition(in_imagePlanePos);
    debugCam.SetDirection(in_imagePlaneDir);
    debugCam.InitRTCamera();
    debugCam.Update();


    initHorizontal = debugCam.horizontal;
    initVertical = debugCam.vertical;
    initlowerLeftCorner = debugCam.lower_left_corner;
    initCamPos = debugCam.GetPosition();
}

void CPUMode::Render() {

    lens = LensData(debugCam.GetPosition());
    ProgramParams::f = lens.calcF();
    ProgramParams::dpt = lens.calcOneByF();
    ProgramParams::g = 10.0;
    ProgramParams::calcImageDist();

    ProgramParams::g = 21.0;
    double dval = lens.calcD();
    lens.d = dval;     
    
    lens.dTolensRad();

    std::cout << "f: " << ProgramParams::f << " dpt: " << ProgramParams::dpt << " dval: " << dval << std::endl;

    glm::vec3 circleCenter = glm::vec3(0.0f, 0.0f, 0.0f), normal = glm::vec3(1.0f, 0.0f, 0.0f);
    
    float circleRad = .15f;
    int numCircles =1 , numPointsPerCircle = 1;    

    addCirclePoints(circleCenter, circleRad, numCircles, numPointsPerCircle, normal);
    
    double picBrightness = 0.0;

    std::cout << "P3\n" << ProgramParams::windowWidth << ' ' << ProgramParams::windowHeight << "\n255\n";
    int index = 0;
    for (int j = ProgramParams::windowHeight - 1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < ProgramParams::windowWidth; ++i) {

            glm::vec3 pixel_color(0, 0, 0);

            double u = float((i) / (ProgramParams::windowWidth - 1.0f));
            double v = float((j) / (ProgramParams::windowHeight - 1.0f));
            Ray r;
            r.id = j * ProgramParams::windowHeight + i;

            auto sphereQueue = DrawObjectManager::GetSphereQueue();
            glm::vec3 rayPos = glm::vec3((initlowerLeftCorner + 0.5f * initHorizontal + 0.5f * initVertical));
            rayPos.x += 1.0;

            rayPos = glm::vec3(0.0, 0.0, 0.0);

            double closestSoFar = CPUMode::infinity;
            double dist = 0.0;
            for (size_t i = 0; i < sphereQueue.size(); i++)
            {
                Sphere sphere = sphereQueue[i];
                glm::vec3 rayDir = glm::vec3(1.0, 0.0, 0.0);
                
                if (CollisionDetection::RaySphereIntersection(sphere.position, sphere.radius, rayPos, rayDir, 0.0, closestSoFar)) {
                    dist = glm::length(glm::vec3(sphere.position-rayPos));
                }
            }

            if (ProgramParams::debugMode) {

                debugCam.SetPosition(glm::vec3(-ProgramParams::b, ProgramParams::posY, ProgramParams::posZ));
                initlowerLeftCorner.x += 1.0;

                debugCam.SetPosition(lens.lensOrigin + glm::normalize(lens.lensOriginLeft) * glm::vec3(ProgramParams::b));                

            }

            
            if (ProgramParams::debugMode) {

                CameraGetRay(u, v, r, 0, false);

                pixel_color += glm::normalize(glm::vec3(ray_color(r, depth)));

                for (int i = 0; i < 3; i++) {
                    r.rgb = i;
                    if (drawColoredRays(r)) {

                        pixel_color = glm::vec3((r.color / float(numPoints)));

                        int brightnessVal = notLitPoints / numPoints;

                        glm::vec3 brightness = glm::vec3(brightnessVal);

                        if (brightnessVal > 0.3) {
                            pixel_color = (pixel_color + ( brightness)) / glm::vec3(2);
                        }

                    }
                }

                pixels[index++] = static_cast<int>(255.999 * (pixel_color.x)/numPoints);
                pixels[index++] = static_cast<int>(255.999 * (pixel_color.y)/numPoints);
                pixels[index++] = static_cast<int>(255.999 * (pixel_color.z)/numPoints);
            }
            else {
                float accBrightness = 0.0;
                
                for (int i = 0; i < numPoints; i++) {
                    CameraGetRay(u, v, r, i, true);
                    Ray initRay;
                    for (int colorCounter = 0; colorCounter < 3; colorCounter++) {
                        initRay = r;
                        initRay.rgb = colorCounter;
                        initRay.rgb = 4;

                        glm::vec3 colorVal = glm::vec3(ray_color(initRay, depth));
                        float brightnessVal = 0.3f;
                        glm::vec3 brightness = glm::vec3(brightnessVal);

                        float brightnessThreshold = .1;
                        accBrightness += initRay.brightness;
                        if (brightnessVal <= brightnessThreshold) {
                            colorVal = glm::vec4(1.0) * initRay.brightness;
                        }
                        else {

                            if (colorCounter == 0) {
                                colorVal.y = 0.0;
                                colorVal.z = 0.0;

                                colorVal.x *= initRay.brightness;
                            }
                           else
                           if (colorCounter == 1) {
                               colorVal.x = 0.0;
                               colorVal.z = 0.0;

                               colorVal.y *= initRay.brightness;

                           }   
                           else
                           if (colorCounter == 2) {
                               colorVal.x = 0.0;
                               colorVal.y = 0.0;

                               colorVal.z *= initRay.brightness;

                           }
                           else {
                               colorVal = glm::vec3(0.0);
                           }

                        }

                        pixel_color += colorVal;     
                        
                        std::cout << "";
                        notLitPoints = 0;
                    }
                    
                }

                accBrightness /= (numPoints*3);

                float pixCx = (pixel_color.x / numPoints),
                    pixCy = (pixel_color.y / numPoints),
                    pixCz = (pixel_color.z / numPoints);

                pixels[index++] = static_cast<int>(255.999 * pixCx);
                pixels[index++] = static_cast<int>(255.999 * pixCy);
                pixels[index++] = static_cast<int>(255.999 * pixCz);

                picBrightness += ((pixCx + pixCy + pixCz) / 3) <= 0 ? 0 : (((pixCx + pixCy + pixCz) / 3));
            }


            fBuffer = pixels;



        }
    }

    picBrightness /= (ProgramParams::windowWidth * ProgramParams::windowHeight);
    std::cout << "brightness: " << picBrightness << std::endl;


}

float CPUMode::lerp(float x, float y, float t) {
    return x * (1.f - t) + y * t;
}

double CPUMode::logInterpolation(double x1, double y1, double x2, double y2, double x) {
    if (y1 <= 0 || y2 <= 0) {
        // Logarithmische Interpolation funktioniert nicht für y-Werte <= 0
        return -1;
    }

    double slope = (log(y2) - log(y1)) / (x2 - x1);
    double intercept = log(y1) - slope * x1;

    return exp(slope * x + intercept);
}

struct Disc {
    glm::vec3 origin, normal;
    double radius;
};

bool doesRayIntersectDisc(Ray ray, Disc disc) {

    double dotProduct = ray.direction.x * disc.normal.x + ray.direction.y * disc.normal.y + ray.direction.z * disc.normal.z;

    // if dot product is 0 ray and disc are parallel and don't intersect
    if (dotProduct == 0) {
        return false;
    }

    // vector from ray origin to disc origin
    glm::vec3 originToDisc = { disc.origin.x - ray.origin.x, disc.origin.y - ray.origin.y, disc.origin.z - ray.origin.z };

    // distance of ray and disc origin
    double t = (originToDisc.x * disc.normal.x + originToDisc.y * disc.normal.y + originToDisc.z * disc.normal.z) / dotProduct;

    // if t negative, disc is behind camera
    if (t < 0) {
        return false;
    }

    // calculate intersection point
    glm::vec3 intersection = { ray.origin.x + t * ray.direction.x, ray.origin.y + t * ray.direction.y, ray.origin.z + t * ray.direction.z };


    double distanceSquared = (intersection.x - disc.origin.x) * (intersection.x - disc.origin.x) +
        (intersection.y - disc.origin.y) * (intersection.y - disc.origin.y) +
        (intersection.z - disc.origin.z) * (intersection.z - disc.origin.z);

    // if distance squared bigger than quadratic radius of disc, intersection located outside of disc
    if (distanceSquared > (disc.radius * disc.radius)) {
        return true;
    }

    return false;
}


bool CPUMode::WorldHit(float t_min, Ray& r, float t_max, hit_record& rec) {

    Disc iris;
    iris.origin = glm::vec3(0.0);
    iris.normal = glm::vec3(1.0, 0.0, 0.0);
    
    iris.radius = iris_radius;

    if (doesRayIntersectDisc(r, iris)) {
    
        rec.geoType = IRIS;
    
        return true;
    
    }

    hit_record temp_rec;
    bool hit_anything = false;
    float closest_so_far = t_max;
    glm::vec3 p;
    auto sphereQueue = DrawObjectManager::GetSphereQueue();
    for (int i = 0; i < DrawObjectManager::GetSphereQueue().size(); i++) {
        if (SphereHit(t_min, r, closest_so_far, temp_rec, DrawObjectManager::GetSphereQueue()[i])) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
            rec.arrIndex = i;
            rec.geoType = SPHERE;
    
        }
    }

    double t = 0.0;
    glm::vec3 normal_out;
    float x = 1.0;
    float y = 4.0;
    float z = 4.0;
    glm::vec3 ellipsoidRadii = glm::vec3(x, y, z);


     if (IntersectRayLens(r.origin, r.direction, lens, normal_out, t, r.inLense)) {
          if (closest_so_far > t) {
              hit_anything = true;
              rec = temp_rec;
     
              rec.t = t;
              closest_so_far = t;
     
              rec.p = ray_at(r.origin, r.direction, t);
     
              glm::vec3 outward_normal = normal_out;
              rec.front_face = check_front_face(r, outward_normal);
              rec.normal = set_face_normal(r, outward_normal);
              rec.mat_ptr = DIELECTRIC;
              rec.arrIndex = 3;
              rec.geoType = SPHERE;
              r.inLense = !r.inLense;
     
     
              r.hit_normal_dir = glm::normalize(rec.normal);
              r.hit_normal_orig = rec.p;
          }
     }

     


    Triangle tri;
    Sphere boundingSphere;
    glm::vec3 boundingSpherePos;
    int offsetI = 0;
    auto goList = GameObjectManager::GetGameObjectList();    


    for (int o = 0; o < goList.size(); o++) {
        break;

        auto dO = DrawObjectManager::GetDrawObjectByID(goList[o].GetDrawObjectID());

        for (int i = int(dO.IndiBufferRangeL);
            i < (int(dO.IndiBufferRangeR)); i += 3) {
           
            auto translation = goList[o].GetTranslation();

            if (o <= 0) {
                tri.position[0] = translation * glm::vec4(dO.shapeData.vertices[dO.shapeData.indices[i] * 6], dO.shapeData.vertices[dO.shapeData.indices[i] * 6 + 1], dO.shapeData.vertices[dO.shapeData.indices[i] * 6 + 2], 1.0f);
                tri.position[1] = translation * glm::vec4(dO.shapeData.vertices[dO.shapeData.indices[i + 1] * 6], dO.shapeData.vertices[dO.shapeData.indices[i + 1] * 6 + 1], dO.shapeData.vertices[dO.shapeData.indices[i + 1] * 6 + 2], 1.0f);
                tri.position[2] = translation * glm::vec4(dO.shapeData.vertices[dO.shapeData.indices[i + 2] * 6], dO.shapeData.vertices[dO.shapeData.indices[i + 2] * 6 + 1], dO.shapeData.vertices[dO.shapeData.indices[i + 2] * 6 + 2], 1.0f);
               
           
            }
            else {
              //  tri.position[0] = translation * glm::vec4(dO.shapeData.vertices[dO.shapeData.indices[i * 3]], dO.shapeData.vertices[dO.shapeData.indices[i * 3] + 1], dO.shapeData.vertices[dO.shapeData.indices[i * 3] + 2], 1.0f);
              //  tri.position[1] = translation * glm::vec4(dO.shapeData.vertices[dO.shapeData.indices[i * 3 + 1]], dO.shapeData.vertices[dO.shapeData.indices[i * 3 + 1] + 1], dO.shapeData.vertices[dO.shapeData.indices[i * 3 + 1] + 2], 1.0f);
              //  tri.position[2] = translation * glm::vec4(dO.shapeData.vertices[dO.shapeData.indices[i * 3 + 2]], dO.shapeData.vertices[dO.shapeData.indices[i * 3 + 2] + 1], dO.shapeData.vertices[dO.shapeData.indices[i * 3 + 2] + 2], 1.0f);
            }

            tri.albedo = glm::vec4(goList[o].GetColor(),1.0f);
            tri.fuzzAndmat_ptr = glm::vec4(1.0f,1.0f,0.0f,0.0f);




            if (RayTriangleIntersect(r, temp_rec, t_min, closest_so_far, tri, p)) {
               
              
                glm::vec3 tmpNorm1 = (dO.shapeData.normals[dO.shapeData.indices[i]]),
                    tmpNorm2 = (dO.shapeData.normals[dO.shapeData.indices[i + 1]]),
                    tmpNorm3 = (dO.shapeData.normals[dO.shapeData.indices[i + 2]]);

                auto offset = 9;

                glm::vec3 tmpCol0 = (glm::vec3(dO.shapeData.vertices[dO.shapeData.indices[i] * offset + 7], dO.shapeData.vertices[dO.shapeData.indices[i] * offset + 7 +1], dO.shapeData.vertices[dO.shapeData.indices[i] * offset + 7 + 2])),
                          tmpCol1 = (glm::vec3(dO.shapeData.vertices[dO.shapeData.indices[i + 1] * offset+ 7], dO.shapeData.vertices[dO.shapeData.indices[i + 1] *offset+ 7 + 1], dO.shapeData.vertices[dO.shapeData.indices[i + 1] * offset+ 7 + 2])),
                          tmpCol2 = (glm::vec3(dO.shapeData.vertices[dO.shapeData.indices[i + 2] * offset+ 7], dO.shapeData.vertices[dO.shapeData.indices[i + 2] *offset+ 7 + 1], dO.shapeData.vertices[dO.shapeData.indices[i + 2] * offset+ 7 + 2]));

                tri.normal = glm::vec4(triInterpolNormal(
                    tri.position[0], tri.position[1], tri.position[2], tmpNorm1, tmpNorm2, tmpNorm3, p), 1.0f);

                tri.albedo = glm::vec4(triInterpolNormal(
                    tri.position[0], tri.position[1], tri.position[2], tmpCol0, tmpCol1, tmpCol2, p), 1.0f);

                std::cout << glm::to_string(tmpCol0) << std::endl << glm::to_string(tmpCol1) << std::endl << glm::to_string(tmpCol2) << std::endl;
                                
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
                rec.arrIndex = i;
                rec.geoType = TRIANGLE;
                rec.tri = tri;

                rec.tri.fuzzAndmat_ptr.y = DIELECTRIC;
                rec.front_face = check_front_face(r, tri.normal);
                rec.normal = set_face_normal(r, tri.normal);

                r.hit_normal_dir = set_face_normal(r, tri.normal);
                r.hit_normal_orig = glm::vec3(r.direction * float(temp_rec.t));

            }

        }
    }


    return hit_anything;
}

glm::vec3 CPUMode::triInterpolNormal(glm::vec3 tri0, glm::vec3 tri1, glm::vec3 tri2, glm::vec3 normal0, glm::vec3 normal1, glm::vec3 normal2, glm::vec3 vertex) {

    glm::vec3 C;
    float u_n, v_n, w_n;

    glm::vec3 v0v1 = tri1 - tri0;
    glm::vec3 v0v2 = tri2 - tri0;

    glm::vec3 N = glm::cross(v0v1, v0v2);

    glm::vec3 edge1 = tri2 - tri1;
    glm::vec3 vp1 = vertex - tri1;

    C = glm::cross(edge1, vp1);

    u_n = glm::dot(N, C);

    glm::vec3 edge2 = tri0 - tri2;
    glm::vec3 vp2 = vertex - tri2;

    C = glm::cross(edge2, vp2);

    v_n = glm::dot(N, C);

    float denom = glm::dot(N, N);

    u_n /= denom;
    v_n /= denom;
    w_n = 1.0f - u_n - v_n;


    return glm::normalize(u_n * normal0 + v_n * normal1 + w_n * normal2);

}

bool CPUMode::RayTriangleIntersect(
    Ray r,
    hit_record& rec,
    float t_min,
    float t_max,
    Triangle tri,
    glm::vec3& outIntersectionPoint)
{


    glm::vec3 vertex0 = tri.position[0];
    glm::vec3 vertex1 = tri.position[1];
    glm::vec3 vertex2 = tri.position[2];

    glm::vec4 normal = tri.normal;

    glm::vec3 edge1, edge2, h, s, q;

    float a, f, u, v;

    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;

    h = glm::cross(r.direction, edge2);
    a = glm::dot(edge1, h);

    if (a > -EPSILON && a < EPSILON)
        return false;    // This ray is parallel to this triangle.

    f = 1.0 / a;
    s = r.origin - vertex0;
    u = f * glm::dot(s, h);
    if (u < 0.0 || u > 1.0)
        return false;
    q = glm::cross(s, edge1);
    v = f * glm::dot(r.direction, q);
    if (v < 0.0 || u + v > 1.0)
        return false;
    // At this stage we can compute t to find out where the intersection point is on the line.
    float t = f * glm::dot(edge2, q);

    if (t < t_min || t_max < t)
        return false;

    rec.t = t;
    rec.p = ray_at(r.origin, r.direction, rec.t);
    
    glm::vec3 outward_normal = glm::normalize(glm::cross(edge1, edge2));
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


bool CPUMode::IntersectRaySphere(glm::vec3 rayOrigin, glm::vec3 rayDirection, glm::vec3 sphereOrigin, double sphereRadius, glm::vec3& iPoint, double& t) {
    glm::vec3 oc = rayOrigin - sphereOrigin;
    double a = glm::dot(rayDirection, rayDirection);
    double b = 2.0 * glm::dot(oc, rayDirection);
    double c = glm::dot(oc, oc) - sphereRadius * sphereRadius;
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return false;
    }
    else {
        double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
        double t2 = (-b + sqrt(discriminant)) / (2.0 * a);
        if (t1 < 0 && t2 < 0) {
            // Beide Intersektionspunkte liegen in der negativen Richtung des Strahls
            return false;
        }
        else if (t1 > 0 && t2 > 0) {
            // Beide Intersektionspunkte liegen in der positiven Richtung des Strahls
            // Wir nehmen den kleineren Wert von t
            iPoint = rayOrigin + float(std::min(t1, t2)) * rayDirection;
            t = std::min(t1, t2);
        }
        else {
            // Ein Intersektionspunkt liegt in der positiven Richtung des Strahls
            // Wir nehmen den größeren Wert von t
            iPoint = rayOrigin + float(std::max(t1, t2)) * rayDirection;
            t = std::max(t1,t2);
        }
        return true;
    }
}



bool CPUMode::IntersectRayLens(glm::vec3 rayOrigin, glm::vec3 rayDirection, LensData lens, glm::vec3& normal, double& t, bool inLense) {
   
    glm::vec3 intersection1;
    glm::vec3 intersection2;
    double t1 = 0.0, t2 = 0.0;

    bool inter1 = IntersectRaySphere(rayOrigin, normalize(rayDirection), lens.lensOriginLeft, lens.lensRadiusLeft, intersection1,t1);
    bool inter2 = IntersectRaySphere(rayOrigin, normalize(rayDirection), lens.lensOriginRight, lens.lensRadiusRight, intersection2,t2);

    auto angle1 = glm::degrees(glm::acos(glm::dot(glm::normalize((lens.lensOrigin - lens.lensOriginRight)), glm::normalize(lens.lensOrigin - intersection1))));
    auto angle2 = glm::degrees(glm::acos(glm::dot(glm::normalize((lens.lensOrigin - lens.lensOriginLeft)), glm::normalize(lens.lensOrigin - intersection2))));    

    if (inter2 && (angle1 <= 90) && inLense) {

        t = t1;
        normal = glm::normalize(glm::vec3(intersection1 - lens.lensOriginLeft));  

        return true;
    }

    if (inter2 && (angle2 <= 90) && !inLense) {

        t = t2;
        normal = glm::normalize(glm::vec3(intersection2 - lens.lensOriginRight));

        return true;
    }

    return false;
}




glm::vec4 CPUMode::ray_color2(Ray& r, bool& hit, bool& no_reflect, bool& no_shadow) {

    hit_record rec;
    Ray scattered;
    glm::vec4 attenuation;

    if (WorldHit(0.01f, r, infinity, rec)) {

        if (rec.geoType == SPHERE) {

            switch (int(DrawObjectManager::GetSphereQueue()[rec.arrIndex].mat_ptr)) {
            case int(LAMBERTIAN):

                if (lambertian_scatter(r, rec, attenuation, scattered, DrawObjectManager::GetSphereQueue()[rec.arrIndex])) {

                    hit = true;

                    r = scattered;

                    return attenuation;


                }


                break;

            case int(DIELECTRIC):
                if (DielectricScatter(r, rec, attenuation, scattered, DrawObjectManager::GetSphereQueue()[rec.arrIndex])) {

                    hit = true;

                    scattered.id = r.id;
                    scattered.rgb = r.rgb;
                    r = scattered;

                    no_shadow = true;

                    return attenuation;


                }

                break;

            case int(METAL):
                if (metal_scatter(r, rec, attenuation, scattered, DrawObjectManager::GetSphereQueue()[rec.arrIndex])) {

                    hit = true;

                    r = scattered;

                    return attenuation;

                }
                break;

            case int(PLAIN):
                if (plain_color(r, rec, attenuation, scattered, DrawObjectManager::GetSphereQueue()[rec.arrIndex])) {

                    hit = true;
                    no_reflect = true;

                    r = scattered;
                    r.brightness = 1.;
                    return attenuation;

                }
                break;
            default:
                glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);               
                break;
            }
        }
        else if (rec.geoType == LENSE) {
            Sphere sphere;
            sphere.fuzz = 1.0f;
            sphere.mat_ptr = DIELECTRIC;
            sphere.color = glm::vec4(1.0f, 0.0f, 0.0f, 1.0f);
            sphere.albedo = glm::vec4(0.0f, 0.0f, 1.0f, 1.0f);

            if (DielectricScatter(r, rec, attenuation, scattered, sphere)) {

                hit = true;

                r = scattered;

                no_shadow = true;

                return attenuation;
                
            }
                       
        }
        else if (rec.geoType == IRIS) {

            hit = true;

            notLitPoints++;

            return glm::vec4(0.0);
        }
        else {
            switch (int(rec.tri.fuzzAndmat_ptr.y)) {
            case int(LAMBERTIAN):

                if (lambertian_scatter(r, rec, attenuation, scattered, rec.tri)) {

                    hit = true;

                    r = scattered;

                    return attenuation;

                }

                break;

            case int(DIELECTRIC):
                if (DielectricScatter(r, rec, attenuation, scattered, rec.tri)) {

                    hit = true;                  

                    r = scattered;

                    no_shadow = true;

                    return attenuation;


                }


                break;

            case int(METAL):
                if (metal_scatter(r, rec, attenuation, scattered, rec.tri)) {

                    hit = true;

                    r = scattered;

                    return attenuation;

                }
                break;
            case int(PLAIN):
                if (plain_color(r, rec, attenuation, scattered, rec.tri)) {

                    hit = true;
                    no_reflect = true;

                    r = scattered;

                    return attenuation;

                }
                break;
            default:
                glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);               
                break;
            }
        }



    }




    no_shadow = true;
    hit = false;
    glm::vec3 unit_direction = glm::normalize(r.direction);
    float t = 0.5 * (unit_direction.y + 1.0);
    r.brightness = 1.;    
    return glm::vec4((1.0f - t) * glm::vec4(1.0f, 1.0f, 1.0f, 1.0f) + t * glm::vec4(0.5f, 0.7f, 1.0f, 1.0f));

}

glm::vec4 CPUMode::ray_color(Ray& r, int depth) {

    glm::vec4 fColor = glm::vec4(1.0f);
    bool hit = false;
    glm::vec4 color;


    for (int i = depth - 1; i >= 0; i--) {
        r.currDepth = depth;
        bool no_reflect = false;
        bool no_shadow = false;
        color = ray_color2(r, hit, no_reflect, no_shadow);       
        fColor *= color;


        if (no_reflect) {
            break;
        }
    }

    if (hit) {
        return fColor;
    }

    glm::vec3 unit_direction = glm::normalize(r.direction);
    float t = 0.5f * (unit_direction.y + 1.0f);
    r.brightness = 1.;    
    return (glm::vec4((1.0f - t) * glm::vec4(1.0f, 1.0f, 1.0f, 1.0f) + t * glm::vec4(0.5f, 0.7f, 1.0f, 1.0f)));
}


bool CPUMode::CameraGetRay(float s, float t, Ray& ray, int o, bool lp) {

    if (lp) {

        ray.origin = glm::vec3((lowerLeftCorner + s * horizontal + t * vertical));
        ray.origin.x -= 1.0;
        ray.direction = glm::normalize(glm::vec3((circlePoints[o])) - ray.origin);
    }
    else {

        //Perspective Rays
        ray.origin = lookFrom;
        ray.direction = glm::normalize(glm::vec3((lowerLeftCorner + s * horizontal + t * vertical) - lookFrom));

    }

    return true;
}

bool CPUMode::camera_get_debug_ray(float s, float t, Ray& ray, int o, bool lp) {

    if (lp) {
       
       glm::vec3 origin;
      
       ray.origin = debugCam.GetPosition();
        
       if(o == 0)
           ray.direction = glm::normalize(glm::vec3(glm::vec3(0.0,.15,0.0)) - ray.origin);
       if (o == 1)
           ray.direction = glm::normalize(glm::vec3(glm::vec3(0.0, 0.0, 0.0)) - ray.origin);
       if (o == 2)
           ray.direction = glm::normalize(glm::vec3(glm::vec3(0.0, -0.15, 0.0)) - ray.origin);
       if (o == 3)
           ray.direction = glm::normalize(glm::vec3(glm::vec3(0.0, 0.30, 0.0)) - ray.origin);
       if (o == 4)
           ray.direction = glm::normalize(glm::vec3(glm::vec3(0.0, 1.0, 0.0)) - ray.origin);
       if (o == 5)
           ray.direction = glm::normalize(glm::vec3(glm::vec3(0.0, -0.30, 0.0)) - ray.origin);     

    }
    else {
        //Perspective Rays
        ray.origin = initCamPos;
        ray.direction = glm::normalize(glm::vec3((initlowerLeftCorner + s * initHorizontal + t * initVertical) - initCamPos));
    }


    return true;
}

glm::vec3 CPUMode::ray_at(glm::vec3 origin, glm::vec3 direction, float t) {
    return glm::vec3(origin + t * direction);
}

glm::vec3 CPUMode::refract_s(glm::vec3& uv, glm::vec3 n, float etai_over_etat) {
    float cos_theta = glm::min(glm::dot(-uv, n), 1.0f);
    glm::vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
    glm::vec3 r_out_parallel = -sqrt(glm::abs(1.0f - r_out_perp.length() * r_out_perp.length())) * n;
    return r_out_perp + r_out_parallel;

}

bool CPUMode::check_front_face(const Ray r, const glm::vec3 outward_normal) {
    return glm::dot(r.direction, outward_normal) < 0;
}

glm::vec3 CPUMode::set_face_normal(const Ray r, const glm::vec3 outward_normal) {
    return (check_front_face(r, outward_normal) ? outward_normal : -outward_normal);
}

bool CPUMode::SphereHit(float t_min, Ray& r, float t_max, hit_record& rec, Sphere& sphere) {

    glm::vec3 oc = r.origin - sphere.position;
    float a = length(r.direction) * length(r.direction);
    float half_b = glm::dot(oc, r.direction);
    float c = length(oc) * length(oc) - sphere.radius * sphere.radius;

    float discriminant = half_b * half_b - a * c;
    if (discriminant < 0)
        return false;
    float sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    float root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = CPUMode::ray_at(r.origin, r.direction, rec.t);

    glm::vec3 outward_normal = (rec.p - sphere.position) / sphere.radius;
    rec.front_face = check_front_face(r, outward_normal);
    rec.normal = set_face_normal(r, outward_normal);
    rec.mat_ptr = int(sphere.fuzz);

    return true;
}


bool CPUMode::lambertian_scatter(
    Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Sphere& sphere)
{

    glm::vec3 scatter_direction = rec.normal + random_unit_vector();

    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = sphere.albedo;
    return true;
}

bool CPUMode::lambertian_scatter(
    Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Triangle& tri)
{

    glm::vec3 scatter_direction = rec.normal + random_unit_vector();

    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = tri.albedo;
    return true;
}

bool CPUMode::plain_color(
    Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Sphere& sphere)
{
    glm::vec3 scatter_direction = rec.normal;

    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = sphere.albedo;
    return true;
}

bool CPUMode::plain_color(
    Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Triangle& tri)
{

    glm::vec3 scatter_direction = rec.normal;

    scattered.origin = rec.p;
    scattered.direction = scatter_direction;

    attenuation = tri.albedo;
    return true;
}

bool CPUMode::metal_scatter(Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Sphere& sphere) {
    glm::vec3 reflected = reflect(glm::normalize(r_in.direction), rec.normal);
    scattered.origin = rec.p;
    scattered.direction = reflected;
    attenuation = sphere.albedo;
    return (glm::dot(scattered.direction, rec.normal) > 0);
}

bool CPUMode::metal_scatter(Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Triangle& tri) {
    glm::vec3 reflected = reflect(glm::normalize(r_in.direction), rec.normal);
    scattered.origin = rec.p;
    scattered.direction = reflected;
    attenuation = tri.albedo;
    return (glm::dot(scattered.direction, rec.normal) > 0);
}


bool CPUMode::DielectricScatter(
    Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Sphere& sphere
) {
    attenuation = glm::vec4(1.0, 1.0, 1.0, 1.0);
    float refraction_ratio;
   


    glm::vec3 unit_direction = glm::normalize(r_in.direction);

    float cos_theta = glm::min(glm::dot(-unit_direction, rec.normal), 1.0f);
    float sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    glm::vec3 direction;
   
    switch (r_in.rgb)
    {
    case 0:
        refraction_ratio = 1.57;
        break;
    case 1:
        refraction_ratio = 1.59;
        break;
    case 2:
        refraction_ratio = 1.61;
        break;


    default:
        refraction_ratio = 1.59;
        break;
    }    

    float ref_val = rec.front_face ? (1.0f / refraction_ratio) : refraction_ratio;
    direction = refract(unit_direction, normalize(rec.normal), float(ref_val));
    
    scattered = r_in;
    scattered.rgb = r_in.rgb;
    scattered.origin = rec.p;
    scattered.direction = glm::normalize(direction);

    return true;
}

bool CPUMode::DielectricScatter(
    Ray& r_in, hit_record& rec, glm::vec4& attenuation, Ray& scattered, Triangle& tri
) {
    attenuation = glm::vec4(1.0, 1.0, 1.0, 1.0);
    
    float refraction_ratio = rec.front_face ? (1.0f / global_refract_index) : ( global_refract_index);
    glm::vec3 normal = rec.front_face ? glm::normalize(rec.normal) : glm::normalize(-rec.normal);
    glm::vec3 unit_direction = glm::normalize(r_in.direction);


    float cos_theta = glm::min(glm::dot(-unit_direction, rec.normal), 1.0f);
    float sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    glm::vec3 direction;    

    dataText << "Refract with index: " << refraction_ratio << " with direction: " << glm::to_string(unit_direction) << " at normal: " << glm::to_string(rec.normal);

    switch (r_in.rgb)
    {
    case 0:
        refraction_ratio = 1.4;
        break;
    case 1:
        refraction_ratio = 1.6;
        break;
    case 2:
        refraction_ratio = 2.0;
        break;


    default:
        break;
    }

    direction = refract(unit_direction, normalize(rec.normal), float(refraction_ratio));

    dataText << " resulting in direction: " << glm::to_string(direction) << std::endl;

    scattered = r_in;

    scattered.origin = rec.p;
    scattered.direction = direction;    
   
    return true;
}


glm::vec3 CPUMode::random_unit_vector() {
    return glm::normalize(random_in_unit_sphere());
}

glm::vec3 CPUMode::random_in_unit_sphere() {

   
    for (int i = 0; i < 100; i++) {
        glm::vec3 p = randomVec3(-1, 1);
        if ((glm::length(p) * length(p)) >= 1) continue;
        return p;
    }
}

glm::vec3 CPUMode::randomVec3(float min, float max) {
    
    glm::uvec2 i = glm::uvec2(0,0);
    return random_pcg3d(glm::uvec3(i.x, i.y, 0));
}

glm::vec3 CPUMode::random_pcg3d(glm::uvec3 v) {
    v = v * 1664525u + 1013904223u;
    v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
    v ^= v >> 16u;
    v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
    return glm::vec3(v) * (1.0f / float(0xffffffffu));
}

void CPUMode::addCirclePoints(glm::vec3 center, float radius, int numCircles, int numPointsPerCircle, glm::vec3 normal) {
    float radiusStep = radius / float(numCircles);
    for (int j = 0; j < numCircles; j++) {
        float currentRadius = radiusStep * float(j + 1);
        for (int i = 0; i < numPointsPerCircle; i++) {
            float angle = float(i) / float(numPointsPerCircle) * 6.28;
            glm::vec3 point = center + glm::vec3(cos(angle), sin(angle), 0.0) * currentRadius;
            point = rotatePoint(point, center, normal, angle);
            circlePoints.push_back(point);
            numPoints++;
        }
    }
}

glm::vec3 CPUMode::rotatePoint(glm::vec3 point, glm::vec3 center, glm::vec3 axis, float angle) {
    glm::vec3 v = point - center;
    glm::vec3 k = axis / length(axis);
    return center + v * cos(angle) + glm::cross(k, v) * sin(angle) + k * glm::dot(k, v) * (1.0f - cos(angle));
}

bool CPUMode::rayIntersectsPlane(glm::vec3 rayOrigin, glm::vec3 rayDirection, glm::vec3 planeOrigin, glm::vec3 planeNormal, glm::vec3 planeU, glm::vec3 planeV) {
    float d = glm::dot(planeNormal, rayDirection);
    if (glm::abs(d) < 1e-6) return false; // No intersection
    float t = glm::dot(planeNormal, planeOrigin - rayOrigin) / d;
    if (t < 0.0) return false; // Intersection behind ray origin
    glm::vec3 intersectionPoint = rayOrigin + rayDirection * t;
    glm::vec3 localIntersectionPoint = intersectionPoint - planeOrigin;
    float u = glm::dot(localIntersectionPoint, glm::normalize(planeU));
    float v = glm::dot(localIntersectionPoint, glm::normalize(planeV));
    float planeWidth = length(planeU);
    float planeHeight = length(planeV);
    return u >= 0.0 && u <= planeWidth && v >= 0.0 && v <= planeHeight;
}

bool CPUMode::CheckDebugRayHit(std::vector<BeginEnd>& rayArr, Ray& r) {

    float t;
    for (BeginEnd& ray:rayArr) {
        if (intersectionRayCylinder(r.origin, r.direction, ray.p1, ray.p2, .01f, t)) {

            r.color = glm::vec4(0.0f, 200.0f, 200.0f, 1.0f);
            if (ray.numHits == 1) {
                r.color = glm::vec4(200.0f, 0.0f, 200.0f, 1.0f);
            }
            if (ray.numHits == 2) {
                r.color = glm::vec4(200.0f, 200.0f, 0.0f, 1.0f);
            }
            if (ray.numHits == 3) {
                r.color = glm::vec4(0.0f, 0.0f, 200.0f, 1.0f);
            }
            if (ray.numHits == 4) {
                r.color = glm::vec4(0.0f, 200.0f, 200.0f, 1.0f);
            }

            return true;
        }
    }
    return false;
}

std::vector<BeginEnd> CPUMode::ComputeDebugRays(glm::vec3 imagePlanePos) {

    int raysPerPixel = 6;
    bool shadow = false;
    float imageWidth = ProgramParams::windowWidth, imageHeight = ProgramParams::windowHeight;
    int div = 1;
    float valX = 0.0f; 
    float valY = 0.0f; 
    bool inLense = false;
    std::vector<BeginEnd> rayArr;

    for (int l = 0; l < div; l++) {
        valX += 0.2f;
        valY = 0.0f;
        for (int j = 0; j < div; j++) {
            valY += 0.2f;
            for (int k = 0; k < raysPerPixel; k++) {
                int step = k;

                Ray debugRay;
                debugRay.inLense = false;
                camera_get_debug_ray(0.5f, 0.5f, debugRay, k, true);

                glm::vec3 rayDir = debugRay.direction;

                Ray camRay;
                camRay.origin = debugRay.origin;
                camRay.direction = rayDir;
                int depth = 5;

                int hitCount = 0;
               
                BeginEnd ray;
                ray.p1 = camRay.origin;
                ray.p2 = camRay.origin + glm::normalize(camRay.direction) * 200.0f;
                ray.inLense = false;
                ray.numHits = 0;

                rayArr.push_back(ray);
                bool no_reflect = false;
                for (int i = 1; i <= depth; i++) {
                    bool hitPoint = false;


                    ray_color2(camRay, hitPoint, shadow, no_reflect);

                    ray.p1 = camRay.origin;
                    ray.p2 = camRay.origin + glm::normalize(camRay.direction) * 200.0f;
                    ray.inLense = camRay.inLense;

                    rayArr.push_back(ray);

                    if (hitPoint) {

                        rayArr[i - 1].p2 = camRay.origin;
                        hitCount++;

                    }
                    else {
                        rayArr[i - 1].p2 = camRay.origin + glm::normalize(camRay.direction) * 200.0f;

                        break;
                    }
                    rayArr[i].numHits = i;

                }                


            }
        }

    }
    return rayArr;
}


std::ofstream CPUMode::OpenData() {
    std::ofstream datei("datei.txt");
    if (datei.is_open()) {
        return datei;
    }
    else {
        std::cout << "Datei konnte nicht geöffnet werden.";
    }
}

bool CPUMode::drawColoredRays(Ray& r) {

    int raysPerPixel = 6; numPoints;
    bool shadow = false;
    float imageWidth = ProgramParams::windowWidth, imageHeight = ProgramParams::windowHeight;
    int div = 1;
    double valX = 0.0f;
    double valY = 0.0f;
    bool inLense = false;
    bool no_shadow = false, reflect = false;   
    double partX = (imageWidth / 4)/ imageWidth;
    double partY = (imageHeight / 4)/ imageHeight;
    double x = (-0.2f+ (0.4f / raysPerPixel)) + partX;
    double y = (-0.2f+ (0.4f / raysPerPixel)) + partY;
    int step = numPoints/3;
    dataText.is_open();

    for (int col = 0; col < 3; col++) {
        valX = 0;
        
        for (int l = 0; valX <= 1.0; l++) {
            valY = 0.0;
                
            for (int j = 0; valY <= 1.0; j++) {

               
                for (int k = 0; k < raysPerPixel; k++) {
                   
                    Ray debugRay;
                    camera_get_debug_ray(0.5, 0.5, debugRay, k, true);
                    glm::vec3 rayDir = debugRay.direction;

                    Ray camRay;
                    camRay.origin = debugRay.origin;
                    camRay.direction = rayDir;
                    camRay.rgb = col;
                    int depth = 5;


                    int hitCount = 0;
                    BeginEnd rayArr[40];
                    rayArr[0].p1 = camRay.origin;
                    rayArr[0].p2 = camRay.origin + normalize(camRay.direction) * 100.0f;
                    rayArr[0].rgb = col;
                    camRay.id = -k;
                    for (int i = 1; i <= depth; i++) {
                        bool hitPoint = false;

                        ray_color2(camRay, hitPoint, no_shadow, reflect);

                        rayArr[i - 1].hit_normal_orig = camRay.hit_normal_orig;
                        rayArr[i - 1].hit_normal_dir = camRay.hit_normal_dir;

                        rayArr[i].p1 = camRay.origin;
                        rayArr[i].p2 = camRay.origin + normalize(camRay.direction) * 100.0f;
                        rayArr[i].rgb = col;

                        if (hitPoint) {
                            hitCount++;
                            rayArr[i - 1].p2 = camRay.origin;
                        }
                        else {
                            rayArr[i - 1].p2 = camRay.origin + normalize(camRay.direction) * 100.0f;
                            break;
                        }

                      

                        rayArr[i].numHits = i;
                    }

                    float t;
                    

                    // Normale
                    for (int i = 0; i <= hitCount; i++) {
                        if (intersectionRayCylinder(r.origin, r.direction, rayArr[i].hit_normal_orig, (rayArr[i].hit_normal_orig + (normalize(rayArr[i].hit_normal_dir) * 0.1f)), .002f, t)) {
                            if (rayArr[i].hit_normal_orig != glm::vec3(0.0f)) {                               

                                if (rayArr[i].hit_normal_dir.x > 0.0f) {
                                    r.color = glm::vec4(0.0f, 255.0f, 0.0f, 1.0f);                                   
                                }
                                else {
                                    r.color = glm::vec4(255.0f, 0.0f, 0.0f, 1.0f);
                                  
                                }
                               
                                dataText.close();
                                return true;
                            }
                        }
                    }

                    for (int i = 0; i <= hitCount; i++) {
                        if (intersectionRayCylinder(r.origin, r.direction, rayArr[i].p1, rayArr[i].p2, .002f, t)) {

                            if (rayArr[i].rgb == 0) {
                                r.color = glm::vec4(255.0f, 0.0f, 0.0f, 1.0f);
                            }
                            if (rayArr[i].rgb == 1) {
                                r.color = glm::vec4(0.0f, 255.0f, 0.0f, 1.0f);
                            }
                            if (rayArr[i].rgb == 2) {
                                r.color = glm::vec4(0.0f, 0.0f, 255.0f, 1.0f);
                            }

                            dataText.close();
                            return true;
                        }
                    }

                }
                valY += partY;
            }
            valX += partX;
        }
    }
    dataText.close();
    return false;
}



bool CPUMode::intersectionRayCylinder(glm::vec3 rayOrigin, glm::vec3 rayDir, glm::vec3 cylinderStart, glm::vec3 cylinderEnd, float cylinderRadius, float& t) {
    glm::vec3 AB = cylinderEnd - cylinderStart;
    glm::vec3 AO = rayOrigin - cylinderStart;
    glm::vec3 AOxAB = glm::cross(AO, AB);
    glm::vec3 VxAB = glm::cross(rayDir, AB);
    float ab2 = glm::dot(AB, AB);
    float a = glm::dot(VxAB, VxAB);
    float b = 2.0 * glm::dot(VxAB, AOxAB);
    float c = glm::dot(AOxAB, AOxAB) - (cylinderRadius * cylinderRadius * ab2);
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
    glm::vec3 P = rayOrigin + rayDir * t;
    float y = glm::dot(P - cylinderStart, AB) / ab2;
    if (y < 0.0 || y > 1.0) {
        return false;
    }

    // Test end caps
    float tCap1, tCap2;
    bool hitCap1 = intersectRayPlane(rayOrigin, rayDir, cylinderStart, glm::normalize(AB), tCap1);
    bool hitCap2 = intersectRayPlane(rayOrigin, rayDir, cylinderEnd, glm::normalize(AB), tCap2);

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


bool CPUMode::intersectRayPlane(glm::vec3 rayOrigin, glm::vec3 rayDir, glm::vec3 planePoint, glm::vec3 planeNormal, float& t) {
    float denom = glm::dot(planeNormal, rayDir);

    if (glm::abs(denom) > 1e-6) {
        glm::vec3 planeToRayOrigin = planePoint - rayOrigin;
        t = glm::dot(planeToRayOrigin, planeNormal) / denom;

        return (t >= 0.0);
    }

    return false;
}