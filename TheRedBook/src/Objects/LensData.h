#pragma once
#include "../ProgramModes/ProgramParams.h"

struct LensData {

    glm::vec3 lensOrigin;

    glm::vec3 lensOriginLeft;
    glm::vec3 lensOriginRight;

    double lensRadiusLeft, lensRadiusRight;
    double d;
    double n;



    LensData(glm::vec3 lensOrigin_in = glm::vec3(0.0), float n_in = 1.59, glm::vec3 lensOriginLeft_in = glm::vec3(0.9 * 2, 0.0, 0.0), glm::vec3 lensOriginRight_in = glm::vec3(0.9 * 2, 0.0, 0.0), double lensRadiusLeft_in = 2.0, double lensRadiusRight_in = 2.0) {
        lensOrigin = lensOrigin_in;

        lensOriginLeft = lensOrigin_in - lensOriginLeft_in;
        lensOriginRight = lensOrigin_in + lensOriginRight_in;

        lensRadiusLeft = lensRadiusLeft_in;
        lensRadiusRight = lensRadiusRight_in;

        n = n_in;
        std::cout << (glm::length(lensRadiusLeft - lensRadiusRight));
        d = glm::abs(lensRadiusLeft) + glm::abs(lensRadiusRight) - glm::length(lensOriginLeft - lensOriginRight);


        std::cout << "d in lens: " << d;
    }

    double calcOneByF() {
        return (n - 1.0) * ((1.0 / lensRadiusLeft) - (1.0 / (-lensRadiusRight)) + ((n - 1.0) * d) / (n * lensRadiusLeft * (-lensRadiusRight)));
    }

    double calcF() {
        return 1.0 / ((n - 1.0) * ((1.0 / lensRadiusLeft) - (1.0 / (-lensRadiusRight)) + ((n - 1.0) * d) / (n * lensRadiusLeft * (-lensRadiusRight))));
    }

    double calcD() {
        double b = ProgramParams::b, g = ProgramParams::g, r1 = lensRadiusLeft, r2 = -lensRadiusRight;
                
        return (n*(b*g*(n-1)*(r1 - r2)+(b*r1*r2)+(g*r1*r2)))/(b*g*pow((n-1),2));
    }

    void dTolensRad() {
        float dOld = (glm::abs(lensRadiusLeft) + glm::abs(lensRadiusRight) - glm::length(lensOriginLeft - lensOriginRight));
        float diff = d - dOld;

        lensRadiusLeft  += (diff/2);
        lensRadiusRight += (diff/2);
    }

};