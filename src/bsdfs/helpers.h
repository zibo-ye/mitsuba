#ifndef HELPERS_H
#define HELPERS_H

#include <cmath>
#include <string>
#include <assert.h>
//#include <random>
//#include <ctime>

//#include <Eigen/Dense>

using namespace std;
//using namespace Eigen;
using namespace mitsuba;

// #define Float float
// #define Vector3 Vector3f

const Float EPS = 1e-6f;

inline Float clamp(Float x, Float a, Float b) {
    if (x < a) return a;
    if (x > b) return b;
    return x;
}

inline int clamp(int x, int a, int b) {
    if (x < a) x = a;
    if (x > b) x = b;
    return x;
}

inline Vector3 clamp(Vector3 x, Float a, Float b) {
    for (int i = 0; i < 3; i++)
        x[i] = clamp(x[i], a, b);
    return x;
}

inline Vector3 max(Vector3 x, Float a) {
    for (int i = 0; i < 3; i++)
        x[i] = std::max(x[i], a);
    return x;
}

inline Vector3 max(Float a, Vector3 x) {
    return max(x, a);
}

inline Vector3 min(Vector3 x, Float a) {
    for (int i = 0; i < 3; i++)
        x[i] = std::min(x[i], a);
    return x;
}

inline Vector3 min(Float a, Vector3 x) {
    return min(x, a);
}

inline Vector3 discardLarger(Vector3 x, Float a) {
    for (int i = 0; i < 3; i++)
        if (x[i] > a){
            x[i] = 0.0f;
        }
    return x;
}

inline Vector3 discardLarger(Float a, Vector3 x) {
    return discardLarger(x, a);
}

inline Float acosSafe(Float cosTheta) {
    return acos((Float) clamp(cosTheta, Float(-1.0 + EPS), Float(1.0 - EPS)));
}

inline Float asinSafe(Float sinTheta) {
    return asin((Float) clamp(sinTheta, Float(-1.0 + EPS), Float(1.0 - EPS)));
}

inline Float G(Float x, Float sigma) {
    if (sigma <= 0.0) {
        //printf("In function G: (sigma = %f) <= 0!\n", sigma);
        return 0.0f;
    }
    return 1.0f / (sqrt(2.0f * M_PI) * sigma) * exp(-0.5f * (x / sigma) * (x / sigma));
}

inline Float GSafe(Float x, Float sigma) {
    if (sigma <= 0.0) {
        return 0.0f;
    }

    Float spreadAngle = sigma / M_PI * 180.0;

    Float gl = 1.0f / (sqrt(2.0f * M_PI) * sigma) * exp(-0.5f * (x / sigma) * (x / sigma));
    Float hf = 1.0f / M_PI;

    if (spreadAngle < 30.0) {
        return gl;
    } else if (spreadAngle > 52.0) {
        return hf;
    } else {
        Float a = (52.0 - spreadAngle) / (52.0 - 30.0);
        return gl * a + hf * (1.0 - a);
    }
}

inline Float laplacian(Float x, Float b) {
    return 1.0 / (2.0 * b) * exp(-fabs(x) / b);
}

inline Float laplacian(Float x, Float mu, Float b) {
    return laplacian(x - mu, b);
}

//Float randFloat() {
//    thread_local std::mt19937 generator(std::random_device{}());
//    std::uniform_real_distribution<Float> distribution(Float(0.0f), Float(1.0f));
//    return distribution(generator);
//}

inline Float randFloat() {
    return rand() / (Float) RAND_MAX;
}

inline Float sampleGaussian(Float mu, Float sigma, Float u1, Float u2) {
    u1 = clamp(u1, EPS, 1.0 - EPS);
    u2 = clamp(u2, EPS, 1.0 - EPS);
    const double two_pi = 2.0 * 3.14159265358979323846;
    double z0, z1;
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    return z0 * sigma + mu;
}

// Float degToRad(Float x) {
//     return x / 180.0 * M_PI;
// }

// Float radToDeg(Float x) {
//     return x / M_PI * 180.0;
// }

inline Float fresnel(Float theta, Float n1, Float n2) {
    if (n1 == n2) return 0.0;

    theta = clamp(theta, -M_PI / 2.0 + EPS, M_PI / 2.0 - EPS);

    Float cosThetaI = cos(theta);
    Float sinThetaI2 = 1.0 - cosThetaI * cosThetaI;
    Float sinThetaT2 = n1 * n1 * sinThetaI2 / (n2 * n2);
    if (sinThetaT2 > 1.0)
        return 1.0;
     
    Float rSFirstTerm = n1 * cosThetaI;
    Float rSSecondTerm = n2 * sqrt(1.0 - sinThetaT2);
    Float rS = (rSFirstTerm - rSSecondTerm) / (rSFirstTerm + rSSecondTerm);
    rS *= rS;
 
    Float rPFirstTerm = n1 * sqrt(1.0 - sinThetaT2);
    Float rPSecondTerm = n2 * cosThetaI;
    Float rP = (rPFirstTerm - rPSecondTerm) / (rPFirstTerm + rPSecondTerm);
    rP *= rP;
 
    return (rS + rP) / 2.0;
}

inline Float fresnelCuticle(Float theta, Float n1, Float n2, Float layers) {
    if (n1 == n2) return 0.0;

    theta = clamp(theta, -M_PI / 2.0 + EPS, M_PI / 2.0 - EPS);

    Float cosThetaI = cos(theta);
    Float sinThetaI2 = 1.0 - cosThetaI * cosThetaI;
    Float sinThetaT2 = n1 * n1 * sinThetaI2 / (n2 * n2);
    if (sinThetaT2 > 1.0)
        return 1.0;

    Float rSFirstTerm = n1 * cosThetaI;
    Float rSSecondTerm = n2 * sqrt(1.0 - sinThetaT2);
    Float rS = (rSFirstTerm - rSSecondTerm) / (rSFirstTerm + rSSecondTerm);
    rS *= rS;

    Float rPFirstTerm = n1 * sqrt(1.0 - sinThetaT2);
    Float rPSecondTerm = n2 * cosThetaI;
    Float rP = (rPFirstTerm - rPSecondTerm) / (rPFirstTerm + rPSecondTerm);
    rP *= rP;

    Float phi_1 = 0.5 * (rS + (1.0 - rS) * (1.0 - rS) * rS / (1.0 - rS * rS)) +
                  0.5 * (rP + (1.0 - rP) * (1.0 - rP) * rP / (1.0 - rP * rP));

    Float phi_m = layers * phi_1 / (1.0 + (layers - 1.0) * phi_1);

    return phi_m;
}

// For non-separable R lobe. All in radians.
inline Float theta_cone_R(Float theta_i, Float alpha, Float phi) {
    Float u = cos(phi / 2.0) * cos(alpha) * cos(theta_i) + sin(alpha) * sin(theta_i);
    Float t = sin(theta_i) - 2.0 * sin(alpha) * u;
    return -asinSafe(t);
}

inline Float beta_cone_R(Float beta, Float phi) {
    return beta * sqrt(2.0) * cos(phi / 2.0);
}

inline void regularizeTheta(Float &theta) {
    if (theta > M_PI / 2.0f)
        theta = M_PI - theta;
    else if (theta < -M_PI / 2.0f)
        theta = -M_PI - theta;
}

inline Float congjugateTheta(Float theta) {
    if (theta > 0.0f)
        theta = M_PI - theta;
    else
        theta = -M_PI - theta;
    return theta;
}

inline void regularizePhi(Float &phi) {
    while (phi > M_PI) phi -= 2.0f * M_PI;
    while (phi < -M_PI) phi += 2.0f * M_PI;
}

inline Vector3 pow3(const Vector3 &v, Float p) {
    return Vector3(pow(v[0], p), pow(v[1], p), pow(v[2], p));
}

inline Vector3 exp3(const Vector3 &v) {
    return Vector3(exp(v[0]), exp(v[1]), exp(v[2]));
}

inline Vector3 mul3(const Vector3 &v1, const Vector3 &v2) {
    return Vector3(v1[0] * v2[0], v1[1] * v2[1], v1[2] * v2[2]);
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}

char* getCmdOption(char **begin, char **end, const std::string & option) {
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
        return *itr;
    }
    return 0;
}

template<typename T>
T getCmdValue(int argc, char *argv[], const std::string & option, T defaultValue) {
    char *s = getCmdOption(argv, argv + argc, option);
    return s? T(atof(s)) : defaultValue;
}

#endif
