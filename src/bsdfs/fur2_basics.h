#ifndef FUR_BASICS_H
#define FUR_BASICS_H

//#include <mitsuba/core/platform.h>
//#include <mitsuba/core/fwd.h>

#include <mitsuba/render/sampler.h>

#include <iostream>
#include <cstdio>
#include <assert.h>
#include "helpers.h"

#define SCATTERED_ATTENUATION 1.0
//#define FINE_PRECOMP

#define MAX_SCATTERING_INNER 16.0
#define MAX_H 1.0
#define MAX_G 1.0
#define MAX_THETA degToRad(90.0)

#define NUM_SCATTERING_INNER 24
#define NUM_H 16
#define NUM_G 16
#define NUM_BINS 720
#define NUM_THETA 16

#define PIXEL_COVERAGE_SCALE 1.0

using namespace std;
using namespace mitsuba;
MTS_NAMESPACE_BEGIN

float scattered[NUM_SCATTERING_INNER][NUM_THETA][NUM_G][NUM_BINS];
float scatteredDist[NUM_SCATTERING_INNER][NUM_THETA][NUM_G][NUM_BINS];
float scatteredM[NUM_SCATTERING_INNER][NUM_THETA][NUM_G][NUM_BINS];


// For dual scattering
Float d_f, d_b;

const Float MAX_ABSORPTION_OUTER = 2.4f;
//const int NUM_SIGMA_SAMPLES = 128;
#if defined(FINE_PRECOMP)
#define NUM_BINS_DS 512
const int NUM_SIGMA_SAMPLES = 63; // mod 3 == 0
const int NUM_THETA_SAMPLES = 512;
const int NUM_PHI_SAMPLES = 512;
int NUM_K_SAMPLES = 32;
#else
#define NUM_BINS_DS 64
const int NUM_SIGMA_SAMPLES = 63; // mod 3 == 0
const int NUM_THETA_SAMPLES = 16;
const int NUM_PHI_SAMPLES = 16;
int NUM_K_SAMPLES = 1;
#endif

float alpha_b[NUM_BINS_DS][NUM_SIGMA_SAMPLES], alpha_f[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float beta_b_sqr[NUM_BINS_DS][NUM_SIGMA_SAMPLES], beta_f_sqr[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float a_b[NUM_BINS_DS][NUM_SIGMA_SAMPLES], a_f[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float a_bu[NUM_BINS_DS][NUM_SIGMA_SAMPLES], a_fu[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float a_bs[NUM_BINS_DS][NUM_SIGMA_SAMPLES], a_fs[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float A_b[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float Delta_b[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float sigma_b_sqr[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float N_R_G[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float N_TT_G[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float N_TRT_G[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float N_TTs_G[NUM_BINS_DS][NUM_SIGMA_SAMPLES];
float N_TRTs_G[NUM_BINS_DS][NUM_SIGMA_SAMPLES];

float interpDS(float matrix[][NUM_SIGMA_SAMPLES], float theta_i, float sigma) {
    Float thetaBin = (theta_i + M_PI / 2.0f) / M_PI * NUM_BINS_DS - 0.5f;
    Float sigmaBin = sigma / MAX_ABSORPTION_OUTER * NUM_SIGMA_SAMPLES - 0.5f;

    int thetaBin1 = clamp((int) thetaBin, 0, NUM_BINS_DS - 1);
    int thetaBin2 = clamp((int) thetaBin + 1, 0, NUM_BINS_DS - 1);
    Float thetaBinRes = thetaBin - floor(thetaBin);

    int sigmaBin1 = clamp((int) sigmaBin, 0, NUM_SIGMA_SAMPLES - 1);
    int sigmaBin2 = clamp((int) sigmaBin + 1, 0, NUM_SIGMA_SAMPLES - 1);
    Float sigmaBinRes = sigmaBin - floor(sigmaBin);

    Float m11 = matrix[thetaBin1][sigmaBin1];
    Float m12 = matrix[thetaBin1][sigmaBin2];
    Float m21 = matrix[thetaBin2][sigmaBin1];
    Float m22 = matrix[thetaBin2][sigmaBin2];

    Float p1 = m11 * (1.0 - sigmaBinRes) + m12 * sigmaBinRes;
    Float p2 = m21 * (1.0 - sigmaBinRes) + m22 * sigmaBinRes;

    return p1 * (1.0 - thetaBinRes) + p2 * thetaBinRes;
}

Vector3 interpDS(float matrix[][NUM_SIGMA_SAMPLES], float theta_i, const Vector3 &sigma) {
    return Vector3(interpDS(matrix, theta_i, sigma[0]),
            interpDS(matrix, theta_i, sigma[1]),
            interpDS(matrix, theta_i, sigma[2]));
}

Vector3 operator*(const Vector3& lhs, const Vector3& rhs) {
    return Vector3(lhs[0] * rhs[0], lhs[1] * rhs[1], lhs[2] * rhs[2]);
}

void operator*=(Vector3 &lhs, const Vector3& rhs) {
    lhs[0] *= rhs[0];
    lhs[1] *= rhs[1];
    lhs[2] *= rhs[2];
}

Vector3 operator/(const Vector3& lhs, const Vector3& rhs) {
    return Vector3(lhs[0] / rhs[0], lhs[1] / rhs[1], lhs[2] / rhs[2]);
}

void operator/=(Vector3 &lhs, const Vector3& rhs) {
    lhs[0] /= rhs[0];
    lhs[1] /= rhs[1];
    lhs[2] /= rhs[2];
}

Float maxComponent(const Vector3& v) {
    return std::max(std::max(v[0], v[1]), v[2]);
}

// Precomputed data for medulla

void initialize(const char *medullaFilename) {
    //    size_t bytesRead;

    //    // Read precomputed medulla profiles (azimuthal)
    //    FILE *fMedullaN = fopen((std::string(medullaFilename) + "Azimuthal.bin").c_str(), "rb");
    //    // Unscattered is not used, so here we read and overwrite it.
    //    bytesRead = fread(scattered, sizeof(float), NUM_SCATTERING_INNER * NUM_H * NUM_G, fMedullaN);
    //    assert(bytesRead != 0);
    //    bytesRead = fread(scattered, sizeof(float), NUM_SCATTERING_INNER * NUM_H * NUM_G * NUM_BINS, fMedullaN);
    //    assert(bytesRead != 0);
    //    bytesRead = fread(scatteredDist, sizeof(float), NUM_SCATTERING_INNER * NUM_H * NUM_G * NUM_BINS, fMedullaN);
    //    assert(bytesRead != 0);
    //    fclose(fMedullaN);

    //    // Read precomputed medulla profiles (longitudinal)
    //    FILE *fMedullaM = fopen((std::string(medullaFilename) + "Longitudinal.bin").c_str(), "rb");
    //    bytesRead = fread(scatteredM, sizeof(float), NUM_SCATTERING_INNER * NUM_THETA * NUM_G * NUM_BINS, fMedullaM);
    //    assert(bytesRead != 0);
    //    fclose(fMedullaM);

    size_t bytesRead;

    // Read compressed scattering profiles.
    FILE *fp = fopen((std::string(medullaFilename) + "Compressed.bin").c_str(), "rb");

    if (!fp) {
        printf("Can't read precomputed scattering profiles!\n");
        return;
    }

    // Azimuthal
    float n1_lambda[16];
    float n1_u0[23][16];
    float n1_u1[16][16];
    float n1_u2[16][16];
    float n1_u3[720][16];

    bytesRead = fread(n1_lambda, sizeof (float), 16, fp);
    bytesRead = fread(n1_u0, sizeof (float), 23 * 16, fp);
    bytesRead = fread(n1_u1, sizeof (float), 16 * 16, fp);
    bytesRead = fread(n1_u2, sizeof (float), 16 * 16, fp);
    bytesRead = fread(n1_u3, sizeof (float), 720 * 16, fp);

    // Azimuthal Dist
    float n2_lambda[16];
    float n2_u0[23][16];
    float n2_u1[16][16];
    float n2_u2[16][16];
    float n2_u3[720][16];

    bytesRead = fread(n2_lambda, sizeof (float), 16, fp);
    bytesRead = fread(n2_u0, sizeof (float), 23 * 16, fp);
    bytesRead = fread(n2_u1, sizeof (float), 16 * 16, fp);
    bytesRead = fread(n2_u2, sizeof (float), 16 * 16, fp);
    bytesRead = fread(n2_u3, sizeof (float), 720 * 16, fp);

    // Longitudinal (backward)
    float m1_lambda[16];
    float m1_u0[16][16];
    float m1_u1[16][16];
    float m1_u2[16][16];
    float m1_u3[360][16];

    bytesRead = fread(m1_lambda, sizeof (float), 16, fp);
    bytesRead = fread(m1_u0, sizeof (float), 16 * 16, fp);
    bytesRead = fread(m1_u1, sizeof (float), 16 * 16, fp);
    bytesRead = fread(m1_u2, sizeof (float), 16 * 16, fp);
    bytesRead = fread(m1_u3, sizeof (float), 360 * 16, fp);

    // Longitudinal (backward)
    float m2_lambda[16];
    float m2_u0[16][16];
    float m2_u1[16][16];
    float m2_u2[16][16];
    float m2_u3[359][16];

    bytesRead = fread(m2_lambda, sizeof (float), 16, fp);
    bytesRead = fread(m2_u0, sizeof (float), 16 * 16, fp);
    bytesRead = fread(m2_u1, sizeof (float), 16 * 16, fp);
    bytesRead = fread(m2_u2, sizeof (float), 16 * 16, fp);
    bytesRead = fread(m2_u3, sizeof (float), 359 * 16, fp);

    fclose(fp);

    // Reconstruct the actual scattering profiles.
    memset(scattered, 0, sizeof (float) * NUM_SCATTERING_INNER * NUM_H * NUM_G * NUM_BINS);
    memset(scatteredDist, 0, sizeof (float) * NUM_SCATTERING_INNER * NUM_H * NUM_G * NUM_BINS);
    memset(scatteredM, 0, sizeof (float) * 16 * NUM_THETA * NUM_G * NUM_BINS);

    for (int i = 1; i < NUM_SCATTERING_INNER; i++)
        for (int j = 0; j < NUM_H; j++)
            for (int k = 0; k < NUM_G; k++)
                for (int l = 0; l < NUM_BINS; l++)
                    for (int t = 0; t < 16; t++) {
                        scattered[i][j][k][l] += n1_lambda[t] * n1_u0[i - 1][t] * n1_u1[j][t] * n1_u2[k][t] * n1_u3[l][t];
                        scatteredDist[i][j][k][l] += n2_lambda[t] * n2_u0[i - 1][t] * n2_u1[j][t] * n2_u2[k][t] * n2_u3[l][t];
                    }

    for (int i = 0; i < 16; i++)
        for (int j = 0; j < NUM_THETA; j++)
            for (int k = 0; k < NUM_G; k++) {
                for (int l = 0; l < NUM_BINS / 2; l++)
                    for (int t = 0; t < 16; t++)
                        scatteredM[i][j][k][l] += m1_lambda[t] * m1_u0[i][t] * m1_u1[j][t] * m1_u2[k][t] * m1_u3[l][t];
                for (int l = 0; l < NUM_BINS / 2 - 1; l++)
                    for (int t = 0; t < 16; t++)
                        scatteredM[i][j][k][l + NUM_BINS / 2 + 1] += m2_lambda[t] * m2_u0[i][t] * m2_u1[j][t] * m2_u2[k][t] * m2_u3[l][t];
            }
}

float interpScatteredN(float binSigma, float binH, float binG, float binPhi, int stage = 0) {
    switch (stage) {
        case 0:
        {
            if (binSigma >= NUM_SCATTERING_INNER - 1)
                return interpScatteredN(NUM_SCATTERING_INNER - 1, binH, binG, binPhi, 1);
            if (binSigma < 1)
                return interpScatteredN(1, binH, binG, binPhi, 1);
            int binSigmaIntA = (int) binSigma;
            int binSigmaIntB = (int) binSigma + 1;
            Float w2 = binSigma - binSigmaIntA;
            Float w1 = 1.0 - w2;
            return interpScatteredN(binSigmaIntA, binH, binG, binPhi, 1) * w1 +
                    interpScatteredN(binSigmaIntB, binH, binG, binPhi, 1) * w2;
            break;
        }
        case 1:
        {
            if (binH >= NUM_H - 1)
                return interpScatteredN(binSigma, NUM_H - 1, binG, binPhi, 2);
            int binHIntA = (int) binH;
            int binHIntB = (int) binH + 1;
            Float w2 = binH - binHIntA;
            Float w1 = 1.0 - w2;
            return interpScatteredN(binSigma, binHIntA, binG, binPhi, 2) * w1 +
                    interpScatteredN(binSigma, binHIntB, binG, binPhi, 2) * w2;
            break;
        }
        case 2:
        {
            if (binG >= NUM_G - 1)
                return interpScatteredN(binSigma, binH, NUM_G - 1, binPhi, 3);
            int binGIntA = (int) binG;
            int binGIntB = (int) binG + 1;
            Float w2 = binG - binGIntA;
            Float w1 = 1.0 - w2;
            return interpScatteredN(binSigma, binH, binGIntA, binPhi, 3) * w1 +
                    interpScatteredN(binSigma, binH, binGIntB, binPhi, 3) * w2;
            break;
        }
        case 3:
        {
            if (binPhi >= NUM_BINS - 1)
                return scattered[(int) binSigma][(int) binH][(int) binG][(int) NUM_BINS - 1];
            int binPhiIntA = (int) binPhi;
            int binPhiIntB = (int) binPhi + 1;
            Float w2 = binPhi - binPhiIntA;
            Float w1 = 1.0 - w2;
            return scattered[(int) binSigma][(int) binH][(int) binG][(int) binPhiIntA] * w1 +
                    scattered[(int) binSigma][(int) binH][(int) binG][(int) binPhiIntB] * w2;
            break;
        }
    }
}

float interpScatteredDist(float binSigma, float binH, float binG, float binPhi, int stage = 0) {
    switch (stage) {
        case 0:
        {
            if (binSigma >= NUM_SCATTERING_INNER - 1)
                return interpScatteredDist(NUM_SCATTERING_INNER - 1, binH, binG, binPhi, 1);
            if (binSigma < 1)
                return interpScatteredDist(1, binH, binG, binPhi, 1);
            int binSigmaIntA = (int) binSigma;
            int binSigmaIntB = (int) binSigma + 1;
            Float w2 = binSigma - binSigmaIntA;
            Float w1 = 1.0 - w2;
            return interpScatteredDist(binSigmaIntA, binH, binG, binPhi, 1) * w1 +
                    interpScatteredDist(binSigmaIntB, binH, binG, binPhi, 1) * w2;
            break;
        }
        case 1:
        {
            if (binH >= NUM_H - 1)
                return interpScatteredDist(binSigma, NUM_H - 1, binG, binPhi, 2);
            int binHIntA = (int) binH;
            int binHIntB = (int) binH + 1;
            Float w2 = binH - binHIntA;
            Float w1 = 1.0 - w2;
            return interpScatteredDist(binSigma, binHIntA, binG, binPhi, 2) * w1 +
                    interpScatteredDist(binSigma, binHIntB, binG, binPhi, 2) * w2;
            break;
        }
        case 2:
        {
            if (binG >= NUM_G - 1)
                return interpScatteredDist(binSigma, binH, NUM_G - 1, binPhi, 3);
            int binGIntA = (int) binG;
            int binGIntB = (int) binG + 1;
            Float w2 = binG - binGIntA;
            Float w1 = 1.0 - w2;
            return interpScatteredDist(binSigma, binH, binGIntA, binPhi, 3) * w1 +
                    interpScatteredDist(binSigma, binH, binGIntB, binPhi, 3) * w2;
            break;
        }
        case 3:
        {
            if (binPhi >= NUM_BINS - 1)
                return scatteredDist[(int) binSigma][(int) binH][(int) binG][(int) NUM_BINS - 1];
            int binPhiIntA = (int) binPhi;
            int binPhiIntB = (int) binPhi + 1;
            Float w2 = binPhi - binPhiIntA;
            Float w1 = 1.0 - w2;
            return scatteredDist[(int) binSigma][(int) binH][(int) binG][(int) binPhiIntA] * w1 +
                    scatteredDist[(int) binSigma][(int) binH][(int) binG][(int) binPhiIntB] * w2;
            break;
        }
    }
}

float interpScatteredM(float binSigma, float binThetaI, float binG, float binThetaR, int stage = 0) {
    switch (stage) {
        case 0:
        {
            if (binSigma >= 16 - 1)
                return interpScatteredM(16 - 1, binThetaI, binG, binThetaR, 1);
            if (binSigma < 1)
                return interpScatteredM(1, binThetaI, binG, binThetaR, 1);
            int binSigmaIntA = (int) binSigma;
            int binSigmaIntB = (int) binSigma + 1;
            Float w2 = binSigma - binSigmaIntA;
            Float w1 = 1.0 - w2;
            return interpScatteredM(binSigmaIntA, binThetaI, binG, binThetaR, 1) * w1 +
                    interpScatteredM(binSigmaIntB, binThetaI, binG, binThetaR, 1) * w2;
            break;
        }
        case 1:
        {
            if (binThetaI >= NUM_THETA - 1)
                return interpScatteredM(binSigma, NUM_THETA - 1, binG, binThetaR, 2);
            int binThetaIIntA = (int) binThetaI;
            int binThetaIIntB = (int) binThetaI + 1;
            Float w2 = binThetaI - binThetaIIntA;
            Float w1 = 1.0 - w2;
            return interpScatteredM(binSigma, binThetaIIntA, binG, binThetaR, 2) * w1 +
                    interpScatteredM(binSigma, binThetaIIntB, binG, binThetaR, 2) * w2;
            break;
        }
        case 2:
        {
            if (binG >= NUM_G - 1)
                return interpScatteredM(binSigma, binThetaI, NUM_G - 1, binThetaR, 3);
            int binGIntA = (int) binG;
            int binGIntB = (int) binG + 1;
            Float w2 = binG - binGIntA;
            Float w1 = 1.0 - w2;
            return interpScatteredM(binSigma, binThetaI, binGIntA, binThetaR, 3) * w1 +
                    interpScatteredM(binSigma, binThetaI, binGIntB, binThetaR, 3) * w2;
            break;
        }
        case 3:
        {
            if (binThetaR >= NUM_BINS - 1)
                return scatteredM[(int) binSigma][(int) binThetaI][(int) binG][(int) NUM_BINS - 1];
            int binThetaRIntA = (int) binThetaR;
            int binThetaRIntB = (int) binThetaR + 1;
            Float w2 = binThetaR - binThetaRIntA;
            Float w1 = 1.0 - w2;
            return scatteredM[(int) binSigma][(int) binThetaI][(int) binG][(int) binThetaRIntA] * w1 +
                    scatteredM[(int) binSigma][(int) binThetaI][(int) binG][(int) binThetaRIntB] * w2;
            break;
        }
    }
}

class FurSolver {
public:
    Float eta;
    Vector3 sigma_a;
    Float alpha;
    Float beta_m, beta_n;
    Float kappa;
    Float sigma_a_m, sigma_s_m;
    Float g;
    Float layers;
    Sampler *sampler;
public:

    FurSolver(Float eta, Vector3 sigma_a, Float alpha, Float beta_m, Float beta_n,
            Float kappa, Float sigma_a_m, Float sigma_s_m, Float g, Float layers, Sampler *sampler_ = NULL) :
    eta(eta), sigma_a(sigma_a), alpha(alpha), beta_m(beta_m), beta_n(beta_n),
    kappa(kappa), sigma_a_m(sigma_a_m), sigma_s_m(sigma_s_m), g(g), layers(layers) {
        if (sampler_)
            sampler = sampler_;
    }

    Float getScatteredN(Float sigma_s, Float h, Float g, Float phi) const {
        if (h < 0.0) {
            h = -h;
            phi = -phi;
        }

        Float binSigma = sigma_s / MAX_SCATTERING_INNER * NUM_SCATTERING_INNER;
        Float binH = h / MAX_H * NUM_H;
        Float binG = g / MAX_G * NUM_G;

        Float binPhi;
        if (phi < 0)
            binPhi = (-phi) / (2.0 * M_PI) * NUM_BINS;
        else
            binPhi = (2.0 * M_PI - phi) / (2.0 * M_PI) * NUM_BINS;

        return interpScatteredN(binSigma, binH, binG, binPhi);
    }

    Float getScatteredDist(Float sigma_s, Float h, Float g, Float phi) const {
        if (h < 0.0) {
            h = -h;
            phi = -phi;
        }

        Float binSigma = sigma_s / MAX_SCATTERING_INNER * NUM_SCATTERING_INNER;
        Float binH = h / MAX_H * NUM_H;
        Float binG = g / MAX_G * NUM_G;
        // if (phi < 0) phi += 2.0 * M_PI;
        Float binPhi;
        if (phi < 0)
            binPhi = (-phi) / (2.0 * M_PI) * NUM_BINS;
        else
            binPhi = (2.0 * M_PI - phi) / (2.0 * M_PI) * NUM_BINS;

        return interpScatteredDist(binSigma, binH, binG, binPhi);
    }

    //    inline Float getScatteredN(Float sigma_s, Float h, Float g, Float phi) const {
    //        if (h < 0.0) {
    //            h = -h;
    //            phi = -phi;
    //        }

    //        Float binSigma = sigma_s / MAX_SCATTERING_INNER * NUM_SCATTERING_INNER;
    //        Float binH = h / MAX_H * NUM_H;
    //        Float binG = g / MAX_G * NUM_G;

    //        Float binPhi;
    //        if (phi < 0)
    //            binPhi = (-phi) / (2.0 * M_PI) * NUM_BINS;
    //        else
    //            binPhi = (2.0 * M_PI - phi) / (2.0 * M_PI) * NUM_BINS;

    //        return scattered[(int) clamp(binSigma + 0.5, 0.0, NUM_SCATTERING_INNER - 1.0)]
    //                        [(int) clamp(binH + 0.5, 0.0, NUM_H - 1.0)]
    //                        [(int) clamp(binG + 0.5, 0.0, NUM_G - 1.0)]
    //                        [(int) clamp(binPhi + 0.5, 0.0, NUM_BINS - 1.0)];
    //    }

    //    inline Float getScatteredDist(Float sigma_s, Float h, Float g, Float phi) const {
    //        Float t = kappa * kappa - h * h;
    //        if (t < 0.0)
    //            return 0.0;
    //        return 2.0 * std::sqrt(t);
    //    }

    // For non-separable R lobe. All in radians.

    Float theta_cone_R(Float theta_i, Float alpha, Float phi) const {
        Float u = cos(phi / 2.0) * cos(alpha) * cos(theta_i) + sin(alpha) * sin(theta_i);
        Float t = sin(theta_i) - 2.0 * sin(alpha) * u;
        return -asinSafe(t);
    }

    Float beta_m_cone_R(Float beta_m, Float phi) const {
        return std::max(beta_m * std::sqrt(2.0) * cos(phi / 2.0), degToRad(0.02));
        //        return beta_m;
    }

    Float Mu(int p, Float theta_i, Float theta_r, Float phi) const {
        Float alpha_p, beta_m_p;

        switch (p) {
            case 0: // R
                alpha_p = alpha;
                beta_m_p = beta_m;
                return G(theta_r - theta_cone_R(theta_i, alpha_p, phi), beta_m_cone_R(beta_m_p, phi));
                break;
            case 1: // TT
                alpha_p = -alpha / 2.0;
                beta_m_p = beta_m; // / 2.0;
                break;
            case 2: // TRT
                alpha_p = -alpha / 2.0 * 3.0;
                beta_m_p = beta_m / 2.0 * 3.0;
                break;
        }

        return G(theta_r - (-theta_i + 2.0 * alpha_p), beta_m_p);
    }

    Float Ms(Float theta_i, Float theta_r, Float phi) const {
        if (theta_i < 0.0) {
            theta_i = -theta_i;
            theta_r = -theta_r;
        }

        Float binSigma = sigma_s_m * kappa / MAX_SCATTERING_INNER * NUM_SCATTERING_INNER;
        Float binThetaI = theta_i / MAX_THETA * NUM_THETA;
        Float binG = g / MAX_G * NUM_G;
        Float binThetaRBack = (theta_r + M_PI / 2.0) / M_PI * (NUM_BINS / 2 - 1);
        Float binThetaRFront = NUM_BINS - 1 - binThetaRBack;
        Float wBack = (M_PI - std::fabs(phi)) / M_PI;
        Float wFront = 1.0 - wBack;

        return interpScatteredM(binSigma, binThetaI, binG, binThetaRBack) * wBack +
                interpScatteredM(binSigma, binThetaI, binG, binThetaRFront) * wFront;
        //        return scatteredM[(int) clamp(binSigma + 0.5, 0.0, NUM_SCATTERING_INNER - 1.0)]
        //                         [(int) clamp(binThetaI + 0.5, 0.0, NUM_THETA - 1.0)]
        //                         [(int) clamp(binG + 0.5, 0.0, NUM_G - 1.0)]
        //                         [(int) clamp(binThetaRBack + 0.5, 0.0, NUM_BINS - 1.0)] * wBack +
        //               scatteredM[(int) clamp(binSigma + 0.5, 0.0, NUM_SCATTERING_INNER - 1.0)]
        //                         [(int) clamp(binThetaI + 0.5, 0.0, NUM_THETA - 1.0)]
        //                         [(int) clamp(binG + 0.5, 0.0, NUM_G - 1.0)]
        //                         [(int) clamp(binThetaRFront + 0.5, 0.0, NUM_BINS - 1.0)] * wFront;
    }

    // For unscattered lobes.

    Vector3 Au(int p, Float h, Float theta_i, Float eta_prime) const {
        Float gamma_i = asinSafe(h);
        Float gamma_t = asinSafe(h / eta_prime);

        Float f = fresnelCuticle(gamma_i, 1.0, eta_prime, layers);

        if (p == 0) {
            return Vector3(f, f, f);
        } else {
            Float theta_t = asinSafe(1.0 / eta * sin(theta_i));
            Vector3 sigma_a_prime = sigma_a / cos(theta_t);
            Float sigma_a_m_prime = sigma_a_m / cos(theta_t);
            Float sigma_s_m_prime = sigma_s_m / cos(theta_t);
            Float sigma_e_m_prime = sigma_a_m_prime + sigma_s_m_prime;

            // Travelling distance in the cortex and medulla, respectively.
            Float d1 = cos(gamma_t), d2 = 0.0;
            Float h_prime = h / eta_prime;
            if (std::fabs(h_prime) < std::fabs(kappa))
                d2 = std::sqrt(kappa * kappa - h_prime * h_prime);
            d1 -= d2;

            return Float(std::pow(1.0 - f, 2.0)) *
                    Float(std::pow(f, p - 1.0)) *
                    pow3(exp3(-2.0 * d1 * sigma_a_prime) *
                    std::exp(-2.0 * d2 * sigma_e_m_prime), p);
        }
    }

    Float Du(Float phi, Float sigma) const {
        Float retVal = 0.0;
        for (int i = -1; i <= 1; i++)
            retVal += G(phi + i * 2.0 * M_PI, sigma);
        return retVal;
    }

    Vector3 Nu(int p, Float h, Float phi, Float theta_i, Float eta_prime) const {
        Float gamma_i = asinSafe(h);
        Float gamma_t = asinSafe(h / eta_prime);

        Vector3 attenuation = Au(p, h, theta_i, eta_prime);
        Float Phi = 2.0 * p * gamma_t - 2.0 * gamma_i + p * M_PI;

        Float beta_n_p_sqr = beta_n * beta_n * (p + 1);
        Float beta_n_p = std::sqrt(beta_n_p_sqr);
        return attenuation * Du(Phi - phi, beta_n_p);
    }

    //    Vector3 Nu(int p, Float phi, Float theta_i) const {
    //        Vector3 ret_val(0.0, 0.0, 0.0);
    //        const int SEGMENTS = 100;

    //        Float dh = 2.0 / SEGMENTS;

    //        Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);

    //        for (int i = 0; i < SEGMENTS; i++) {
    //            Float h;
    //            h = (i + 0.5) / SEGMENTS * 2.0 - 1.0;

    //            ret_val += Nu(p, h, phi, theta_i, eta_prime) * dh;
    //        }

    //        ret_val *= 0.5;

    //        return ret_val;
    //    }

    Vector3 Nu(int p, Float h, Float phi, Float theta_i) const {
        Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);
        return Nu(p, h, phi, theta_i, eta_prime);
    }

    //    // For scattered lobes. (This doesn't consider absorption from medulla, but we don't need it anyway.)
    //    Vector3 As(int p, Float h, Float theta_i, Float eta_prime) const {
    //        assert(p > 0);

    //        Float gamma_i = asinSafe(h);
    //        Float gamma_t = asinSafe(h / eta_prime);

    //        Float f = fresnelCuticle(gamma_i, 1.0, eta_prime, layers);

    //        Float theta_t = asinSafe(1.0 / eta * sin(theta_i));
    //        Vector3 sigma_a_prime = sigma_a / cos(theta_t);
    //        Float sigma_a_m_prime = sigma_a_m / cos(theta_t);
    //        Float sigma_s_m_prime = sigma_s_m / cos(theta_t);
    //        Float sigma_e_m_prime = sigma_a_m_prime + sigma_s_m_prime;

    //        Float d1 = cos(gamma_t), d2 = 0.0;
    //        Float h_prime = h / eta_prime;
    //        if (std::fabs(h_prime) >= std::fabs(kappa))
    //            return Vector3(0.0, 0.0, 0.0);
    //        d2 = std::sqrt(kappa * kappa - h_prime * h_prime);
    //        d1 -= d2;

    //        Vector3 retVal = (1.0 - f) * std::pow(f, p - 1.0) *
    //                mul3(exp3(-(2 * p - 1) * d1 * sigma_a_prime),
    //                     exp3(-(1.0 - kappa) * sigma_a_prime)) *
    //                std::exp(-2 * d2 * (p - 1) * sigma_e_m_prime);

    //        return retVal;
    //    }

    // For scattered lobes.

    Vector3 As(int p, Float h, Float phi, Float theta_i, Float eta_prime) const {
        assert(p > 0);

        Float gamma_i = asinSafe(h);
        Float gamma_t = asinSafe(h / eta_prime);

        Float f = fresnelCuticle(gamma_i, 1.0, eta_prime, layers);

        Float theta_t = asinSafe(1.0 / eta * sin(theta_i));
        Vector3 sigma_a_prime = sigma_a / cos(theta_t);
        Float sigma_a_m_prime = sigma_a_m / cos(theta_t);
        Float sigma_s_m_prime = sigma_s_m / cos(theta_t);
        Float sigma_e_m_prime = sigma_a_m_prime + sigma_s_m_prime;

        Float d1 = cos(gamma_t), d2 = 0.0;
        Float h_prime = h / eta_prime;
        if (std::fabs(h_prime) >= std::fabs(kappa))
            return Vector3(0.0, 0.0, 0.0);
        d2 = std::sqrt(kappa * kappa - h_prime * h_prime);
        d1 -= d2;

        Vector3 retVal = (1.0 - f) * std::pow(f, p - 1.0) *
                mul3(exp3(-(2 * p - 1) * d1 * sigma_a_prime),
                exp3(-(1.0 - kappa) * sigma_a_prime)) *
                std::exp(-2 * d2 * (p - 1) * sigma_e_m_prime);

        Float phi_query = phi - (gamma_t - gamma_i) - (M_PI + 2.0 * gamma_t) * (p - 1);
        regularizePhi(phi_query);

        Float avgDist = kappa * getScatteredDist(sigma_s_m_prime * kappa, h_prime / kappa, g, phi_query);
        		retVal = retVal * std::exp(-avgDist * sigma_a_m_prime)*(Vector3(1.) -
			pow3(exp3(-2.0 * d1 * sigma_a_prime) *
			std::exp(-2.0 * d2 * sigma_e_m_prime), p));

        return retVal;
    }

    Float Ds(int p, Float h, Float phi, Float theta_i, Float eta_prime) const {
        assert(p > 0);

        Float gamma_i = asinSafe(h);
        Float gamma_t = asinSafe(h / eta_prime);

        Float theta_t = asinSafe(1.0 / eta * sin(theta_i));
        Float sigma_s_m_prime = sigma_s_m / cos(theta_t);

        Float h_prime = h / eta_prime;

        Float phi_query = phi - (gamma_t - gamma_i) - (M_PI + 2.0 * gamma_t) * (p - 1);
        regularizePhi(phi_query);
        return getScatteredN(sigma_s_m_prime * kappa, h_prime / kappa, g, phi_query);
    }

    Vector3 Ns(int p, Float h, Float phi, Float theta_i, Float eta_prime) const {
        if (p == 0)
            return Vector3(0.0, 0.0, 0.0);

        Float h_prime = h / eta_prime;
        if (std::fabs(h_prime) >= std::fabs(kappa))
            return Vector3(0.0, 0.0, 0.0);

        Vector3 attenuation = As(p, h, phi, theta_i, eta_prime);
        return attenuation * Ds(p, h, phi, theta_i, eta_prime);
    }

    //    Vector3 Ns(int p, Float phi, Float theta_i) const {
    //        Vector3 ret_val(0.0, 0.0, 0.0);
    //        const int SEGMENTS = 50;

    //        Float dh = 2.0 / SEGMENTS;

    //        Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);
    //        for (int i = 0; i < SEGMENTS; i++) {
    //            Float h = (i + 0.5) / SEGMENTS * 2.0 - 1.0;
    //            ret_val += Ns(p, h, phi, theta_i, eta_prime) * dh;
    //        }

    //        ret_val *= 0.5;
    //        return ret_val;
    //    }

    Vector3 Ns(int p, Float h, Float phi, Float theta_i) const {
        Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);
        return Ns(p, h, phi, theta_i, eta_prime);
    }



    // Choose a lobe to sample.

    void chooseLobe(Float theta_i, Float h, int &lobe, Float &pdfLobe) const {
        Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);

        Float ASpec[5];
        ASpec[0] = Au(0, h, theta_i, eta_prime).length();
        ASpec[1] = Au(1, h, theta_i, eta_prime).length();
        ASpec[2] = Au(2, h, theta_i, eta_prime).length();
        ASpec[3] = As(1, h, M_PI / 2.0, theta_i, eta_prime).length();
        ASpec[4] = As(2, h, M_PI / 2.0, theta_i, eta_prime).length();

        Float sumPdf = 0;
        for (int i = 0; i < 5; i++)
            sumPdf += ASpec[i];
        if (sumPdf == 0.0) {
            lobe = -1;
            pdfLobe = 0.0;
            return;
        }
        for (int i = 0; i < 5; i++)
            ASpec[i] /= sumPdf;

        int choice = 4;
        Float randNum = sampler->next1D();
        for (int i = 0; i < 5; i++) {
            randNum -= ASpec[i];
            if (randNum <= 0.0) {
                choice = i;
                break;
            }
        }

        lobe = choice;
        pdfLobe = ASpec[choice];
    }

    void sampleMr(Float theta_i, Float phi, Float &theta_o, Float &pdfM, Float &weightM) const {
        Float theta_cone = theta_cone_R(theta_i, alpha, phi);
        Float beta_cone = beta_m_cone_R(beta_m, phi);

        theta_o = sampleGaussian(theta_cone, beta_cone, sampler->next1D(), sampler->next1D());
        pdfM = G(theta_o - theta_cone, beta_cone);
        weightM = 1.0;
    }

    void sampleMu(int p, Float theta_i, Float &theta_o, Float &pdfM, Float &weightM) const {
        Float alpha_p, beta_m_p;

        switch (p) {
            case 1: // TT
                alpha_p = -alpha / 2.0;
                beta_m_p = beta_m; // / 2.0;
                break;
            case 2: // TRT
                alpha_p = -alpha / 2.0 * 3.0;
                beta_m_p = beta_m / 2.0 * 3.0;
                break;
        }

        theta_o = sampleGaussian(-theta_i + 2.0 * alpha_p, beta_m_p, sampler->next1D(), sampler->next1D());
        pdfM = G(theta_o - (-theta_i + 2.0 * alpha_p), beta_m_p);
        weightM = 1.0;
    }

    void sampleNu(int p, Float h, Float theta_i, Float phi_i, Float &phi_o, Float &pdfN, Vector3 &weightN) const {
        Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);

        Float gamma_i = asinSafe(h);
        Float gamma_t = asinSafe(h / eta_prime);

        Vector3 attenuation = Au(p, h, theta_i, eta_prime);
        Float Phi = 2.0 * p * gamma_t - 2.0 * gamma_i + p * M_PI;

        Float beta_n_p_sqr = beta_n * beta_n * (p + 1);
        Float beta_n_p = std::sqrt(beta_n_p_sqr);

        Float phi = sampleGaussian(Phi, beta_n_p, sampler->next1D(), sampler->next1D());

        phi_o = phi_i + phi;
        pdfN = G(Phi - phi, beta_n_p);
        weightN = attenuation;
    }

    void sampleMs(Float theta_i, Float phi, Float &theta_o, Float &pdfM, Float &weightM) const {
        theta_o = asinSafe(sampler->next1D() * 2.0 - 1.0);
        pdfM = 0.5 * cos(theta_o);
        //        theta_o = (sampler->next1D() - 0.5) * M_PI;
        //        pdfM = 1.0 / M_PI;
        weightM = Ms(theta_i - alpha, theta_o - alpha, phi) / pdfM;
    }

    void sampleNs(int p, Float h, Float theta_i, Float phi_i, Float &phi_o, Float &pdfN, Vector3 &weightN) const {
        Float phi = (sampler->next1D() * 2.0 - 1.0) * M_PI;
        pdfN = 1.0 / (2.0 * M_PI);
        phi_o = phi_i + phi;
        weightN = max(discardLarger(Ns(p, h, phi, theta_i), 5.0), 0.0) / pdfN;
    }

    // Sample the chosen lobe. Return theta_o, phi_o, and corresponding sampling weights.

    void sampleLobe(int lobe, Float h, Float theta_i, Float phi_i, Float &theta_o, Float &pdfM, Float &weightM, Float &phi_o, Float &pdfN, Vector3 &weightN) const {
        if (lobe == 0) {
            // R lobe.
            sampleNu(lobe, h, theta_i, phi_i, phi_o, pdfN, weightN);
            Float phi = phi_o - phi_i;
            regularizePhi(phi);
            sampleMr(theta_i, phi, theta_o, pdfM, weightM);
        } else if (lobe <= 2) {
            // Unscattered TT and TRT lobes.
            sampleNu(lobe, h, theta_i, phi_i, phi_o, pdfN, weightN);
            sampleMu(lobe, theta_i, theta_o, pdfM, weightM);
        } else {
            // Scattered TT and TRT lobes.
            sampleNs(lobe - 2, h, theta_i, phi_i, phi_o, pdfN, weightN);
            Float phi = phi_o - phi_i;
            regularizePhi(phi);
            sampleMs(theta_i, phi, theta_o, pdfM, weightM);
        }
        //        regularizePhi(phi_o);
        //        regularizeTheta(theta_o);
    }

    // Longitudinal PDF of a chosen lobe

    inline Float pdfMLobe(int lobe, Float theta_i, Float theta_o, Float phi) const {
        if (lobe <= 2) {
            // Unscattered R, TT and TRT lobes.
            regularizePhi(phi);
            return Mu(lobe, theta_i, theta_o, phi);
        } else {
            // Scattered TT and TRT lobes.
            //            return 1.0 / M_PI;
            return 0.5 * cos(theta_o);
        }
    }

    // Azimuthal PDF of a chosen lobe

    inline Float pdfNLobe(int lobe, Float h, Float theta_i, Float phi) const {
        if (lobe <= 2) {
            // Unscattered R, TT and TRT lobes.
            int p = lobe;

            Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);

            Float gamma_i = asinSafe(h);
            Float gamma_t = asinSafe(h / eta_prime);

            Float Phi = 2.0 * p * gamma_t - 2.0 * gamma_i + p * M_PI;

            Float beta_n_p_sqr = beta_n * beta_n * (p + 1);
            Float beta_n_p = std::sqrt(beta_n_p_sqr);

            Float phiDiff = Phi - phi;
            regularizePhi(phiDiff);

            return Du(phiDiff, beta_n_p);
        } else {
            // Scattered TT and TRT lobes.
            return 1.0 / (2.0 * M_PI);
        }
    }

    Float pdfMN(Float h, Float theta_i, Float theta_o, Float phi) const {
        Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);

        Float ASpec[5];
        ASpec[0] = Au(0, h, theta_i, eta_prime).length();
        ASpec[1] = Au(1, h, theta_i, eta_prime).length();
        ASpec[2] = Au(2, h, theta_i, eta_prime).length();
        ASpec[3] = As(1, h, M_PI / 2.0, theta_i, eta_prime).length();
        ASpec[4] = As(2, h, M_PI / 2.0, theta_i, eta_prime).length();

        Float sumPdf = 0;
        for (int i = 0; i < 5; i++)
            sumPdf += ASpec[i];
        if (sumPdf == 0.0)
            return 0.0;
        for (int i = 0; i < 5; i++)
            ASpec[i] /= sumPdf;


        Float retVal = 0.0;
        for (int i = 0; i < 5; i++)
            retVal += pdfMLobe(i, theta_i, theta_o, phi) *
            pdfNLobe(i, h, theta_i, phi) *
            ASpec[i];
        retVal /= cos(theta_o);
        return retVal;
    }

    inline long double integrateLinear(long double a, long double b, long double c, long double d, long double e, long double f, long double h1, long double h2) const {
        // Common computations.
        long double exp1 = std::exp(-a * h1 + b);
        long double exp2 = std::exp(-a * h2 + b);
        long double a2 = a * a;
        long double a3 = a * a * a;
        long double a4 = a2 * a2;

        // Integrate cx^3*exp(-ax+b).
        long double p31 = -c * (a3 * h1 * h1 * h1 + 3.0 * a2 * h1 * h1 + 6.0 * a * h1 + 6) * exp1 / a4;
        long double p32 = -c * (a3 * h2 * h2 * h2 + 3.0 * a2 * h2 * h2 + 6.0 * a * h2 + 6) * exp2 / a4;

        // Integrate dx^2*exp(-ax+b).
        long double p21 = -d * (a2 * h1 * h1 + 2.0 * a * h1 + 2) * exp1 / a3;
        long double p22 = -d * (a2 * h2 * h2 + 2.0 * a * h2 + 2) * exp2 / a3;

        // Integrate ex*exp(-ax+b).
        long double p11 = -e * (a * h1 + 1) * exp1 / a2;
        long double p12 = -e * (a * h2 + 1) * exp2 / a2;

        // Integrate f*exp(-ax+b).
        long double p01 = -f * exp1 / a;
        long double p02 = -f * exp2 / a;

        return (p02 + p12 + p22 + p32) - (p01 + p11 + p21 + p31);
    }

    //    long double integrateLinear(long double a, long double b, long double d, long double e, long double f, long double h1, long double h2) const {
    //        if (std::fabs(a) < EPS * 5.0)
    //            return 0.0;

    //        long double a2 = a * a;
    //        long double a3 = std::pow(a, (long double) 3.0);
    //        long double t1a = d * (h1 * h1 / a + 2.0 * h1 / a2 + 2.0 / a3);
    //        long double t2a = (f + e * h1) / a + e / a2;
    //        long double t3a = -std::exp(-a * h1);

    //        long double t1b = d * (h2 * h2 / a + 2.0 * h2 / a2 + 2.0 / a3);
    //        long double t2b = (f + e * h2) / a + e / a2;
    //        long double t3b = -std::exp(-a * h2);

    //        long double retVal = std::exp(b) * (t3b * (t1b + t2b) - t3a * (t1a + t2a));

    //        return retVal;
    //    }


    //    long double integrateLinear(long double a, long double b, long double c, long double d, long double e, long double f, long double h1, long double h2) const {
    //        if (std::fabs(a) < EPS * 5.0)
    //            return 0.0;

    //        long double a4 = std::pow(a, (long double) 4.0);
    //        long double a3 = std::pow(a, (long double) 3.0);
    //        long double a2 = a * a;
    //        long double h13 = std::pow(h1, (long double) 3.0);
    //        long double h12 = h1 * h1;
    //        long double h23 = std::pow(h2, (long double) 3.0);
    //        long double h22 = h2 * h2;
    //        long double x1 = std::exp(-a * h1);
    //        long double x2 = std::exp(-a * h2);

    //        long double t = -1.0 * std::exp(b);

    //        long double aa = (x2 * h23 - x1 * h13) * c / a;
    //        long double bb = (x2 * h22 - x1 * h12) * d / a;
    //        long double cc = (x2 * h2 - x1 * h1) * e / a;
    //        long double dd = (x2 - x1) * f / a;
    //        long double ee = (x2 * h22 - x1 * h12) * 3.0 * c / a2;
    //        long double ff = (x2 * h2 - x1 * h1) * 2.0 * d / a2;
    //        long double gg = (x2 - x1) * e / a2;
    //        long double hh = (x2 * h2 - x1 * h1) * 6.0 * c / a3;
    //        long double ii = (x2 - x1) * 2.0 * d / a3;
    //        long double jj = (x2 - x1) * 6.0 * c / a4;

    //        return t * (aa + bb + cc + dd + ee + ff + gg + hh + ii + jj);
    //    }

    // Calculate 1/2 * int(A * D).
    // A = x * y *
    //     exp3(-((2.0 * p - 1.0) * s1 + (1.0 - kappa)) * sigma_a_prime) *
    //     std::exp(-2.0 * (p - 1) * s2 * sigma_e_m_prime) *
    //     std::exp(-avgDist * sigma_a_m_prime);
    // x = std::pow(f, p - 1.0);
    // y = (1.0 - f);

    Vector3 integrateS(int p, long double kappa, long double f1, long double f2, long double s11, long double s12, long double s21, long double s22, long double avgDist1, long double avgDist2, Vector3 sigma_a_prime, long double sigma_e_m_prime, long double sigma_a_m_prime, long double D1, long double D2, long double h1, long double h2) const {
        long double x1 = std::pow(f1, (long double) (p - 1.0));
        long double y1 = 1.0 - f1;

        long double x2 = std::pow(f2, (long double) (p - 1.0));
        long double y2 = 1.0 - f2;

        // Now compute each coefficient in representation
        // (c*x^3+d*x^2+e*x+f)*std::exp(-a*x+b)

        long double k_x = (x2 - x1) / (h2 - h1);
        long double b_x = x1 - h1 * k_x;

        long double k_y = (y2 - y1) / (h2 - h1);
        long double b_y = y1 - h1 * k_y;

        long double k_D = (D2 - D1) / (h2 - h1);
        long double b_D = D1 - h1 * k_D;

        long double c = k_x * k_y * k_D;
        long double d = k_x * k_y * b_D + k_x * b_y * k_D + b_x * k_y * k_D;
        long double e = k_x * b_y * b_D + b_x * b_y * k_D + b_x * k_y * b_D;
        ;
        long double f = b_x * b_y * b_D;


        long double k_s1 = (s12 - s11) / (h2 - h1);
        long double b_s1 = s11 - h1 * k_s1;

        long double k_s2 = (s22 - s21) / (h2 - h1);
        long double b_s2 = s21 - h1 * k_s2;

        long double k_avgDist = (avgDist2 - avgDist1) / (h2 - h1);
        long double b_avgDist = avgDist1 - h1 * k_avgDist;

        Vector3 retVal;
        for (int i = 0; i < 3; i++) {
            // Note there's already a negative symbol before a.
            long double a = (2.0 * p - 1.0) * k_s1 * sigma_a_prime[i] +
                    2.0 * (p - 1.0) * k_s2 * sigma_e_m_prime +
                    k_avgDist * sigma_a_m_prime;
            long double b = -(((2.0 * p - 1.0) * b_s1 + (1.0 - kappa)) * sigma_a_prime[i] +
                    2.0 * (p - 1.0) * b_s2 * sigma_e_m_prime +
                    b_avgDist * sigma_a_m_prime);

            retVal[i] = integrateLinear(a, b, c, d, e, f, h1, h2);
            //printf("%f %f %f %f %f %f %f %f\n", a, b, c, d, e, f, h1, h2);
            //            if (c != 0.0)
            //                retVal[i] = integrateLinear(a, b, c, d, e, f, h1, h2);
            //            else
            //                retVal[i] = integrateLinear(a, b, d, e, f, h1, h2);

            //            else {
            //                long double a2 = std::pow(a, (long double) 2.0);

            //                d /= a2;
            //                e = 2.0 * b * d / a2 + e / a;
            //                f += d * b * b / a2 + e * b / a;
            //                h1 = a * h1 - b;
            //                h2 = a * h2 - b;

            //                long double t = a;
            //                a = 1.0;
            //                b = 0.0;
            //                retVal[i] = integrateLinear(a, b, d, e, f, h1, h2) / t;


            //            retVal[i] = std::min(retVal[i], Float(5.0f));
            //            if (retVal[i] > 5.0)
            //                printf("%f %f %f %f %f %f %f %f\n%f %f %f %f %f %f\n\n",
            //                       a, b, c, d, e, f, h1, h2,
            //                       avgDist1, avgDist2, D1, D2, h1, h2);
            //            }
        }

        return Float(0.5) * retVal;
    }


    //    Vector3 Ns(int p, Float phi, Float theta_i) const {
    //        Vector3 retVal(0.0f, 0.0f, 0.0f);
    //        int SEGMENTS = 10;
    //        for (int i = 0; i < SEGMENTS; i++) {
    //            Float h = (i + sampler->next1D()) / SEGMENTS * 2.0 - 1.0;
    //            retVal += Ns(p, h, phi, theta_i) * (2.0 / SEGMENTS);
    //        }

    //        return Float(0.5f) * retVal;
    //    }

    Vector3 Ns(int p, Float phi, Float theta_i) const {
        if (kappa == 0.0)
            return Vector3(0.0, 0.0, 0.0);

        // Now using our analytic method.
        // Assume h is in [0, 1].
        const int SEGMENTS = 4;

        Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);
        Float hMax = min(Float(1.0), kappa * eta_prime) - EPS;

        Vector3 retVal(0.0, 0.0, 0.0);

        for (int i = 0; i < SEGMENTS; i++) {
            Float h1 = (i + 0.0) / SEGMENTS;
            Float h2 = (i + 1.0) / SEGMENTS;
            h1 = std::sqrt(h1);
            h2 = std::sqrt(h2);

            //            Float h1, h2;
            //            h1 = 1.0 - std::pow(1.5, -(i + 0.0));
            //            if (i < SEGMENTS - 1) {
            //                h2 = 1.0 - std::pow(1.5, -(i + 1.0));
            //            } else {
            //                h2 = 1.0;
            //            }

            h1 *= hMax;
            h2 *= hMax;

            Float gamma_i1 = asinSafe(h1);
            Float gamma_t1 = asinSafe(h1 / eta_prime);
            Float f1 = fresnelCuticle(gamma_i1, 1.0, eta_prime, layers);

            Float gamma_i2 = asinSafe(h2);
            Float gamma_t2 = asinSafe(h2 / eta_prime);
            Float f2 = fresnelCuticle(gamma_i2, 1.0, eta_prime, layers);

            Float theta_t = asinSafe(1.0 / eta * sin(theta_i));
            Vector3 sigma_a_prime = sigma_a / cos(theta_t);
            Float sigma_a_m_prime = sigma_a_m / cos(theta_t);
            Float sigma_s_m_prime = sigma_s_m / cos(theta_t);
            Float sigma_e_m_prime = sigma_a_m_prime + sigma_s_m_prime;

            // Travelling distance in the cortex and medulla, respectively.
            Float d11 = cos(gamma_t1), d21 = 0.0;
            Float h_prime1 = h1 / eta_prime;
            if (std::fabs(h_prime1) < std::fabs(kappa))
                d21 = std::sqrt(kappa * kappa - h_prime1 * h_prime1);
            d11 -= d21;
            Float d12 = cos(gamma_t2), d22 = 0.0;
            Float h_prime2 = h2 / eta_prime;
            if (std::fabs(h_prime2) < std::fabs(kappa))
                d22 = std::sqrt(kappa * kappa - h_prime2 * h_prime2);
            d12 -= d22;

            Float phi_query1, phi_query2, avgDist1, avgDist2, D1, D2;

            // Use phi.
            phi_query1 = phi - (gamma_t1 - gamma_i1) - (M_PI + 2.0 * gamma_t1) * (p - 1);
            regularizePhi(phi_query1);
            avgDist1 = kappa * getScatteredDist(sigma_s_m_prime * kappa, h_prime1 / kappa, g, phi_query1);
            D1 = getScatteredN(sigma_s_m_prime * kappa, h_prime1 / kappa, g, phi_query1);

            phi_query2 = phi - (gamma_t2 - gamma_i2) - (M_PI + 2.0 * gamma_t2) * (p - 1);
            regularizePhi(phi_query2);
            avgDist2 = kappa * getScatteredDist(sigma_s_m_prime * kappa, h_prime2 / kappa, g, phi_query2);
            D2 = getScatteredN(sigma_s_m_prime * kappa, h_prime2 / kappa, g, phi_query2);

            retVal += integrateS(p, kappa, f1, f2, d11, d12, d21, d22, avgDist1, avgDist2, sigma_a_prime, sigma_e_m_prime, sigma_a_m_prime, D1, D2, h1, h2);

            // Use -phi.
            phi_query1 = -phi - (gamma_t1 - gamma_i1) - (M_PI + 2.0 * gamma_t1) * (p - 1);
            regularizePhi(phi_query1);
            avgDist1 = kappa * getScatteredDist(sigma_s_m_prime * kappa, h_prime1 / kappa, g, phi_query1);
            D1 = getScatteredN(sigma_s_m_prime * kappa, h_prime1 / kappa, g, phi_query1);

            phi_query2 = -phi - (gamma_t2 - gamma_i2) - (M_PI + 2.0 * gamma_t2) * (p - 1);
            regularizePhi(phi_query2);
            avgDist2 = kappa * getScatteredDist(sigma_s_m_prime * kappa, h_prime2 / kappa, g, phi_query2);
            D2 = getScatteredN(sigma_s_m_prime * kappa, h_prime2 / kappa, g, phi_query2);

            retVal += integrateS(p, kappa, f1, f2, d11, d12, d21, d22, avgDist1, avgDist2, sigma_a_prime, sigma_e_m_prime, sigma_a_m_prime, D1, D2, h1, h2);
        }

        //        if (retVal[0] > 1.0f)
        //            printf("%f %f %f\n", retVal[0], retVal[1], retVal[2]);

        return retVal;
    }

    Vector3 NsMid(int p, Float hMinInput, Float hMaxInput, Float phi, Float theta_i) const {
        if (kappa == 0.0)
            return Vector3(0.0, 0.0, 0.0);

        // Now using our analytic method.
        // Assume h is in [0, 1].
        const int SEGMENTS = 4;

        Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);
        Float hMax = min(Float(1.0), kappa * eta_prime) - EPS;

        Vector3 retVal(0.0, 0.0, 0.0);

        for (int i = 0; i < SEGMENTS; i++) {
            Float h1 = (i + 0.0) / SEGMENTS;
            Float h2 = (i + 1.0) / SEGMENTS;
            h1 = std::sqrt(h1);
            h2 = std::sqrt(h2);

            h1 *= hMax;
            h2 *= hMax;


            if (h2 <= hMinInput || h1 >= hMaxInput)
                continue;

            //            printf("%f %f\n", h1, h2);


            h1 = max(h1, hMinInput);
            h2 = min(h2, hMaxInput);


            Float gamma_i1 = asinSafe(h1);
            Float gamma_t1 = asinSafe(h1 / eta_prime);
            Float f1 = fresnelCuticle(gamma_i1, 1.0, eta_prime, layers);

            Float gamma_i2 = asinSafe(h2);
            Float gamma_t2 = asinSafe(h2 / eta_prime);
            Float f2 = fresnelCuticle(gamma_i2, 1.0, eta_prime, layers);

            Float theta_t = asinSafe(1.0 / eta * sin(theta_i));
            Vector3 sigma_a_prime = sigma_a / cos(theta_t);
            Float sigma_a_m_prime = sigma_a_m / cos(theta_t);
            Float sigma_s_m_prime = sigma_s_m / cos(theta_t);
            Float sigma_e_m_prime = sigma_a_m_prime + sigma_s_m_prime;

            // Travelling distance in the cortex and medulla, respectively.
            Float d11 = cos(gamma_t1), d21 = 0.0;
            Float h_prime1 = h1 / eta_prime;
            if (std::fabs(h_prime1) < std::fabs(kappa))
                d21 = std::sqrt(kappa * kappa - h_prime1 * h_prime1);
            d11 -= d21;
            Float d12 = cos(gamma_t2), d22 = 0.0;
            Float h_prime2 = h2 / eta_prime;
            if (std::fabs(h_prime2) < std::fabs(kappa))
                d22 = std::sqrt(kappa * kappa - h_prime2 * h_prime2);
            d12 -= d22;

            Float phi_query1, phi_query2, avgDist1, avgDist2, D1, D2;

            // Use phi.
            phi_query1 = phi - (gamma_t1 - gamma_i1) - (M_PI + 2.0 * gamma_t1) * (p - 1);
            regularizePhi(phi_query1);
            avgDist1 = kappa * getScatteredDist(sigma_s_m_prime * kappa, h_prime1 / kappa, g, phi_query1);
            D1 = getScatteredN(sigma_s_m_prime * kappa, h_prime1 / kappa, g, phi_query1);

            phi_query2 = phi - (gamma_t2 - gamma_i2) - (M_PI + 2.0 * gamma_t2) * (p - 1);
            regularizePhi(phi_query2);
            avgDist2 = kappa * getScatteredDist(sigma_s_m_prime * kappa, h_prime2 / kappa, g, phi_query2);
            D2 = getScatteredN(sigma_s_m_prime * kappa, h_prime2 / kappa, g, phi_query2);

            retVal += integrateS(p, kappa, f1, f2, d11, d12, d21, d22, avgDist1, avgDist2, sigma_a_prime, sigma_e_m_prime, sigma_a_m_prime, D1, D2, h1, h2);

            //            // Use -phi.
            //            phi_query1 = -phi - (gamma_t1 - gamma_i1) - (M_PI + 2.0 * gamma_t1) * (p - 1);
            //            regularizePhi(phi_query1);
            //            avgDist1 = kappa * getScatteredDist(sigma_s_m_prime * kappa, h_prime1 / kappa, g, phi_query1);
            //            D1 = getScatteredN(sigma_s_m_prime * kappa, h_prime1 / kappa, g, phi_query1);

            //            phi_query2 = -phi - (gamma_t2 - gamma_i2) - (M_PI + 2.0 * gamma_t2) * (p - 1);
            //            regularizePhi(phi_query2);
            //            avgDist2 = kappa * getScatteredDist(sigma_s_m_prime * kappa, h_prime2 / kappa, g, phi_query2);
            //            D2 = getScatteredN(sigma_s_m_prime * kappa, h_prime2 / kappa, g, phi_query2);

            //            retVal += integrateS(p, kappa, f1, f2, d11, d12, d21, d22, avgDist1, avgDist2, sigma_a_prime, sigma_e_m_prime, sigma_a_m_prime, D1, D2, h1, h2);
        }

        return retVal / (hMaxInput - hMinInput) * 2.0;
    }

    inline long double erfMine(long double x) const {
        long double tau;

        long double t = 1.0 / (1.0 + 0.5 * std::fabs(x));
        long double k = -x * x -
                1.26551223 +
                1.00002368 * t +
                0.37409196 * std::pow(t, (long double) 2.0) +
                0.09678418 * std::pow(t, (long double) 3.0) -
                0.18628806 * std::pow(t, (long double) 4.0) +
                0.27886807 * std::pow(t, (long double) 5.0) -
                1.13520398 * std::pow(t, (long double) 6.0) +
                1.48851587 * std::pow(t, (long double) 7.0) -
                0.82215223 * std::pow(t, (long double) 8.0) +
                0.17087277 * std::pow(t, (long double) 9.0);

        tau = t * std::exp(k);

        if (x >= 0)
            return 1.0 - tau;
        else
            return tau - 1.0;
    }

    inline long double integrateR(long double f1, long double f2, long double d1, long double d2, long double h1, long double h2) const {
        const long double beta = beta_n;

        long double k_d = (d2 - d1) / (h2 - h1);
        long double b_d = d1 - h1 * k_d;

        long double k_A = (f2 - f1) / (h2 - h1);
        long double b_A = f1 - h1 * k_A;

        long double c = 1.0 / 2.0 / (std::sqrt(2.0 * M_PI) * beta);
        long double t = beta / (2.0 * k_d * k_d);

        long double erf_h2 = erfMine(d2 / (std::sqrt(2.0) * beta));
        long double exp_h2 = std::exp(-d2 * d2 / (2.0 * beta * beta));

        long double erf_h1 = erfMine(d1 / (std::sqrt(2.0) * beta));
        long double exp_h1 = std::exp(-d1 * d1 / (2.0 * beta * beta));

        long double I = t * (std::sqrt(2.0 * M_PI) * (erf_h2 - erf_h1) * (b_A * k_d - b_d * k_A) - 2.0 * k_A * beta * (exp_h2 - exp_h1));

        long double retVal = c * I;

        if (retVal < 0.0)
            return 0.0;

        //        if (retVal < 0.0)
        //            printf("BAD %Lf %Lf %Lf %Lf %Lf %Lf\n", f1, f2, d1, d2, h1, h2);

        if (std::isinf(retVal) || std::isinf(retVal) || std::isnan(retVal) || std::isnan(retVal))
            return 0.0;

        //        if (retVal > 0.1)
        //            printf("BAD %Lf: %Lf %Lf %Lf %Lf %Lf %Lf\n", retVal, f1, f2, d1, d2, h1, h2);

        return retVal;
    }

    // For R lobe only

    Float integrateR(Float f1, Float f2, Float Phi1, Float Phi2, Float phi, Float h1, Float h2) const {
        Float d1 = phi - Phi1;
        Float d2 = phi - Phi2;
        Float dd = d2 - d1;
        regularizePhi(d1);
        d2 = d1 + dd;

        if (d2 > M_PI) {
            return integrateR(f1, f2, d1, d2, h1, h2) +
                    integrateR(f1, f2, d1 - M_PI * 2.0, d2 - M_PI * 2.0, h1, h2);
        } else if (d2 < -M_PI) {
            return integrateR(f1, f2, d1, d2, h1, h2) +
                    integrateR(f1, f2, d1 + M_PI * 2.0, d2 + M_PI * 2.0, h1, h2);
        }
        return integrateR(f1, f2, d1, d2, h1, h2);
    }

    inline long double integrateQuad(long double a, long double b, long double c, long double d, long double e, long double f, long double h1, long double h2) const {
        if (a <= 0.0)
            return 0.0;

        if (std::fabs((2.0 * a * h1 - b) / (2.0 * std::sqrt(a))) > 5.0 &&
                std::fabs((2.0 * a * h2 - b) / (2.0 * std::sqrt(a))) > 5.0)
            return 0.0;


        long double x1 = erfMine((2.0 * a * h1 - b) / (2.0 * std::sqrt(a)));
        long double x2 = erfMine((2.0 * a * h2 - b) / (2.0 * std::sqrt(a)));

        long double y1 = std::exp(-a * h1 * h1 + b * h1 + c) * (2.0 * a * (d * h1 + e) + b * d);
        long double y2 = std::exp(-a * h2 * h2 + b * h2 + c) * (2.0 * a * (d * h2 + e) + b * d);

        long double p = std::sqrt(M_PI) / (8.0 * std::pow(a, (long double) 2.5)) * std::exp(c + (b * b) / (4.0 * a)) * (4.0 * a * a * f + 2.0 * a * (b * e + d) + b * b * d);
        long double q = 1.0 / (4.0 * a * a);

        if (std::isinf(p) || std::isinf(q) || std::isnan(p) || std::isnan(q))
            return 0.0;
        long double retVal = (p * x2 - q * y2) - (p * x1 - q * y1);


        return retVal < 0.0 ? 0.0 : retVal;
    }

    inline long double integrateQuad(long double a, long double b, long double c, long double e, long double f, long double h1, long double h2) const {
        if (a <= 0.0)
            return 0.0;

        if (std::fabs((2.0 * a * h1 - b) / (2.0 * std::sqrt(a))) > 5.0 &&
                std::fabs((2.0 * a * h2 - b) / (2.0 * std::sqrt(a))) > 5.0)
            return 0.0;


        long double x1 = erfMine((2.0 * a * h1 - b) / (2.0 * std::sqrt(a)));
        long double x2 = erfMine((2.0 * a * h2 - b) / (2.0 * std::sqrt(a)));

        long double y1 = std::exp(-a * h1 * h1 + b * h1 + c) * (2.0 * a * e);
        long double y2 = std::exp(-a * h2 * h2 + b * h2 + c) * (2.0 * a * e);

        long double p = std::sqrt(M_PI) / (8.0 * std::pow(a, (long double) 2.5)) * std::exp(c + (b * b) / (4.0 * a)) * (4.0 * a * a * f + 2.0 * a * b * e);
        long double q = 1.0 / (4.0 * a * a);

        if (std::isinf(p) || std::isinf(q) || std::isnan(p) || std::isnan(q))
            return 0.0;
        long double retVal = (p * x2 - q * y2) - (p * x1 - q * y1);


        return retVal < 0.0 ? 0.0 : retVal;
    }

    inline long double integrateQuad(long double a, long double b, long double c, long double f, long double h1, long double h2) const {
        if (a <= 0.0)
            return 0.0;

        if (std::fabs((2.0 * a * h1 - b) / (2.0 * std::sqrt(a))) > 5.0 &&
                std::fabs((2.0 * a * h2 - b) / (2.0 * std::sqrt(a))) > 5.0)
            return 0.0;


        long double x1 = erfMine((2.0 * a * h1 - b) / (2.0 * std::sqrt(a)));
        long double x2 = erfMine((2.0 * a * h2 - b) / (2.0 * std::sqrt(a)));

        long double y1 = 0.0;
        long double y2 = 0.0;

        long double p = std::sqrt(M_PI) / (8.0 * std::sqrt(a)) * std::exp(c + (b * b) / (4.0 * a)) * (4.0 * f);
        long double q = 1.0 / (4.0 * a * a);

        if (std::isinf(p) || std::isinf(q) || std::isnan(p) || std::isnan(q))
            return 0.0;
        long double retVal = (p * x2 - q * y2) - (p * x1 - q * y1);


        return retVal < 0.0 ? 0.0 : retVal;
    }

    inline Float integrateU(Float beta_p, int p, Float d, Float e, Float f, Float k_s1, Float b_s1, Float k_s2, Float b_s2, Float d1, Float d2, Float sigma_aa, Float sigma_ee, Float h1, Float h2) const {
        // Judge if \phi is between [\Phi_{min} - 3 \beta, \Phi_{max} + 3 \beta].
        Float dMin = std::min(d1, d2) - 3.0 * beta_p;
        Float dMax = std::max(d1, d2) + 3.0 * beta_p;
        if (!(dMin <= 0.0 && dMax >= 0.0))
            return 0.0;

        // Now calculate a, b, c;
        Float cc = 1.0 / 2.0 / (std::sqrt(2.0 * M_PI) * beta_p);

        Float k_d = (d2 - d1) / (h2 - h1);
        Float b_d = d1 - h1 * k_d;

        //        k_d = std::min(1000.0, std::max(1e-3, k_d));

        Float a = k_d * k_d / (2.0 * beta_p * beta_p);
        Float b = -(k_s1 * 2.0 * p * sigma_aa + k_s2 * 2.0 * p * sigma_ee + k_d * b_d / (beta_p * beta_p));
        Float c = -(b_s1 * 2.0 * p * sigma_aa + b_s2 * 2.0 * p * sigma_ee + b_d * b_d / (2.0 * beta_p * beta_p));

        //        if (p == 2 && h1 < 0.01) {
        //            printf("%f %f %f %f %f %f\n", d1, d2, h1, h2, k_d, b_d);
        //            printf("%f %f %f %f %f %f %f %f\n", a, b, c, d, e, f, h1, h2);
        //        }
        // Now calculate cc * int[(d*x^2+e*x+f)*std::exp(-a*x^2+b*x+c)]
        //        h2 = std::min(h2, (Float) 0.9);
        //        return cc * integrateQuad(a, b, c, d, e, f, h1, h2);

        //        Float t = 2.0 * b / (a * (h1 + h2));
        //        a *= t * t;
        //        b *= t;
        //        d *= t * t;
        //        e *= t;
        //        h1 /= t;
        //        h2 /= t;
        //        cc *= t;

        //        if (a < 1e-4)
        //            return 0.0;
        if (a < 0.0) {
            //            printf("%f %f %f %f %f %f\n", d1, d2, h1, h2, k_d, b_d);
            //            printf("%f %f %f %f %f %f %f %f\n", a, b, c, d, e, f, h1, h2);
            return 0.0;
        }

        if (d == 0.0 && e == 0.0)
            return cc * integrateQuad(a, b, c, f, h1, h2);

        if (d == 0.0)
            return cc * integrateQuad(a, b, c, e, f, h1, h2);

        //        Float t = std::sqrt(2.0) * beta_p * (h2 - h1) / (d2 - d1);
        //        Float tInv = k_d / (std::sqrt(2.0) * beta_p);
        //        a = 1.0;
        //        b *= t;
        //        d *= t * t;
        //        e *= t;
        //        h1 *= tInv;
        //        h2 *= tInv;
        //        cc *= t;

        //        Float t = 2.0 * b / (a * (h1 + h2));
        //        a *= t * t;
        //        b *= t;
        //        d *= t * t;
        //        e *= t;
        //        h1 /= t;
        //        h2 /= t;
        //        cc *= t;

        return cc * integrateQuad(a, b, c, d, e, f, h1, h2);
    }

    // Calculate 1/2 * int(A * D).
    // A = x * y *
    //     exp3(-2.0 * p * d1 * sigma_a_prime) *
    //     std::exp(-2.0 * p * d2 * sigma_e_m_prime);
    // D = G(Phi - phi);
    // x = std::pow(1.0 - f, 2.0);
    // y = std::pow(f, p - 1.0);

    Vector3 integrateU(int p, Float f1, Float f2, Float s11, Float s12, Float s21, Float s22, Vector3 sigma_a_prime, Float sigma_e_m_prime, Float Phi1, Float Phi2, Float phi, Float h1, Float h2) const {
        Float x1 = std::pow(f1, p - 1.0);
        Float y1 = std::pow(1.0 - f1, 2.0);

        Float x2 = std::pow(f2, p - 1.0);
        Float y2 = std::pow(1.0 - f2, 2.0);

        // Now compute each coefficient in representation
        // (d*x^2+e*x+f)*std::exp(-a*x^2+b*x+c)

        Float k_x = (x2 - x1) / (h2 - h1);
        Float b_x = x1 - h1 * k_x;

        Float k_y = (y2 - y1) / (h2 - h1);
        Float b_y = y1 - h1 * k_y;

        Float d = k_x * k_y;
        Float e = k_x * b_y + k_y * b_x;
        Float f = b_x * b_y;

        Float k_s1 = (s12 - s11) / (h2 - h1);
        Float b_s1 = s11 - h1 * k_s1;

        Float k_s2 = (s22 - s21) / (h2 - h1);
        Float b_s2 = s21 - h1 * k_s2;

        //        if (x1 == 1.0 && x2 == 1.0)
        //            printf("%f %f %f %f\n", f1, f2, h1, h2);
        //        printf("%f %f %f %f %f %f %f %f\n", x1, y1, x2, y2, k_x, b_x, k_y, b_y);


        const Float &beta = beta_n;
        Float beta_p_sqr = beta * beta * (p + 1);
        Float beta_p = std::sqrt(beta_p_sqr);

        Float d1 = phi - Phi1;
        Float d2 = phi - Phi2;
        Float dd = d2 - d1;
        regularizePhi(d1);
        d2 = d1 + dd;

        if (d2 > M_PI) {
            return Vector3(integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1, d2, sigma_a_prime[0], sigma_e_m_prime, h1, h2) +
                    integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1 - M_PI * 2.0, d2 - M_PI * 2.0, sigma_a_prime[0], sigma_e_m_prime, h1, h2),
                    integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1, d2, sigma_a_prime[1], sigma_e_m_prime, h1, h2) +
                    integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1 - M_PI * 2.0, d2 - M_PI * 2.0, sigma_a_prime[1], sigma_e_m_prime, h1, h2),
                    integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1, d2, sigma_a_prime[2], sigma_e_m_prime, h1, h2) +
                    integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1 - M_PI * 2.0, d2 - M_PI * 2.0, sigma_a_prime[2], sigma_e_m_prime, h1, h2));
        } else if (d2 < -M_PI) {
            return Vector3(integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1, d2, sigma_a_prime[0], sigma_e_m_prime, h1, h2) +
                    integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1 + M_PI * 2.0, d2 + M_PI * 2.0, sigma_a_prime[0], sigma_e_m_prime, h1, h2),
                    integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1, d2, sigma_a_prime[1], sigma_e_m_prime, h1, h2) +
                    integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1 + M_PI * 2.0, d2 + M_PI * 2.0, sigma_a_prime[1], sigma_e_m_prime, h1, h2),
                    integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1, d2, sigma_a_prime[2], sigma_e_m_prime, h1, h2) +
                    integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1 + M_PI * 2.0, d2 + M_PI * 2.0, sigma_a_prime[2], sigma_e_m_prime, h1, h2));
        }
        return Vector3(integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1, d2, sigma_a_prime[0], sigma_e_m_prime, h1, h2),
                integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1, d2, sigma_a_prime[1], sigma_e_m_prime, h1, h2),
                integrateU(beta_p, p, d, e, f, k_s1, b_s1, k_s2, b_s2, d1, d2, sigma_a_prime[2], sigma_e_m_prime, h1, h2));
    }

    Vector3 Nu(int p, Float phi, Float theta_i) const {
        // Now using our analytic method.
        // Assume h is in [0, 1].
        const int SEGMENTS = 5;

        Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);

        Vector3 retVal(0.0, 0.0, 0.0);

        for (int i = 0; i < SEGMENTS; i++) {
            Float h1 = (i + 0.0) / SEGMENTS;
            Float h2 = (i + 1.0) / SEGMENTS;
            h1 *= 1.0 - EPS;
            h2 *= 1.0 - EPS;
            h1 = std::sqrt(h1);
            h2 = std::sqrt(h2);

            Float gamma_i1 = asinSafe(h1);
            Float gamma_t1 = asinSafe(h1 / eta_prime);
            Float f1 = fresnelCuticle(gamma_i1, 1.0, eta_prime, layers);
            Float Phi1 = 2.0 * p * gamma_t1 - 2.0 * gamma_i1 + p * M_PI;

            Float gamma_i2 = asinSafe(h2);
            Float gamma_t2 = asinSafe(h2 / eta_prime);
            Float f2 = fresnelCuticle(gamma_i2, 1.0, eta_prime, layers);
            Float Phi2 = 2.0 * p * gamma_t2 - 2.0 * gamma_i2 + p * M_PI;

            if (p == 0) {
                retVal += Vector3(1.0, 1.0, 1.0) * integrateR(f1, f2, Phi1, Phi2, phi, h1, h2);
                retVal += Vector3(1.0, 1.0, 1.0) * integrateR(f1, f2, Phi1, Phi2, -phi, h1, h2);
            } else {
                Float theta_t = asinSafe(1.0 / eta * sin(theta_i));
                Vector3 sigma_a_prime = sigma_a / cos(theta_t);
                Float sigma_a_m_prime = sigma_a_m / cos(theta_t);
                Float sigma_s_m_prime = sigma_s_m / cos(theta_t);
                Float sigma_e_m_prime = sigma_a_m_prime + sigma_s_m_prime;

                // Travelling distance in the cortex and medulla, respectively.
                Float d11 = cos(gamma_t1), d21 = 0.0;
                Float h_prime1 = h1 / eta_prime;
                if (std::fabs(h_prime1) < std::fabs(kappa))
                    d21 = std::sqrt(kappa * kappa - h_prime1 * h_prime1);
                d11 -= d21;
                Float d12 = cos(gamma_t2), d22 = 0.0;
                Float h_prime2 = h2 / eta_prime;
                if (std::fabs(h_prime2) < std::fabs(kappa))
                    d22 = std::sqrt(kappa * kappa - h_prime2 * h_prime2);
                d12 -= d22;

                retVal += integrateU(p, f1, f2, d11, d12, d21, d22, sigma_a_prime, sigma_e_m_prime, Phi1, Phi2, phi, h1, h2);
                retVal += integrateU(p, f1, f2, d11, d12, d21, d22, sigma_a_prime, sigma_e_m_prime, Phi1, Phi2, -phi, h1, h2);
            }
        }
        if (std::isnan(retVal[0]) || std::isnan(retVal[1]) || std::isnan(retVal[2]))
            retVal = Vector3(0.0f, 0.0f, 0.0f);
        if (std::isinf(retVal[0]) || std::isinf(retVal[1]) || std::isinf(retVal[2]))
            retVal = Vector3(0.0f, 0.0f, 0.0f);
        for (int i = 0; i < 3; i++)
            retVal[i] = std::min(std::max(Float(0.0f), retVal[i]), Float(1e3f));

        return retVal;
    }

    Vector3 NuMid(int p, Float hMinInput, Float hMaxInput, Float phi, Float theta_i) const {
        // Now using our analytic method.
        // Assume h is in [0, 1].
        //        Assert(hMinInput >= 0.0 && hMaxInput >= 0.0);

        const int SEGMENTS = 5;

        Float eta_prime = std::sqrt(eta * eta - sin(theta_i) * sin(theta_i)) / cos(theta_i);

        Vector3 retVal(0.0, 0.0, 0.0);

        for (int i = 0; i < SEGMENTS; i++) {
            Float h1 = (i + 0.0) / SEGMENTS;
            Float h2 = (i + 1.0) / SEGMENTS;
            h1 *= 1.0 - EPS;
            h2 *= 1.0 - EPS;
            h1 = std::sqrt(h1);
            h2 = std::sqrt(h2);

            if (h2 <= hMinInput || h1 >= hMaxInput)
                continue;

            h1 = max(h1, hMinInput);
            h2 = min(h2, hMaxInput);

            //            printf("%f %f\n", h1, h2);

            Float gamma_i1 = asinSafe(h1);
            Float gamma_t1 = asinSafe(h1 / eta_prime);
            Float f1 = fresnelCuticle(gamma_i1, 1.0, eta_prime, layers);
            Float Phi1 = 2.0 * p * gamma_t1 - 2.0 * gamma_i1 + p * M_PI;

            Float gamma_i2 = asinSafe(h2);
            Float gamma_t2 = asinSafe(h2 / eta_prime);
            Float f2 = fresnelCuticle(gamma_i2, 1.0, eta_prime, layers);
            Float Phi2 = 2.0 * p * gamma_t2 - 2.0 * gamma_i2 + p * M_PI;

            if (p == 0) {
                retVal += Vector3(1.0, 1.0, 1.0) * integrateR(f1, f2, Phi1, Phi2, phi, h1, h2);
                //                retVal += Vector3(1.0, 1.0, 1.0) * integrateR(f1, f2, Phi1, Phi2, -phi, h1, h2);
            } else {
                Float theta_t = asinSafe(1.0 / eta * sin(theta_i));
                Vector3 sigma_a_prime = sigma_a / cos(theta_t);
                Float sigma_a_m_prime = sigma_a_m / cos(theta_t);
                Float sigma_s_m_prime = sigma_s_m / cos(theta_t);
                Float sigma_e_m_prime = sigma_a_m_prime + sigma_s_m_prime;

                // Travelling distance in the cortex and medulla, respectively.
                Float d11 = cos(gamma_t1), d21 = 0.0;
                Float h_prime1 = h1 / eta_prime;
                if (std::fabs(h_prime1) < std::fabs(kappa))
                    d21 = std::sqrt(kappa * kappa - h_prime1 * h_prime1);
                d11 -= d21;
                Float d12 = cos(gamma_t2), d22 = 0.0;
                Float h_prime2 = h2 / eta_prime;
                if (std::fabs(h_prime2) < std::fabs(kappa))
                    d22 = std::sqrt(kappa * kappa - h_prime2 * h_prime2);
                d12 -= d22;

                retVal += integrateU(p, f1, f2, d11, d12, d21, d22, sigma_a_prime, sigma_e_m_prime, Phi1, Phi2, phi, h1, h2);
                //                retVal += integrateU(p, f1, f2, d11, d12, d21, d22, sigma_a_prime, sigma_e_m_prime, Phi1, Phi2, -phi, h1, h2);
            }
        }
        if (std::isnan(retVal[0]) || std::isnan(retVal[1]) || std::isnan(retVal[2]))
            retVal = Vector3(0.0f, 0.0f, 0.0f);
        if (std::isinf(retVal[0]) || std::isinf(retVal[1]) || std::isinf(retVal[2]))
            retVal = Vector3(0.0f, 0.0f, 0.0f);
        for (int i = 0; i < 3; i++)
            retVal[i] = std::min(std::max(Float(0.0f), retVal[i]), Float(1e3f));

        //        for (int i = 0; i < 3; i++)
        //            printf("%f %f %f\n", retVal[0], retVal[1], retVal[2]);

        return retVal / (hMaxInput - hMinInput) * 2.0;
    }

};

#define MASK_MEDULLA

void precomputeDS(const char *filename,
        Float eta, Vector3 sigma_a, Float alpha, Float beta_m, Float beta_n,
        Float kappa, Float sigma_a_m, Float sigma_s_m, Float g, Float layers) {

    Vector3 alpha_b_pre[NUM_BINS_DS], alpha_f_pre[NUM_BINS_DS];
    Vector3 beta_b_sqr_pre[NUM_BINS_DS], beta_f_sqr_pre[NUM_BINS_DS];
    Vector3 a_b_pre[NUM_BINS_DS], a_f_pre[NUM_BINS_DS];
    Vector3 a_bu_pre[NUM_BINS_DS], a_fu_pre[NUM_BINS_DS];
    Vector3 a_bs_pre[NUM_BINS_DS], a_fs_pre[NUM_BINS_DS];
    Vector3 A_b_pre[NUM_BINS_DS];
    Vector3 Delta_b_pre[NUM_BINS_DS];
    Vector3 sigma_b_pre_sqr[NUM_BINS_DS];
    Vector3 N_R_G_pre[NUM_BINS_DS];
    Vector3 N_TT_G_pre[NUM_BINS_DS];
    Vector3 N_TRT_G_pre[NUM_BINS_DS];
    Vector3 N_TTs_G_pre[NUM_BINS_DS];
    Vector3 N_TRTs_G_pre[NUM_BINS_DS];

    Vector3 abR[NUM_BINS_DS], afR[NUM_BINS_DS];
    Vector3 abTT[NUM_BINS_DS], afTT[NUM_BINS_DS];
    Vector3 abTRT[NUM_BINS_DS], afTRT[NUM_BINS_DS];
    Vector3 abTTs[NUM_BINS_DS], afTTs[NUM_BINS_DS];
    Vector3 abTRTs[NUM_BINS_DS], afTRTs[NUM_BINS_DS];

    for (int batch = 0; batch < NUM_SIGMA_SAMPLES / 3; batch++) {
        printf("Precomputing batch %d\n", batch);

        Float sigma_a1 = ((batch * 3 + 0) + 0.5) / NUM_SIGMA_SAMPLES * MAX_ABSORPTION_OUTER;
        Float sigma_a2 = ((batch * 3 + 1) + 0.5) / NUM_SIGMA_SAMPLES * MAX_ABSORPTION_OUTER;
        Float sigma_a3 = ((batch * 3 + 2) + 0.5) / NUM_SIGMA_SAMPLES * MAX_ABSORPTION_OUTER;
        Vector3 sigma_a_pre(sigma_a1, sigma_a2, sigma_a3);

        // Precomputations for dual scattering.
        FurSolver furSolver(eta, sigma_a_pre, alpha, beta_m, beta_n, kappa, sigma_a_m, sigma_s_m, g, layers);

#pragma omp parallel for schedule(dynamic)
        for (int k = 0; k < NUM_BINS_DS; k++) {
            Float theta_i = (k + 0.5) / NUM_BINS_DS * M_PI - M_PI / 2.0;

            // Calculate alpha and beta.
            Float alpha_R = alpha;
            Float alpha_TT = -alpha / 2.0;
            Float alpha_TRT = -alpha / 2.0 * 3.0;
            Float alpha_TTs = alpha_TT;
            Float alpha_TRTs = alpha_TRT;

            Float beta_R = beta_m;
            Float beta_TT = beta_m;
            Float beta_TRT = beta_m / 2.0 * 3.0;
            Float beta_TTs, beta_TRTs;

            // Calculate beta_TRTs and beta_TTs.

            Float EX, EX2, P;
            EX = EX2 = P = 0.0;
            for (int i = 0; i < NUM_THETA_SAMPLES; i++) {
                Float theta_r = (i + 0.5) / NUM_THETA_SAMPLES * M_PI - M_PI / 2.0;

                Float theta_i_compute = theta_i;
                Float theta_r_compute = theta_r;
                if (theta_i_compute < 0.0) {
                    theta_i_compute = -theta_i_compute;
                    theta_r_compute = -theta_r_compute;
                }

                Float binSigma = sigma_s_m * kappa / MAX_SCATTERING_INNER * NUM_SCATTERING_INNER;
                Float binThetaI = theta_i_compute / MAX_THETA * NUM_THETA;
                Float binG = g / MAX_G * NUM_G;
                Float binThetaRBack = (theta_r_compute + M_PI / 2.0) / M_PI * (NUM_BINS / 2 - 1);
                Float p = interpScatteredM(binSigma, binThetaI, binG, binThetaRBack);

                EX += theta_r * p;
                EX2 += theta_r * theta_r * p;
                P += p;
            }
            EX /= P;
            EX2 /= P;
            beta_TRTs = std::sqrt(EX2 - EX * EX);

            EX = EX2 = P = 0.0;
            for (int i = 0; i < NUM_THETA_SAMPLES; i++) {
                Float theta_r = (i + 0.5) / NUM_THETA_SAMPLES * M_PI - M_PI / 2.0;

                Float theta_i_compute = theta_i;
                Float theta_r_compute = theta_r;
                if (theta_i_compute < 0.0) {
                    theta_i_compute = -theta_i_compute;
                    theta_r_compute = -theta_r_compute;
                }

                Float binSigma = sigma_s_m * kappa / MAX_SCATTERING_INNER * NUM_SCATTERING_INNER;
                Float binThetaI = theta_i_compute / MAX_THETA * NUM_THETA;
                Float binG = g / MAX_G * NUM_G;
                Float binThetaRBack = (theta_r_compute + M_PI / 2.0) / M_PI * (NUM_BINS / 2 - 1);
                Float binThetaRFront = NUM_BINS - 1 - binThetaRBack;
                Float p = interpScatteredM(binSigma, binThetaI, binG, binThetaRFront);

                EX += theta_r * p;
                EX2 += theta_r * theta_r * p;
                P += p;
            }
            EX /= P;
            EX2 /= P;
            beta_TTs = std::sqrt(EX2 - EX * EX);

            //            printf("%f %f, %f\n", beta_TRTs, beta_TTs, M_PI / sqrt(12.0));


            //            beta_TTs = beta_TRTs = M_PI / sqrt(12.0);

#if defined(MASK_MEDULLA)
            beta_TTs = beta_TRTs = alpha_TTs = alpha_TRTs = 0.f;
#endif


            Float dTheta = M_PI / NUM_THETA_SAMPLES;
            Float dPhiR = M_PI / NUM_PHI_SAMPLES;
            Float dPhiI = M_PI / NUM_K_SAMPLES;

            // Back hemisphere
            abR[k] = Vector3(0.0f, 0.0f, 0.0f);
            abTT[k] = Vector3(0.0f, 0.0f, 0.0f);
            abTRT[k] = Vector3(0.0f, 0.0f, 0.0f);
            abTTs[k] = Vector3(0.0f, 0.0f, 0.0f);
            abTRTs[k] = Vector3(0.0f, 0.0f, 0.0f);

            for (int j = 0; j < NUM_PHI_SAMPLES; j++) {
                Float phi_r = (j + 0.5) / NUM_PHI_SAMPLES * M_PI - M_PI / 2.0;
                for (int l = 0; l < NUM_K_SAMPLES; l++) {
                    Float phi_i = (l + 0.5) / NUM_K_SAMPLES * M_PI - M_PI / 2.0;
                    Float phi = phi_r - phi_i;
                    regularizePhi(phi);

                    abR[k] += furSolver.Nu(0, phi, theta_i);
                    abTT[k] += furSolver.Nu(1, phi, theta_i);
                    abTRT[k] += furSolver.Nu(2, phi, theta_i);
                    abTTs[k] += max(discardLarger(furSolver.Ns(1, phi, theta_i), 10.0f), 0.0f) * SCATTERED_ATTENUATION;
                    abTRTs[k] += max(discardLarger(furSolver.Ns(2, phi, theta_i), 10.0f), 0.0f) * SCATTERED_ATTENUATION;
                }
            }

            abR[k] *= dPhiR * dPhiI / M_PI;
            abTT[k] *= dPhiR * dPhiI / M_PI;
            abTRT[k] *= dPhiR * dPhiI / M_PI;
            abTTs[k] *= dPhiR * dPhiI / M_PI;
            abTRTs[k] *= dPhiR * dPhiI / M_PI;

            a_bu_pre[k] = abR[k] + abTT[k] + abTRT[k];
            a_bs_pre[k] = abTTs[k] + abTRTs[k];
            a_bu_pre[k] = clamp(a_bu_pre[k], EPS, 1.0f - EPS);
            a_bs_pre[k] = clamp(a_bs_pre[k], EPS, 1.0f - EPS);
            a_b_pre[k] = clamp(a_bs_pre[k] + a_bu_pre[k], EPS, 1.0f - EPS);


            // Front hemisphere
            afR[k] = Vector3(0.0f, 0.0f, 0.0f);
            afTT[k] = Vector3(0.0f, 0.0f, 0.0f);
            afTRT[k] = Vector3(0.0f, 0.0f, 0.0f);
            afTTs[k] = Vector3(0.0f, 0.0f, 0.0f);
            afTRTs[k] = Vector3(0.0f, 0.0f, 0.0f);

            for (int j = 0; j < NUM_PHI_SAMPLES; j++) {
                Float phi_r = (j + 0.5) / NUM_PHI_SAMPLES * M_PI + M_PI / 2.0;
                for (int l = 0; l < NUM_K_SAMPLES; l++) {
                    Float phi_i = (l + 0.5) / NUM_K_SAMPLES * M_PI - M_PI / 2.0;
                    Float phi = phi_r - phi_i;
                    regularizePhi(phi);

                    afR[k] += furSolver.Nu(0, phi, theta_i);
                    afTT[k] += furSolver.Nu(1, phi, theta_i);
                    afTRT[k] += furSolver.Nu(2, phi, theta_i);
                    afTTs[k] += max(discardLarger(furSolver.Ns(1, phi, theta_i), 10.0f), 0.0f) * SCATTERED_ATTENUATION;
                    afTRTs[k] += max(discardLarger(furSolver.Ns(2, phi, theta_i), 10.0f), 0.0f) * SCATTERED_ATTENUATION;
                }
            }

            afR[k] *= dPhiR * dPhiI / M_PI;
            afTT[k] *= dPhiR * dPhiI / M_PI;
            afTRT[k] *= dPhiR * dPhiI / M_PI;
            afTTs[k] *= dPhiR * dPhiI / M_PI;
            afTRTs[k] *= dPhiR * dPhiI / M_PI;

            a_fu_pre[k] = afR[k] + afTT[k] + afTRT[k];
            a_fs_pre[k] = afTTs[k] + afTRTs[k];
            a_fu_pre[k] = clamp(a_fu_pre[k], EPS, 1.0f - EPS);
            a_fs_pre[k] = clamp(a_fs_pre[k], EPS, 1.0f - EPS);
            a_f_pre[k] = clamp(a_fs_pre[k] + a_fu_pre[k], EPS, 1.0f - EPS);

            alpha_b_pre[k] = (abR[k] * alpha_R + abTT[k] * alpha_TT + abTRT[k] * alpha_TRT +
                    abTTs[k] * alpha_TTs + abTRTs[k] * alpha_TRTs) / a_b_pre[k];
            beta_b_sqr_pre[k] = (abR[k] * beta_R + abTT[k] * beta_TT + abTRT[k] * beta_TRT +
                    abTTs[k] * beta_TTs + abTRTs[k] * beta_TRTs) / a_b_pre[k];
            beta_b_sqr_pre[k] *= beta_b_sqr_pre[k];


            alpha_f_pre[k] = (afR[k] * alpha_R + afTT[k] * alpha_TT + afTRT[k] * alpha_TRT +
                    afTTs[k] * alpha_TTs + afTRTs[k] * alpha_TRTs) / a_f_pre[k];
            beta_f_sqr_pre[k] = (afR[k] * beta_R + afTT[k] * beta_TT + afTRT[k] * beta_TRT +
                    afTTs[k] * beta_TTs + afTRTs[k] * beta_TRTs) / a_f_pre[k];
            beta_f_sqr_pre[k] *= beta_f_sqr_pre[k];

            //            printf("%f %f %f\n", beta_b_sqr_pre[k][0], beta_b_sqr_pre[k][1], beta_b_sqr_pre[k][2]);

            // Calculate A_b_pre.
            Vector3 A_1, A_3;
            for (int i = 0; i < 3; i++) {
                A_1[i] = a_b_pre[k][i] * std::pow(a_f_pre[k][i], 2.0) / (1.0 - std::pow(a_f_pre[k][i], 2.0));
                A_3[i] = std::pow(a_b_pre[k][i], 3.0) * std::pow(a_f_pre[k][i], 2.0) / std::pow(1.0 - std::pow(a_f_pre[k][i], 2.0), 3.0);
            }
            A_b_pre[k] = A_1 + A_3;

            // Calculate Delta_b_pre.
            Delta_b_pre[k] = Vector3(0.0);
            for (int i = 0; i < 3; i++) {
                Delta_b_pre[k][i] += alpha_b_pre[k][i] * (1.0 - 2.0 * std::pow(a_b_pre[k][i], 2.0) / std::pow(1.0 - std::pow(a_f_pre[k][i], 2.0), 2.0));
                Delta_b_pre[k][i] += alpha_f_pre[k][i] * (2.0 * std::pow(1.0 - std::pow(a_f_pre[k][i], 2.0), 2.0) + 4.0 * a_f_pre[k][i] * a_f_pre[k][i] * a_b_pre[k][i] * a_b_pre[k][i]) / std::pow(1.0 - std::pow(a_f_pre[k][i], 2.0), 3.0);
            }

            // Calculate sigma_b_pre_sqr.
            Vector3 sigma_b_pre;
            for (int i = 0; i < 3; i++) {
                sigma_b_pre[i] = (1.0 + 0.7 * a_f_pre[k][i] * a_f_pre[k][i]);
                sigma_b_pre[i] *= (a_b_pre[k][i] * std::sqrt(2.0 * beta_f_sqr_pre[k][i] + beta_b_sqr_pre[k][i]) + a_b_pre[k][i] * a_b_pre[k][i] * a_b_pre[k][i] * std::sqrt(2.0 * beta_f_sqr_pre[k][i] + 3.0 * beta_b_sqr_pre[k][i])) /
                        (a_b_pre[k][i] + a_b_pre[k][i] * a_b_pre[k][i] * a_b_pre[k][i] * (2.0 * std::sqrt(beta_f_sqr_pre[k][i]) + 3.0 * std::sqrt(beta_b_sqr_pre[k][i])));
                sigma_b_pre_sqr[k][i] = sigma_b_pre[i] * sigma_b_pre[i];
            }

            // Calculate N_X_G.
            N_R_G_pre[k] = Vector3(0.0, 0.0, 0.0);
            N_TT_G_pre[k] = Vector3(0.0, 0.0, 0.0);
            N_TRT_G_pre[k] = Vector3(0.0, 0.0, 0.0);
            N_TTs_G_pre[k] = Vector3(0.0, 0.0, 0.0);
            N_TRTs_G_pre[k] = Vector3(0.0, 0.0, 0.0);

            for (int j = 0; j < NUM_PHI_SAMPLES; j++) {
                Float phiPrime = (j + 0.5) / NUM_PHI_SAMPLES * M_PI / 2.0 + M_PI / 2.0;
                regularizePhi(phiPrime);
                Float dPhi = M_PI / 2.0 / NUM_PHI_SAMPLES;

                N_R_G_pre[k] += furSolver.Nu(0, phiPrime, theta_i) * dPhi;
                N_TT_G_pre[k] += furSolver.Nu(1, phiPrime, theta_i) * dPhi;
                N_TRT_G_pre[k] += furSolver.Nu(2, phiPrime, theta_i) * dPhi;

                N_TTs_G_pre[k] += max(discardLarger(furSolver.Ns(1, phiPrime, theta_i), 5.0f), 0.0f) * dPhi;
                N_TRTs_G_pre[k] += max(discardLarger(furSolver.Ns(2, phiPrime, theta_i), 5.0f), 0.0f) * dPhi;
            }
            N_R_G_pre[k] *= 2 / M_PI;
            N_TT_G_pre[k] *= 2 / M_PI;
            N_TRT_G_pre[k] *= 2 / M_PI;
            N_TTs_G_pre[k] *= 2 / M_PI;
            N_TRTs_G_pre[k] *= 2 / M_PI;
        }

        for (int channel = 0; channel < 3; channel++) {
            // Output alpha_b_pre.
            for (int k = 0; k < NUM_BINS_DS; k++) {
                alpha_b[k][batch * 3 + channel] = alpha_b_pre[k][channel];
                alpha_f[k][batch * 3 + channel] = alpha_f_pre[k][channel];
                beta_b_sqr[k][batch * 3 + channel] = beta_b_sqr_pre[k][channel];
                beta_f_sqr[k][batch * 3 + channel] = beta_f_sqr_pre[k][channel];
                a_b[k][batch * 3 + channel] = a_b_pre[k][channel];
                a_bu[k][batch * 3 + channel] = a_bu_pre[k][channel];
                a_bs[k][batch * 3 + channel] = a_bs_pre[k][channel];
                a_f[k][batch * 3 + channel] = a_f_pre[k][channel];
                a_fu[k][batch * 3 + channel] = a_fu_pre[k][channel];
                a_fs[k][batch * 3 + channel] = a_fs_pre[k][channel];
                A_b[k][batch * 3 + channel] = A_b_pre[k][channel];
                sigma_b_sqr[k][batch * 3 + channel] = sigma_b_pre_sqr[k][channel];

                N_R_G[k][batch * 3 + channel] = N_R_G_pre[k][channel];
                N_TT_G[k][batch * 3 + channel] = N_TT_G_pre[k][channel];
                N_TRT_G[k][batch * 3 + channel] = N_TRT_G_pre[k][channel];
                N_TTs_G[k][batch * 3 + channel] = N_TTs_G_pre[k][channel];
                N_TRTs_G[k][batch * 3 + channel] = N_TRTs_G_pre[k][channel];
            }
        }
    }
}

MTS_NAMESPACE_END

#endif
