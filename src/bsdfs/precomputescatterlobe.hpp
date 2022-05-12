#ifndef PRECOMPUTESCATTER_HPP_
#define PRECOMPUTESCATTER_HPP_

#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

#define NUM_SCATTERING_INNER 64
#define NUM_H 64
#define NUM_THETA 64
#define NUM_G 16
#define NUM_BINS 720

class PrecomputedScatterLobe
{

private:
    float scattered[NUM_SCATTERING_INNER][NUM_H][NUM_G][NUM_BINS];

public:
    PrecomputedScatterLobe(std::string path)
    {
        FILE *fMedulla = fopen((std::string(path)).c_str(), "rb");
        fread(scattered, sizeof(float), NUM_SCATTERING_INNER * NUM_H * NUM_G * NUM_BINS, fMedulla);
        fclose(fMedulla);
    }

    Float evalM(Float sigma, Float thetaI, Float g, Float thetaR) const
    {
        // thetaI normalize from [-PI/2,PI/2]
        Float thetaI_normalized = thetaI / M_PI + 0.5;
        // thetaI normalize from [0,2PI]
        Float thetaR_normalized = thetaR / (2 * M_PI);
        // Float thetaR_normalized = thetaR / M_PI + 0.5;
        return eval(sigma, math::clamp(thetaI_normalized, 0.0f, 1.0f), g, math::clamp(thetaR_normalized, 0.0f, 1.0f));
    }

    Float evalN(Float sigma, Float h, Float g, Float phi) const
    {
        // h normalize between [-1,1]
        Float h_normalized = h * 0.5 + 0.5;
        phi = std::abs(phi);
        phi = phi > M_PI ? 2 * M_PI - phi : phi;
        Float phi_normalized = phi / (M_PI);
        return eval(sigma, math::clamp(h_normalized, 0.0f, 1.0f), g, math::clamp(phi_normalized, 0.0f, 1.0f));
    }

    // private:
    Float eval(Float sigma, Float thetaI_or_h_normlized, Float g, Float thetaR_or_phi_normlized) const
    {
        // sigma normalized from [0,20]
        int scatter_index = math::clamp(sigma / 20.0f, 0.0f, 0.9999f) * NUM_SCATTERING_INNER;
        // g normalize from [0, 0.8]
        int g_index = math::clamp(g / 0.8f, 0.0f, 0.9999f) * NUM_G;
        // phi normalized from [0, 2PI]
        int bin_index = thetaR_or_phi_normlized * 720;
        int theta_index = thetaI_or_h_normlized * 64;
        auto val = scattered[scatter_index][theta_index][g_index][bin_index];
        return math::clamp(val, 0.0f, 1000.0f);
    }
};
MTS_NAMESPACE_END
#endif