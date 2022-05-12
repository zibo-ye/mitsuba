/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include "gausslegendre.hpp"
#include "precomputeazimuthallobe.hpp"

MTS_NAMESPACE_BEGIN

template <int n>
static Float Pow(Float v)
{
    static_assert(n > 0, "Power can’t be negative");
    Float n2 = Pow<n / 2>(v);
    return n2 * n2 * Pow<n & 1>(v);
}

template <>
Float Pow<1>(Float v) { return v; }

template <>
Float Pow<0>(Float v) { return 1; }

class Marschner : public BSDF
{
public:
    Marschner(const Properties &props)
        : BSDF(props)
    {
        m_eta = props.getFloat("eta", 1.55);
        m_sigma = props.getSpectrum("sigma", Spectrum(0.5f));

        // Float eu[3] = { 0.419f, 0.697f, 1.37f };
        // Float ph[3] = { 0.187f, 0.4f, 1.05f };
        // const Spectrum eumelaninSigmaA = Spectrum(eu);
        // const Spectrum pheomelaninSigmaA = Spectrum(ph);
        // Float melanin_ratio = props.getFloat("melanin_ratio", 1.0f);
        // Float melanin_concentration = props.getFloat("melanin_concentration", 8.0f);
        // m_sigma = melanin_concentration * (eumelaninSigmaA * (1.f - melanin_ratio) + pheomelaninSigmaA * melanin_ratio);

        m_alpha = props.getFloat("alpha", 2.5f);
        m_alpha = degToRad(m_alpha);
        m_roughness = props.getFloat("roughness", 0.2f);
        m_roughness = props.getFloat("roughness", 0.2f);

        m_betaR = std::max(M_PI_2 * m_roughness, 0.04);
        m_betaTT = m_betaR * 0.5f;
        m_betaTRT = m_betaR * 2.0f;

        m_disableR = props.getBoolean("disableR", false);
        m_disableTT = props.getBoolean("disableTT", false);
        m_disableTRT = props.getBoolean("disableTRT", false);

        precomputeAzimuthalDistributions();
        std::cout << toString() << endl;
    }

    Marschner(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager)
    {

        m_eta = stream->readFloat();
        m_sigma = Spectrum(stream);

        m_alpha = stream->readFloat();
        m_roughness = stream->readFloat();

        configure();
    }

    void configure()
    {

        m_components.clear();

        m_components.push_back(EDiffuseReflection | EFrontSide);
        m_usesRayDifferentials = false;

        precomputeAzimuthalDistributions();

        BSDF::configure();
    }

    Spectrum brdf(Float thetaI, Float thetaR, Float Phi) const
    {

        if (Phi < 0.0f)
            Phi += M_PI * 2.0f;
        if (Phi > M_PI * 2.0f)
            Phi -= M_PI * 2.0f;

        Float thetaD = std::abs(thetaR - thetaI) * 0.5f;
        Float cosThetaD = std::cos(thetaD);

        Float sinThetaO = std::sin(thetaR);
        Float cosThetaO = std::cos(thetaR);

        auto Mps = M(std::sin(thetaI), std::cos(thetaI), sinThetaO, cosThetaO);

        Spectrum result(0.0f);

        if (!m_disableR)
            result += Mps[0] * _nR->eval(Phi, cosThetaD);
        if (!m_disableTT)
            result += Mps[1] * _nTT->eval(Phi, cosThetaD);
        if (!m_disableTRT)
            result += Mps[2] * _nTRT->eval(Phi, cosThetaD);

        return result;
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure = ESolidAngle) const
    {
        if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle)
            return Spectrum(0.0f);

        // geoFrame defined as X: axis, Z: normal
        // bRec.wi and wo is defined at local space of shFrame
        auto wiworld = normalize(bRec.its.toWorld(bRec.its.wi));
        auto woworld = normalize(bRec.its.toWorld(bRec.wo));
        auto wi = bRec.its.geoFrame.toLocal(wiworld);
        auto wo = bRec.its.geoFrame.toLocal(woworld);

        Float thetaI = std::asin(math::clamp(wi.x, -1.0f, 1.0f));
        Float thetaR = std::asin(math::clamp(wo.x, -1.0f, 1.0f));
        Float phiI = std::atan2(wi.y, wi.z);
        Float phiR = std::atan2(wo.y, wo.z);
        Float phi = phiR - phiI;
        Spectrum result(0.0f);

        return brdf(thetaI, thetaR, phi);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const
    {
        if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle)
            return 0.0f;

        return warp::squareToUniformSpherePdf();
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const
    {
        if (!(bRec.typeMask & EDiffuseReflection))
            return Spectrum(0.0f);

        Float _pdf;
        return this->sample(bRec, _pdf, sample);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const
    {
        if (!(bRec.typeMask & EDiffuseReflection))
            return Spectrum(0.0f);

        bRec.wo = warp::squareToUniformSphere(sample);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;
        pdf = warp::squareToUniformSpherePdf();
        return eval(bRec);
    }

    void addChild(const std::string &name, ConfigurableObject *child)
    {
        BSDF::addChild(name, child);
    }

    void serialize(Stream *stream, InstanceManager *manager) const
    {
        BSDF::serialize(stream, manager);

        stream->writeFloat(m_eta);
        m_sigma.serialize(stream);

        // manager->serialize(stream, m_diffuseReflectance.get());
    }

    std::string toString() const
    {
        std::ostringstream oss;
        oss << "Marschner Hair[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  alpha = \"" << m_alpha << "\"," << endl
            << "  beta = \"" << m_betaR << "\"," << endl
            << "  eta = \"" << m_eta << "\"," << endl
            << "  sigma = \"" << m_sigma.toString() << "\"," << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()

private:
    float sampleM(float v, float sinThetaI, float cosThetaI, float xi1, float xi2) const
    {
        // Version from the paper (very unstable)
        // float cosTheta = v*std::log(std::exp(1.0f/v) - 2.0f*xi1*std::sinh(1.0f/v));
        // More stable version from "Numerically stable sampling of the von Mises Fisher distribution on S2 (and other tricks)"
        float cosTheta = 1.0f + v * std::log(xi1 + (1.0f - xi1) * std::exp(-2.0f / v));
        float sinTheta = math::trigInverse(cosTheta);
        float cosPhi = std::cos(M_PI * 2.0f * xi2);

        return -cosTheta * sinThetaI + sinTheta * cosPhi * cosThetaI;
    }

    Float Phi(float gammaI, float gammaT, int p)
    {
        return 2.0f * p * gammaT - 2.0f * gammaI + p * M_PI;
    }

    Float I0(Float x) const
    {
        Float result = 1.0f;
        Float xSq = x * x;
        Float xi = xSq;
        Float denom = 4.0f;
        for (int i = 1; i <= 10; ++i)
        {
            result += xi / denom;
            xi *= xSq;
            denom *= 4.0f * Float((i + 1) * (i + 1));
        }
        return result;
    }

    Float LogI0(Float x) const
    {
        if (x > 12.0f)
            // More stable evaluation of log(I0(x))
            // See also https://publons.com/discussion/12/
            return x + 0.5f * (std::log(1.0f / (M_PI * 2.0f * x)) + 1.0f / (8.0f * x));
        else
            return std::log(I0(x));
    }

    // calculate M
    std::vector<Float> M(Float sinThetaI, Float cosThetaI, Float sinThetaO, Float cosThetaO) const
    {
        Float alpha[3], sinThetaIp[3], cosThetaIp[3];
        alpha[0] = -2 * m_alpha;
        alpha[1] = m_alpha;
        alpha[2] = 4 * m_alpha;
        for (int i = 0; i < 3; ++i)
        {
            Float sinAlpha = std::sin(alpha[i]);
            Float cosAlpha = math::safe_sqrt(1 - math::sqr(alpha[i]));
            sinThetaIp[i] = sinThetaI * cosAlpha + cosThetaI * sinAlpha;
            cosThetaIp[i] = cosThetaI * cosAlpha - sinThetaI * sinAlpha;
        }
        Float v[3];
        v[0] = m_betaR * m_betaR;
        v[1] = .25 * v[0];
        v[2] = 4 * v[0];
        std::vector<Float> result(3);
        for (int i = 0; i < 3; ++i)
        {
            result[i] = Mp(cosThetaIp[i], cosThetaO, sinThetaIp[i], sinThetaO, v[i]);
        }
        return result;
    }

    Float Mp(Float cosThetaI, Float cosThetaO, Float sinThetaI, Float sinThetaO, Float v) const
    {
        Float a = cosThetaI * cosThetaO / v;
        Float b = sinThetaI * sinThetaO / v;
        Float mp = (v <= .1) ? (std::exp(LogI0(a) - b - 1 / v + 0.6931f + std::log(1 / (2 * v)))) : (std::exp(-b) * I0(a)) / (std::sinh(1 / v) * 2 * v);
        return mp;
    }

    // Standard normalized Gaussian
    float g(float beta, float theta)
    {
        return std::exp(-theta * theta / (2.0f * beta * beta)) / (std::sqrt(2.0f * M_PI) * beta);
    }

    float D(float beta, float phi)
    {
        float result = 0.0f;
        float delta;
        float shift = 0.0f;
        do
        {
            delta = g(beta, phi + shift) + g(beta, phi - shift - 2.0f * M_PI);
            result += delta;
            shift += 2.0f * M_PI;
        } while (delta > 1e-4f);
        return result;
    }

    static inline Float dielectricReflectance(Float eta, Float cosThetaI)
    {
        Float cosThetaT;
        if (cosThetaI < 0.0f)
        {
            eta = 1.0f / eta;
            cosThetaI = -cosThetaI;
        }
        Float sinThetaTSq = eta * eta * (1.0f - cosThetaI * cosThetaI);
        if (sinThetaTSq > 1.0f)
        {
            cosThetaT = 0.0f;
            return 1.0f;
        }
        cosThetaT = std::sqrt(std::max(1.0f - sinThetaTSq, 0.0f));

        Float Rs = (eta * cosThetaI - cosThetaT) / (eta * cosThetaI + cosThetaT);
        Float Rp = (eta * cosThetaT - cosThetaI) / (eta * cosThetaT + cosThetaI);

        return (Rs * Rs + Rp * Rp) * 0.5f;
    }

    // Precompute LUT
    void precomputeAzimuthalDistributions()
    {
        const int Resolution = PrecomputedAzimuthalLobe::AzimuthalResolution;
        std::unique_ptr<Spectrum[]> valuesR(new Spectrum[Resolution * Resolution]);
        std::unique_ptr<Spectrum[]> valuesTT(new Spectrum[Resolution * Resolution]);
        std::unique_ptr<Spectrum[]> valuesTRT(new Spectrum[Resolution * Resolution]);

// Ideally we could simply make this a constexpr, but MSVC does not support that yet (boo!)
#define NumPoints 140

        GaussLegendre<NumPoints> integrator;
        const auto points = integrator.points();
        const auto weights = integrator.weights();

        // Cache the gammaI across all integration points
        std::array<float, NumPoints> gammaIs;
        for (int i = 0; i < NumPoints; ++i)
            gammaIs[i] = std::asin(points[i]);

        // Precompute the Gaussian detector and sample it into three 1D tables.
        // This is the only part of the precomputation that is actually approximate.
        // 2048 samples are enough to support the lowest roughness that the BCSDF
        // can reliably simulate
        const int NumGaussianSamples = 2048;
        std::unique_ptr<float[]> Ds[3];
        for (int p = 0; p < 3; ++p)
        {
            Ds[p].reset(new float[NumGaussianSamples]);
            for (int i = 0; i < NumGaussianSamples; ++i)
                Ds[p][i] = D(m_betaR, i / (NumGaussianSamples - 1.0f) * M_PI * 2.0f);
        }

        // Simple wrapped linear interpolation of the precomputed table
        auto approxD = [&](int p, float phi)
        {
            float u = std::abs(phi * (0.5f * M_1_PI * (NumGaussianSamples - 1)));
            int x0 = int(u);
            int x1 = x0 + 1;
            u -= x0;
            return Ds[p][x0 % NumGaussianSamples] * (1.0f - u) + Ds[p][x1 % NumGaussianSamples] * u;
        };

        // Here follows the actual precomputation of the azimuthal scattering functions
        // The scattering functions are parametrized with the azimuthal angle phi,
        // and the cosine of the half angle, cos(thetaD).
        // This parametrization makes the azimuthal function relatively smooth and allows using
        // really low resolutions for the table (64x64 in this case) without any visual
        // deviation from ground truth, even at the lowest supported roughness setting
        for (int y = 0; y < Resolution; ++y)
        {
            float cosHalfAngle = y / (Resolution - 1.0f);

            // Precompute reflection Fresnel factor and reduced absorption coefficient
            float iorPrime = std::sqrt(m_eta * m_eta - (1.0f - cosHalfAngle * cosHalfAngle)) / cosHalfAngle;
            float cosThetaT = std::sqrt(1.0f - (1.0f - cosHalfAngle * cosHalfAngle) * math::sqr(1.0f / m_eta));
            Spectrum sigmaAPrime = m_sigma / cosThetaT;

            // Precompute gammaT, f_t and internal absorption across all integration points
            std::array<float, NumPoints> fresnelTerms, gammaTs;
            std::array<Spectrum, NumPoints> absorptions;
            for (int i = 0; i < NumPoints; ++i)
            {
                gammaTs[i] = std::asin(math::clamp(points[i] / iorPrime, -1.0f, 1.0f));
                fresnelTerms[i] = dielectricReflectance(1.0f / m_eta, cosHalfAngle * std::cos(gammaIs[i]));
                // std::cout<< gammaIs[i] << " " << m_eta << " " << fresnelTerms[i] << std::endl;
                absorptions[i] = (-sigmaAPrime * 2.0f * std::cos(gammaTs[i])).exp();
            }

            for (int phiI = 0; phiI < Resolution; ++phiI)
            {
                float phi = M_PI * 2.0f * phiI / (Resolution - 1.0f);

                float integralR = 0.0f;
                Spectrum integralTT(0.0f);
                Spectrum integralTRT(0.0f);

                // Here follows the integration across the fiber width, h.
                // Since we were able to precompute most of the factors that
                // are constant w.r.t phi for a given h,
                // we don't have to do much work here.
                for (int i = 0; i < integrator.numSamples(); ++i)
                {
                    float fR = fresnelTerms[i];
                    Spectrum T = absorptions[i];

                    float AR = fR;
                    Spectrum ATT = (1.0f - fR) * (1.0f - fR) * T;
                    Spectrum ATRT = ATT * fR * T;

                    integralR += weights[i] * approxD(0, phi - Phi(gammaIs[i], gammaTs[i], 0)) * AR;
                    integralTT += weights[i] * approxD(1, phi - Phi(gammaIs[i], gammaTs[i], 1)) * ATT;
                    integralTRT += weights[i] * approxD(2, phi - Phi(gammaIs[i], gammaTs[i], 2)) * ATRT;
                }
                valuesR[phiI + y * Resolution] = Spectrum(0.5f * integralR);
                valuesTT[phiI + y * Resolution] = 0.5f * integralTT;
                valuesTRT[phiI + y * Resolution] = 0.5f * integralTRT;
            }
        }

        // Hand the values off to the helper class to construct sampling CDFs and so forth
        _nR = new PrecomputedAzimuthalLobe(std::move(valuesR));
        _nTT = new PrecomputedAzimuthalLobe(std::move(valuesTT));
        _nTRT = new PrecomputedAzimuthalLobe(std::move(valuesTRT));
    }

    Float m_eta;      // index of refraction
    Spectrum m_sigma; // absorption coefficient
    Float m_alpha;    // the angle that the small scales on the surface of hair are offset from the base cylinder, expressed in degree | longitudinal shift: R lobe (typ 5◦ to 10◦)
    Float m_beta;     // longitudinal width (stdev.): R lobe (typ 5◦ to 10◦)
    Float m_roughness;

    bool m_disableR;
    bool m_disableTT;
    bool m_disableTRT;

    Float m_betaR, m_betaTT, m_betaTRT;
    Float _alpha[3];
    Float _beta[3];

    PrecomputedAzimuthalLobe *_nR;
    PrecomputedAzimuthalLobe *_nTT;
    PrecomputedAzimuthalLobe *_nTRT;

    // ref<Texture> m_diffuseReflectance;
};

MTS_IMPLEMENT_CLASS_S(Marschner, false, BSDF)
MTS_EXPORT_PLUGIN(Marschner, "Marschner's Hair Model")
MTS_NAMESPACE_END
