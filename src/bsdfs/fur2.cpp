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
#include "fur2_basics.h"
#include "mitsuba/render/scene.h"
#include <vector>
#include <map>
using namespace std;

MTS_NAMESPACE_BEGIN

class Fur2 : public BSDF {
public:
    Fur2(const Properties &props)
        : BSDF(props) {
        // If color texture is given, read and convert it into absorption on the fly
        colorTexture = new ConstantSpectrumTexture(props.getSpectrum("color", Spectrum(1.0f)));

        useColorTexture = props.getBoolean("useColorTexture", false);
        medulla_ratio = props.getFloat("medulla_ratio", 0.0);
        eta = props.getFloat("eta", 1.55);
        beta_m = props.getFloat("beta_m", 0.055);
        beta_n = props.getFloat("beta_n", 0.055);
        alpha = props.getFloat("alpha", 3.0);
        g = props.getFloat("g", 0.0);
        layers = props.getFloat("layers", 0.5);
        absorptionFactor = props.getFloat("absorptionFactor", 4.0);
        absorption_outer = props.getVector("absorption_outer", Vector3(0.0f, 0.0f, 0.0f));
        absorption_inner = props.getFloat("absorption_inner", 0.0f);
        scattering_inner = props.getFloat("scattering_inner", 0.5f);
        furMaskCode = props.getInteger("furMaskCode", 31);

        // Render mode:
        //     "near":      near field only
        //     "far":       far field only
        //     "nearfar":   near field for first bounce, far field for others
        renderMode = props.getString("renderMode", "far");
        if (renderMode != "near" && renderMode != "far" && renderMode != "nearfar" && renderMode != "midfar")
            printf("Invalid input for render mode!\n");


        // Color scale:
        //     > 0:             scale by colorScale
        // Color gamma:
        //     > 0:             gamma by colorScale
        colorScale = props.getFloat("colorScale", 2.0);
        colorGamma = props.getFloat("colorGamma", 1.0);

        // Tilt should always be negative for physical correctness.
        alpha = -degToRad(alpha);

        initialize("medulla");
    }

    Fur2(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        colorTexture = static_cast<Texture *>(manager->getInstance(stream));

//        perspectiveCamera = static_cast<PerspectiveCamera *>(manager->getInstance(stream));

        useColorTexture = stream->readBool();
        medulla_ratio = stream->readFloat();
        eta = stream->readFloat();
        beta_m = stream->readFloat();
        beta_n = stream->readFloat();
        alpha = stream->readFloat();
        g = stream->readFloat();
        layers = stream->readFloat();
        absorptionFactor = stream->readFloat();
        absorption_outer.x = stream->readFloat();
        absorption_outer.y = stream->readFloat();
        absorption_outer.z = stream->readFloat();
        absorption_inner = stream->readFloat();
        scattering_inner = stream->readFloat();

        renderMode = stream->readString();
        colorScale = stream->readFloat();
        colorGamma = stream->readFloat();

        initialize("medulla");

        configure();
    }

    void configure() {
        /* Verify the input parameter and fix them if necessary */
        colorTexture = ensureEnergyConservation(colorTexture, "color", 1.0f);

        m_components.clear();
        //m_components.push_back(EFUR2);
        if (colorTexture->getMaximum().max() > 0)
            m_components.push_back(EDiffuseReflection | EFrontSide
                | (colorTexture->isConstant() ? 0 : ESpatiallyVarying));
            m_usesRayDifferentials = colorTexture->usesRayDifferentials();

        BSDF::configure();
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        return colorTexture->eval(its);
    }

    Vector3 normalizeSafe(Vector3 v) const {
        if (v.length() == 0.0)
            return v;
        else
            return normalize(v);
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        Vector3 absorption = absorption_outer;
        if (useColorTexture && medulla_ratio < 1) {
            Spectrum color = colorTexture->eval(bRec.its, false);
            Float color_r, color_g, color_b;
            color.toLinearRGB(color_r, color_g, color_b);

            color_r = clamp(std::pow(color_r, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
            color_g = clamp(std::pow(color_g, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
            color_b = clamp(std::pow(color_b, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));

            Float medulla_occlusion = std::pow(medulla_ratio, 1.5);
            Float absorption_r = (-log(color_r) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
            Float absorption_g = (-log(color_g) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
            Float absorption_b = (-log(color_b) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);

            absorption_r = clamp(absorption_r, 0.01, MAX_ABSORPTION_OUTER);
            absorption_g = clamp(absorption_g, 0.01, MAX_ABSORPTION_OUTER);
            absorption_b = clamp(absorption_b, 0.01, MAX_ABSORPTION_OUTER);

            absorption = Vector3(absorption_r, absorption_g, absorption_b);
        }

        Vector3 wi = bRec.its.toWorld(bRec.wi);
        Vector3 wo = bRec.its.toWorld(bRec.wo);
        const Vector3 &u = bRec.its.geoFrame.s;
        const Vector3 &v = bRec.its.geoFrame.t;
        const Vector3 &w = bRec.its.geoFrame.n;

        Float theta_i = M_PI / 2.0f - acosSafe(dot(u, wi));
        Float theta_r = M_PI / 2.0f - acosSafe(dot(u, wo));
        theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));
        theta_r = clamp(theta_r, degToRad(-89.5), degToRad(89.5));

        Float phi_i = acosSafe(dot(normalizeSafe(wi - dot(u, wi) * u), v)) * (dot(w, wi) > 0 ? 1 : -1);
        Float phi_r = acosSafe(dot(normalizeSafe(wo - dot(u, wo) * u), v)) * (dot(w, wo) > 0 ? 1 : -1);
        Float phi = phi_r - phi_i;
        regularizePhi(phi);

        FurSolver furSolver(eta, absorption, alpha, beta_m, beta_n,
                            medulla_ratio, absorption_inner, scattering_inner, g, layers, bRec.sampler);

        Vector3 R, TT, TRT, TTs, TRTs;

        Float ms = furSolver.Ms(theta_i - alpha, theta_r - alpha, phi);

        if (renderMode == "near" || (renderMode == "nearfar" && bRec.depth == 1)) {
            Float h = bRec.its.storage[1];
            R = furSolver.Mu(0, theta_i, theta_r, phi) * furSolver.Nu(0, h, phi, theta_i);
            TT = furSolver.Mu(1, theta_i, theta_r, phi) * furSolver.Nu(1, h, phi, theta_i);
            TRT = furSolver.Mu(2, theta_i, theta_r, phi) * furSolver.Nu(2, h, phi, theta_i);
            TTs = ms * furSolver.Ns(1, h, phi, theta_i);
            TRTs = ms * furSolver.Ns(2, h, phi, theta_i);
        } else if (renderMode == "midfar" && bRec.depth == 1) {
            Point origin = perspectiveCamera->s_origin;
            Vector3 lookAt = perspectiveCamera->s_lookAt;

            Point p = bRec.its.p;
            Vector3 pVec = p - origin;
            Float projDist = dot(pVec, lookAt);
            Float fovX = perspectiveCamera->getXFov() / 180.0 * M_PI;
            Float covX = tan(fovX / 2.0) * projDist * 2.0 * PIXEL_COVERAGE_SCALE;
            int resX = perspectiveCamera->getFilm()->getSize().x;
            Float pixCov = covX / resX;
            Float furRadius = bRec.its.storage[2];
            Vector3 pVecNorm = normalize(pVec);
            Float pixCovProj = pixCov * dot(pVecNorm, lookAt);
            Float hRange = pixCovProj / furRadius;

            Float h = bRec.its.storage[1];
            Float hMin = clamp(h - hRange / 2.0, -1.0 + EPS, 1.0 - EPS);
            Float hMax = clamp(h + hRange / 2.0, -1.0 + EPS, 1.0 - EPS);

//            if (hMin != -0.999 || hMax != 0.999)
//                printf("%f %f\n", hMin, hMax);

            if (hMin > 0.0) {
                R = furSolver.Mu(0, theta_i, theta_r, phi) * furSolver.NuMid(0, hMin, hMax, phi, theta_i);
                TT = furSolver.Mu(1, theta_i, theta_r, phi) * furSolver.NuMid(1, hMin, hMax, phi, theta_i);
                TRT = furSolver.Mu(2, theta_i, theta_r, phi) * furSolver.NuMid(2, hMin, hMax, phi, theta_i);
                TTs = ms * max(discardLarger(furSolver.NsMid(1, hMin, hMax, phi, theta_i), 10.0f), 0.0f);
                TRTs = ms * max(discardLarger(furSolver.NsMid(2, hMin, hMax, phi, theta_i), 10.0f), 0.0f);
            } else if (hMax < 0.0) {
                R = furSolver.Mu(0, theta_i, theta_r, phi) * furSolver.NuMid(0, -hMax, -hMin, -phi, theta_i);
                TT = furSolver.Mu(1, theta_i, theta_r, phi) * furSolver.NuMid(1, -hMax, -hMin, -phi, theta_i);
                TRT = furSolver.Mu(2, theta_i, theta_r, phi) * furSolver.NuMid(2, -hMax, -hMin, -phi, theta_i);
                TTs = ms * max(discardLarger(furSolver.NsMid(1, -hMax, -hMin, -phi, theta_i), 10.0f), 0.0f);
                TRTs = ms * max(discardLarger(furSolver.NsMid(2, -hMax, -hMin, -phi, theta_i), 10.0f), 0.0f);
            } else {
                R = furSolver.Mu(0, theta_i, theta_r, phi) * (furSolver.NuMid(0, 0, hMax, phi, theta_i) * hMax / (hMax - hMin) + furSolver.NuMid(0, 0, -hMin, -phi, theta_i) * -hMin / (hMax - hMin));
                TT = furSolver.Mu(1, theta_i, theta_r, phi) * (furSolver.NuMid(1, 0, hMax, phi, theta_i) * hMax / (hMax - hMin) + furSolver.NuMid(1, 0, -hMin, -phi, theta_i) * -hMin / (hMax - hMin));
                TRT = furSolver.Mu(2, theta_i, theta_r, phi) * (furSolver.NuMid(2, 0, hMax, phi, theta_i) * hMax / (hMax - hMin) + furSolver.NuMid(2, 0, -hMin, -phi, theta_i) * -hMin / (hMax - hMin));
                TTs = ms * (max(discardLarger(furSolver.NsMid(1, 0, hMax, phi, theta_i) * hMax / (hMax - hMin), 10.0f), 0.0f) + max(discardLarger(furSolver.NsMid(1, 0, -hMin, -phi, theta_i) * -hMin / (hMax - hMin), 10.0f), 0.0f));
                TRTs = ms * (max(discardLarger(furSolver.NsMid(2, 0, hMax, phi, theta_i) * hMax / (hMax - hMin), 10.0f), 0.0f) + max(discardLarger(furSolver.NsMid(2, 0, -hMin, -phi, theta_i) * -hMin / (hMax - hMin), 10.0f), 0.0f));
            }
////printf("%f %f %f\n", R.x, R.y, R.z);
//            TT = Vector3(0.0, 0.0, 0.0);
//            TRT = Vector3(0.0, 0.0, 0.0);
//            TTs = Vector3(0.0, 0.0, 0.0);
//            TRTs = Vector3(0.0, 0.0, 0.0);
        } else {
            R = furSolver.Mu(0, theta_i, theta_r, phi) * furSolver.Nu(0, phi, theta_i);
            TT = furSolver.Mu(1, theta_i, theta_r, phi) * furSolver.Nu(1, phi, theta_i);
            TRT = furSolver.Mu(2, theta_i, theta_r, phi) * furSolver.Nu(2, phi, theta_i);
            TTs = ms * max(discardLarger(furSolver.Ns(1, phi, theta_i), 10.0f), 0.0f);
            TRTs = ms * max(discardLarger(furSolver.Ns(2, phi, theta_i), 10.0f), 0.0f);

//            if (max(discardLarger(furSolver.Ns(1, phi, theta_i), 10.0f), 0.0f).length() > 1.0 ||
//                max(discardLarger(furSolver.Ns(2, phi, theta_i), 10.0f), 0.0f).length() > 1.0)
//            printf("R = %f %f %f\nTT=%f %f %f\nTRT=%f %f %f\nTTs=%f %f %f\nTRTs=%f %f %f\n\n",
//                   R.x, R.y, R.z,
//                   TT.x, TT.y, TT.z,
//                   TRT.x, TRT.y, TRT.z,
//                   TTs.x, TTs.y, TTs.z,
//                   TRTs.x, TRTs.y, TRTs.z);
        }
        R = furMaskCode & 1 ? R : Vector3(Float(0));
        TT = (furMaskCode >> 1) & 1 ? TT : Vector3(Float(0));
        TRT = (furMaskCode >> 2) & 1 ? TRT : Vector3(Float(0));
        TTs = (furMaskCode >> 3) & 1 ? TTs : Vector3(Float(0));
        TRTs = (furMaskCode >> 4) & 1 ? TRTs : Vector3(Float(0));

        Vector3 bsdfCosine;
        if (bRec.component == -1)
            bsdfCosine = (R + TT + TRT + TTs + TRTs) / cos(theta_r);
        else if (bRec.component == 0)
            bsdfCosine = (R + TT + TRT) / cos(theta_r);
        else
            bsdfCosine = (TTs + TRTs) / cos(theta_r);

        if (std::isnan(bsdfCosine[0]) || std::isnan(bsdfCosine[1]) || std::isnan(bsdfCosine[2]))
            return Spectrum(0.0f);

        if (std::isinf(bsdfCosine[0]) || std::isinf(bsdfCosine[1]) || std::isinf(bsdfCosine[2]))
            return Spectrum(0.0f);

        for (int i = 0; i < 3; i++)
            bsdfCosine[i] = std::min(std::max(Float(0.0f), bsdfCosine[i]), Float(1e3f));

        Spectrum retVal(0.0f);
        retVal.fromLinearRGB(bsdfCosine[0], bsdfCosine[1], bsdfCosine[2]);

        //Vector3 ret = Vector3(1, 1, 1)*1e5; retVal.fromLinearRGB(ret[0], ret[1], ret[2]); // for debug
        
        
        return retVal;
    }
    
    Float normrnd(const BSDFSamplingRecord &bRec, Float mu, Float sigma) const {
        double U = bRec.sampler->next1D();
        double V = bRec.sampler->next1D();
        double X = sqrt(-2.0 * log(U)) * cos(2.0 * M_PI * V);
        return X * sigma + mu;
    }

//    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
//        Vector3 wi = bRec.its.toWorld(bRec.wi);
//        Vector3 wo = bRec.its.toWorld(bRec.wo);
//        const Vector3 &u = bRec.its.geoFrame.s;
//        const Vector3 &v = bRec.its.geoFrame.t;
//        const Vector3 &w = bRec.its.geoFrame.n;

//        Float theta_o = M_PI / 2.0f - acosSafe(dot(u, wo));
//        theta_o = clamp(theta_o, degToRad(-89.5), degToRad(89.5));
//        Float theta_pdf = 0.5 * cos(theta_o);

//        Float phi_i = acosSafe(dot(normalizeSafe(wi - dot(u, wi) * u), v)) * (dot(w, wi) > 0 ? 1 : -1);
//        Float phi_o = 2.0f * asinSafe(bRec.sampler->next1D() * 2.0 - 1.0) + phi_i;
//        Float phi_pdf = cos((phi_o - phi_i) / 2.0f) / 4.0f;

//        return theta_pdf * phi_pdf / cos(theta_o);
//    }

//    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
//        Vector3 wi = bRec.its.toWorld(bRec.wi);
//        const Vector3 &u = bRec.its.geoFrame.s;
//        const Vector3 &v = bRec.its.geoFrame.t;
//        const Vector3 &w = bRec.its.geoFrame.n;

//        Float theta_i = M_PI / 2.0f - acosSafe(dot(u, wi));
//        theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));

//        Float phi_i = acosSafe(dot(normalizeSafe(wi - dot(u, wi) * u), v)) * (dot(w, wi) > 0 ? 1 : -1);


//        Float theta_o = asinSafe(bRec.sampler->next1D() * 2.0 - 1.0);
//        Float theta_pdf = 0.5 * cos(theta_o);

//        Float phi_o = 2.0f * asinSafe(bRec.sampler->next1D() * 2.0 - 1.0) + phi_i;
//        Float phi_pdf = cos((phi_o - phi_i) / 2.0f) / 4.0f;


//        Float ul = sin(theta_o);
//        Float vl = cos(theta_o) * cos(phi_o);
//        Float wl = cos(theta_o) * sin(phi_o);
//        bRec.wo = bRec.its.toLocal(ul * u + vl * v + wl * w);

//        bRec.sampledComponent = 0;
//        bRec.sampledType = EDiffuseReflection;
//        pdf = Fur2::pdf(bRec, ESolidAngle);
//        return eval(bRec, ESolidAngle) / pdf;
//    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        Vector3 absorption = absorption_outer;
        if (useColorTexture && medulla_ratio < 1) {
            Spectrum color = colorTexture->eval(bRec.its, false);
            Float color_r, color_g, color_b;
            color.toLinearRGB(color_r, color_g, color_b);

            color_r = clamp(std::pow(color_r, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
            color_g = clamp(std::pow(color_g, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
            color_b = clamp(std::pow(color_b, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));

            Float medulla_occlusion = std::pow(medulla_ratio, 1.5);
            Float absorption_r = (-log(color_r) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
            Float absorption_g = (-log(color_g) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
            Float absorption_b = (-log(color_b) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);

            absorption_r = clamp(absorption_r, 0.01, MAX_ABSORPTION_OUTER);
            absorption_g = clamp(absorption_g, 0.01, MAX_ABSORPTION_OUTER);
            absorption_b = clamp(absorption_b, 0.01, MAX_ABSORPTION_OUTER);

            absorption = Vector3(absorption_r, absorption_g, absorption_b);
        }

        Vector3 wi = bRec.its.toWorld(bRec.wi);
        Vector3 wo = bRec.its.toWorld(bRec.wo);
        const Vector3 &u = bRec.its.geoFrame.s;
        const Vector3 &v = bRec.its.geoFrame.t;
        const Vector3 &w = bRec.its.geoFrame.n;

        Float theta_i = M_PI / 2.0f - acosSafe(dot(u, wi));
        Float theta_r = M_PI / 2.0f - acosSafe(dot(u, wo));
        theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));
        theta_r = clamp(theta_r, degToRad(-89.5), degToRad(89.5));

        Float phi_i = acosSafe(dot(normalizeSafe(wi - dot(u, wi) * u), v)) * (dot(w, wi) > 0 ? 1 : -1);
        Float phi_r = acosSafe(dot(normalizeSafe(wo - dot(u, wo) * u), v)) * (dot(w, wo) > 0 ? 1 : -1);
        Float phi = phi_r - phi_i;
        regularizePhi(phi);

        FurSolver furSolver(eta, absorption, alpha, beta_m, beta_n,
                            medulla_ratio, absorption_inner, scattering_inner, g, layers, bRec.sampler);

        Float h;
        if (renderMode == "near" || (renderMode == "nearfar" && bRec.depth == 1)) {
            h = bRec.its.storage[1];
        } else if (renderMode == "midfar") {
            Point origin = perspectiveCamera->s_origin;
            Vector3 lookAt = perspectiveCamera->s_lookAt;
            Point p = bRec.its.p;
            Vector3 pVec = p - origin;
            Float projDist = dot(pVec, lookAt);
            Float fovX = perspectiveCamera->getXFov() / 180.0 * M_PI;
            Float covX = tan(fovX / 2.0) * projDist * 2.0 * PIXEL_COVERAGE_SCALE;
            int resX = perspectiveCamera->getFilm()->getSize().x;
            Float pixCov = covX / resX;
            Float furRadius = bRec.its.storage[2];
            Vector3 pVecNorm = normalize(pVec);
            Float pixCovProj = pixCov * dot(pVecNorm, lookAt);
            Float hRange = pixCovProj / furRadius;

            Float hMid = bRec.its.storage[1];
            Float hMin = clamp(hMid - hRange / 2.0, -1.0 + EPS, 1.0 - EPS);
            Float hMax = clamp(hMid + hRange / 2.0, -1.0 + EPS, 1.0 - EPS);

            h = bRec.randNum * (hMax - hMin) + hMin;
        } else {
            h = (bRec.randNum * 2.0 - 1.0) * (1.0 - EPS);
        }

        return std::max(Float(0.0f), furSolver.pdfMN(h, theta_i, theta_r, phi));
    }

    // Choose a lobe and sample it.
    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        //bRec.eta = Float(1.f);
        Vector3 absorption = absorption_outer;
        if (useColorTexture && medulla_ratio < 1) {
            Spectrum color = colorTexture->eval(bRec.its, false);
            Float color_r, color_g, color_b;
            color.toLinearRGB(color_r, color_g, color_b);

            color_r = clamp(std::pow(color_r, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
            color_g = clamp(std::pow(color_g, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
            color_b = clamp(std::pow(color_b, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));

            Float medulla_occlusion = std::pow(medulla_ratio, 1.5);
            Float absorption_r = (-log(color_r) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
            Float absorption_g = (-log(color_g) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
            Float absorption_b = (-log(color_b) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);

            absorption_r = clamp(absorption_r, 0.01, MAX_ABSORPTION_OUTER);
            absorption_g = clamp(absorption_g, 0.01, MAX_ABSORPTION_OUTER);
            absorption_b = clamp(absorption_b, 0.01, MAX_ABSORPTION_OUTER);

            absorption = Vector3(absorption_r, absorption_g, absorption_b);
        }

        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;

        Vector3 wi = bRec.its.toWorld(bRec.wi);
        const Vector3 &u = bRec.its.geoFrame.s;
        const Vector3 &v = bRec.its.geoFrame.t;
        const Vector3 &w = bRec.its.geoFrame.n;

        Float theta_i = M_PI / 2.0f - acosSafe(dot(u, wi));
        theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));

        Float phi_i = acosSafe(dot(normalizeSafe(wi - dot(u, wi) * u), v)) * (dot(w, wi) > 0 ? 1 : -1);

        FurSolver furSolver(eta, absorption, alpha, beta_m, beta_n,
                            medulla_ratio, absorption_inner, scattering_inner, g, layers, bRec.sampler);

        // Choose a lobe to sample.
        Float h;
        if (renderMode == "near" || (renderMode == "nearfar" && bRec.depth == 1)) {
            h = bRec.its.storage[1];
        } else if (renderMode == "midfar") {
            Point origin = perspectiveCamera->s_origin;
            Vector3 lookAt = perspectiveCamera->s_lookAt;
            Point p = bRec.its.p;
            Vector3 pVec = p - origin;
            Float projDist = dot(pVec, lookAt);
            Float fovX = perspectiveCamera->getXFov() / 180.0 * M_PI;
            Float covX = tan(fovX / 2.0) * projDist * 2.0 * PIXEL_COVERAGE_SCALE;
            int resX = perspectiveCamera->getFilm()->getSize().x;
            Float pixCov = covX / resX;
            Float furRadius = bRec.its.storage[2];
            Vector3 pVecNorm = normalize(pVec);
            Float pixCovProj = pixCov * dot(pVecNorm, lookAt);
            Float hRange = pixCovProj / furRadius;

            Float hMid = bRec.its.storage[1];
            Float hMin = clamp(hMid - hRange / 2.0, -1.0 + EPS, 1.0 - EPS);
            Float hMax = clamp(hMid + hRange / 2.0, -1.0 + EPS, 1.0 - EPS);

            h = bRec.randNum * (hMax - hMin) + hMin;
        } else {
            h = (bRec.randNum * 2.0 - 1.0) * (1.0 - EPS);
        }

        int lobe;
        Float pdfLobe;
        furSolver.chooseLobe(theta_i, h, lobe, pdfLobe);
        
        bRec.sampledComponent = lobe < 3 ? 0 : 1;
        
        if (lobe == -1) {
            pdf = 0.0f;
            return Spectrum(0.0f);
        }
        
        if (((furMaskCode >> lobe) & 1) == 0) {
            pdf = 0.0f;
            return Spectrum(0.0f);
        }

        // Sample the chosen lobe.
        Float theta_o, phi_o;
        Float pdfM, pdfN;
        Float weightM;
        Vector3 weightN;
        furSolver.sampleLobe(lobe, h, theta_i, phi_i, theta_o, pdfM, weightM, phi_o, pdfN, weightN);

        regularizePhi(phi_o);
        regularizeTheta(theta_o);
        theta_o = clamp(theta_o, degToRad(-89.5), degToRad(89.5));

        // Convert angles to outgoing direction.
        Float ul = sin(theta_o);
        Float vl = cos(theta_o) * cos(phi_o);
        Float wl = cos(theta_o) * sin(phi_o);
        bRec.wo = bRec.its.toLocal(u * ul + v * vl + w * wl);

        // Calculate the final sampling weight and pdf.
        Vector3 weight = weightM * weightN / pdfLobe;

        if (std::isnan(weight[0]) || std::isnan(weight[1]) || std::isnan(weight[2]))
            return Spectrum(0.0f);

        if (std::isinf(weight[0]) || std::isinf(weight[1]) || std::isinf(weight[2]))
            return Spectrum(0.0f);
        

        weight[0] = std::max(Float(0.0f), std::min(weight[0], Float(2.0f)));
        weight[1] = std::max(Float(0.0f), std::min(weight[1], Float(2.0f)));
        weight[2] = std::max(Float(0.0f), std::min(weight[2], Float(2.0f)));

        // Compute the pdf. This can be more efficient.
        Float phi = phi_o - phi_i;
        pdf = pdfM * pdfN * pdfLobe / cos(theta_o);

        Spectrum retVal(0.f);
        retVal.fromLinearRGB(weight[0], weight[1], weight[2]);
        return retVal;

//        // Calculate the final sampling weight and pdf.
//        Float phi = phi_o - phi_i;
////        pdf = pdfM * pdfN * pdfLobe / cos(theta_o);
//        pdf = furSolver.pdfMN(h, theta_i, theta_o, phi);

//        Spectrum retVal = eval(bRec, ESolidAngle) / pdf;

//        Vector3 weight;
//        retVal.toLinearRGB(weight[0], weight[1], weight[2]);

//        if (std::isnan(weight[0]) || std::isnan(weight[1]) || std::isnan(weight[2]))
//            return Spectrum(0.0f);

//        if (std::isinf(weight[0]) || std::isinf(weight[1]) || std::isinf(weight[2]))
//            return Spectrum(0.0f);

//        weight[0] = std::max(Float(0.0f), std::min(weight[0], Float(2.0f)));
//        weight[1] = std::max(Float(0.0f), std::min(weight[1], Float(2.0f)));
//        weight[2] = std::max(Float(0.0f), std::min(weight[2], Float(2.0f)));

//        retVal.fromLinearRGB(weight[0], weight[1], weight[2]);
//        return retVal;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        Float pdf;
        return Fur2::sample(bRec, pdf, sample);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))
                && name == "color") {
            colorTexture = static_cast<Texture *>(child);
        } else if (child->getClass()->derivesFrom(MTS_CLASS(Sensor))) {
//            printf("Yes\n");
            perspectiveCamera = static_cast<PerspectiveCamera *>(child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        manager->serialize(stream, colorTexture.get());
//        perspectiveCamera->serialize(stream, manager);

        stream->writeBool(useColorTexture);
        stream->writeFloat(medulla_ratio);
        stream->writeFloat(eta);
        stream->writeFloat(beta_m);
        stream->writeFloat(beta_n);
        stream->writeFloat(alpha);
        stream->writeFloat(g);
        stream->writeFloat(layers);
        stream->writeFloat(absorptionFactor);
        stream->writeFloat(absorption_outer.x);
        stream->writeFloat(absorption_outer.y);
        stream->writeFloat(absorption_outer.z);
        stream->writeFloat(absorption_inner);
        stream->writeFloat(scattering_inner);

        stream->writeString(renderMode);
        stream->writeFloat(colorScale);
        stream->writeFloat(colorGamma);
    }

    Float getRoughness(const Intersection &its, int component) const {
        return std::numeric_limits<Float>::infinity();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "SmoothDiffuse[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  reflectance = " << indent(colorTexture->toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    std::string renderMode;
    Float colorScale, colorGamma;
    ref<Texture> colorTexture;
    bool useColorTexture;
    Float medulla_ratio;
    Float eta;
    Vector3 absorption_outer;
    Float beta_m, beta_n;
    Float scattering_inner, absorption_inner;
    Float alpha;
    Float g;
    Float layers;
    Float absorptionFactor;
    int furMaskCode;
    PerspectiveCamera *perspectiveCamera;
};

// ================ Hardware shader implementation ================

class Fur2Shader : public Shader {
public:
    Fur2Shader(Renderer *renderer, const Texture *reflectance)
        : Shader(renderer, EBSDFShader), colorTexture(reflectance) {
        colorTextureShader = renderer->registerShaderForResource(colorTexture.get());
    }

    bool isComplete() const {
        return colorTextureShader.get() != NULL;
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(colorTexture.get());
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(colorTextureShader.get());
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
            << "        return vec3(0.0);" << endl
            << "    return " << depNames[0] << "(uv) * inv_pi * cosTheta(wo);" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return " << evalName << "(uv, wi, wo);" << endl
            << "}" << endl;
    }

    MTS_DECLARE_CLASS()
private:
    ref<const Texture> colorTexture;
    ref<Shader> colorTextureShader;
};

Shader *Fur2::createShader(Renderer *renderer) const {
    return new Fur2Shader(renderer, colorTexture.get());
}

MTS_IMPLEMENT_CLASS(Fur2Shader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Fur2, false, BSDF)
MTS_EXPORT_PLUGIN(Fur2, "Our fur2 BSDF")
MTS_NAMESPACE_END
