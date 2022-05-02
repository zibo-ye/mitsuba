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
#include "precomputescatterlobe.hpp"

MTS_NAMESPACE_BEGIN

template <int n>static Float Pow(Float v) {
    static_assert(n > 0, "Power can’t be negative");
    Float n2 = Pow<n / 2>(v);
    return n2 * n2 * Pow<n & 1>(v);
}
        
template <> Float Pow<1>(Float v) { return v; }

template <> Float Pow<0>(Float v) { return 1; }

class YanFurRT : public BSDF {
public:
    YanFurRT(const Properties &props)
        : BSDF(props) {
        m_kappa  = props.getFloat("kappa", 0.88f);
		m_eta = props.getFloat("eta", 1.55);

        m_alpha = props.getFloat("alpha", 5.48f);
        m_alpha = degToRad(m_alpha);
        m_beta_m = props.getFloat("beta_m", 11.64f);
        m_beta_n = props.getFloat("beta_n", 7.49f);
        m_beta_m = degToRad(m_beta_m);
        m_beta_n = degToRad(m_beta_n);

        m_sigma_ca = props.getSpectrum("sigma_ca", Spectrum(0.64f));
        m_sigma_ms = props.getSpectrum("sigma_ms", Spectrum(1.69f));
        m_sigma_ma = props.getSpectrum("sigma_ma", Spectrum(0.17f));

        _sigma_ca = m_sigma_ca.getLuminance();
        _sigma_ms = m_sigma_ms.getLuminance();
        _sigma_ma = m_sigma_ma.getLuminance();
		
        m_g  = props.getFloat("g", 0.44f);
        m_l  = props.getFloat("l", 0.47f);

        _alpha[0] = m_alpha;
        _alpha[1] = -m_alpha/2;
        _alpha[2] = -3*m_alpha/2;
        _beta_m[0] = m_beta_m;
        _beta_m[1] = m_beta_m/2;
        _beta_m[2] = 3*m_beta_m/2;
        _beta_n[0] = m_beta_n;
        _beta_n[1] = m_beta_n * std::sqrt(2.0);
        _beta_n[2] = m_beta_n * std::sqrt(3.0);

        precomputeAzimuthalDistributions();

        _ms = new PrecomputedScatterLobe("data/fur/medulla_longitudinal.bin");

        sample_profile();

        /*for(float i = 0.1f; i<3.0f; i = i + 0.2f){
            std::cout<<"pre-compute: " << i  << ":" << _nR->eval(i, 0.1f).toString() << " " << 
            _nTT->eval(i, 0.1f).toString() << " " <<
            _nTRT->eval(i, 0.1f).toString() << std::endl;
        }*/
    }

    YanFurRT(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {

        m_kappa = stream->readFloat();
		m_eta = stream->readFloat();
        m_alpha = stream->readFloat();
        m_beta_m = stream->readFloat();
        m_beta_n = stream->readFloat();

		m_sigma_ca = Spectrum(stream);
        m_sigma_ms = Spectrum(stream);
        m_sigma_ma = Spectrum(stream);

        m_g = stream->readFloat();
        m_l = stream->readFloat();

        configure();
    }

    void configure() {

        m_components.clear();

        m_components.push_back(EDiffuseReflection | EFrontSide);
        m_usesRayDifferentials = false;

        precomputeAzimuthalDistributions();

        BSDF::configure();
    }

    Spectrum brdf(Float thetaI, Float thetaR, Float Phi) const {

        if (Phi < 0.0f)
            Phi += M_PI*2.0f;
        if (Phi > M_PI*2.0f)
            Phi -= M_PI*2.0f;

        Float thetaD = std::abs(thetaR - thetaI)*0.5f;
        Float cosThetaD = std::cos(thetaD);

        auto Mp = M(thetaR, thetaI);
        auto Np = N(Phi, cosThetaD);

        auto unscatter =  Mp[0] *  Np[0]
                        + Mp[1] *   Np[1]
                        + Mp[2] *   Np[2];

        auto scatter = Mp[3] * (Np[3] + Np[4]);

        Float cosThetaI = std::cos(thetaI);
        Float cosSqrtThetaI = std::max(cosThetaI*cosThetaI, 0.1f);
        return  (unscatter+scatter)/(cosSqrtThetaI);
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure=ESolidAngle) const {
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
        
        return brdf(thetaI, thetaR, phi);
    }

    unsigned char* remap_pseudo_color(float val){
        static unsigned char color[3] = {0,0,0};
        if (val <= 0.5f)
        {
            val *= 2.0f;
            color[0] = (unsigned char)(unsigned int)(0);
            color[1] = (unsigned char)(unsigned int)(255 * (val) + 0.5f);
            color[2] = (unsigned char)(unsigned int)(255 * (1.0f - val) + 0.5f);
        }
        else
        {
            val = val * 2.0f - 1.0f;
            color[0] = (unsigned char)(unsigned int)(255 * (val) + 0.5f);
            color[1] = (unsigned char)(unsigned int)(255 * (1.0f - val) + 0.5f);
            color[2] = (unsigned char)(unsigned int)(0);
        }
        return color;
    }

    void sample_profile(){

        int stepTheta = 256;
        int stepPhi = 512;

        Float thetaI = degToRad(-40);
        Float thetaInterval = 40.0f/stepTheta;
        Float phiInterval = 180.0f/stepPhi;

        FILE* fp = fopen("profile_fur.ppm", "wb");
        (void)fprintf(fp, "P6\n%d %d\n255\n", stepPhi,stepTheta);
        for(int i = stepTheta; i >0; i--){
            for(int j = stepPhi; j >0; j--){
                Float thetaO = degToRad(10.0f+thetaInterval*i);
                Float phi = degToRad(phiInterval*j);
                Spectrum observe = brdf(thetaI, thetaO, phi);
                Float lum = observe.getLuminance();
                lum = math::clamp(0.0f, 1.0f, lum*10.0f);
                auto color = remap_pseudo_color(lum);
                fwrite(color, 1, 3, fp);
            }
        }
        fclose(fp);   

        PrecomputedScatterLobe* _ns = new PrecomputedScatterLobe("data/fur/medulla_azimuthal.bin");
        int step = 256;
        Float stepInterval = 1.0f/step;
        FILE* fp2 = fopen("profile_CN.ppm", "wb");
        (void)fprintf(fp2, "P6\n%d %d\n255\n", step,step);
        for(int i = step; i >0; i--){
            for(int j = step; j >0; j--){
                Float val = _ns->eval(_sigma_ms, stepInterval*i, m_g, stepInterval*j);
                auto color = remap_pseudo_color(val);
                fwrite(color, 1, 3, fp2);
            }
        }
        fclose(fp2);   

        Float plotCNSigmaInterval = 1.0f/8;
        Float phiCNPlotInterval = 1.0f/100;
        FILE* fp_cn_txt = fopen("profile_cn.txt", "w");
        for(int j = 8; j >=0; j--){
            for(int i = 100; i >=0; i--){
                Float val = _ns->eval(plotCNSigmaInterval*j, 0.8f, 0.5f, phiCNPlotInterval*i);
                fprintf(fp_cn_txt, "sigma\t%f\tphi\t%f\tval\t%f\n", 20*plotCNSigmaInterval*j,M_PI*2*phiCNPlotInterval*i, val);
            }
            fprintf(fp_cn_txt, "\n");
        }
        fclose(fp_cn_txt);   

        FILE* fp3 = fopen("profile_CM.ppm", "wb");
        (void)fprintf(fp3, "P6\n%d %d\n255\n", step,step);
        for(int i = step; i >0; i--){
            for(int j = step; j >0; j--){
                Float val = _ms->eval(_sigma_ms, stepInterval*i, m_g, stepInterval*j);
                auto color = remap_pseudo_color(val);
                fwrite(color, 1, 3, fp3);
            }
        }
        fclose(fp3);   

        Float plotCMSigmaInterval = 1.0f/8;
        Float thetaOCMPlotInterval = 1.0f/100;
        FILE* fp_cm_txt = fopen("profile_cm.txt", "w");
        for(int j = 8; j >=0; j--){
            for(int i = 100; i >=0; i--){
                Float val = _ms->eval(plotCMSigmaInterval*j, 0.2f, 0.5f, thetaOCMPlotInterval*i);
                fprintf(fp_cm_txt, "sigma\t%f\tthetaR\t%f\tval\t%f\n", 20*plotCMSigmaInterval*j,M_PI*2*thetaOCMPlotInterval*i, val);
            }
            fprintf(fp_cm_txt, "\n");
        }
        fclose(fp_cm_txt);   

        Float pi2step = 2*M_PI/step;
        FILE* fp4 = fopen("profile_ntrt.ppm", "wb");
        (void)fprintf(fp4, "P6\n%d %d\n255\n", step,step);
        for(int i = step; i >0; i--){
            for(int j = step; j >0; j--){
                Float val = _nTRT->eval(pi2step*i, stepInterval*j).getLuminance();
                auto color = remap_pseudo_color(val);
                fwrite(color, 1, 3, fp4);
            }
        }
        fclose(fp4);   

        Float plotInterval = 1.0f/4;
        Float phiPlotInterval = 2*M_PI/100;
        FILE* fp_ntt_txt = fopen("profile_ntrt.txt", "w");
        for(int j = 4; j >=0; j--){
            for(int i = 100; i >=0; i--){
                Float val = _nTRT->eval(phiPlotInterval*i, plotInterval*j).getLuminance();
                fprintf(fp_ntt_txt, "Phi\t%f\tcosThetaD\t%f\tval\t%f\n", phiPlotInterval*i,plotInterval*j, val);
            }
            fprintf(fp_ntt_txt, "\n");
        }
        fclose(fp_ntt_txt);   

        FILE* fp5 = fopen("profile_ntts.ppm", "wb");
        (void)fprintf(fp5, "P6\n%d %d\n255\n", step,step);
        for(int i = step; i >0; i--){
            for(int j = step; j >0; j--){
                Float val = _nTTs->eval(phiPlotInterval*i, stepInterval*j).getLuminance();
                auto color = remap_pseudo_color(val);
                fwrite(color, 1, 3, fp5);
            }
        }
        fclose(fp5);    


        FILE* fp_ntts_txt = fopen("profile_ntts.txt", "w");
        for(int j = 4; j >=0; j--){
            for(int i = 100; i >=0; i--){
                Float val = _nTTs->eval(phiPlotInterval*i, plotInterval*j).getLuminance();
                fprintf(fp_ntts_txt, "Phi\t%f\tcosThetaD:%f\tval\t%f\n", phiPlotInterval*i,plotInterval*j, val);
            }
            fprintf(fp_ntts_txt, "\n");
        }
        fclose(fp_ntts_txt);   

        FILE* fp6 = fopen("profile_ntrts.ppm", "wb");
        (void)fprintf(fp6, "P6\n%d %d\n255\n", step,step);
        for(int i = step; i >0; i--){
            for(int j = step; j >0; j--){
                Float val = _nTRTs->eval(phiPlotInterval*i, stepInterval*j).getLuminance();
                auto color = remap_pseudo_color(val);
                fwrite(color, 1, 3, fp6);
            }
        }
        fclose(fp6);   

        FILE* fp_ntrts_txt = fopen("profile_ntrts.txt", "w");
        for(int j = 4; j >=0; j--){
            for(int i = 100; i >=0; i--){
                Float val = _nTRTs->eval(phiPlotInterval*i, plotInterval*j).getLuminance();
                fprintf(fp_ntrts_txt, "Phi\t%f\tval\t%f\n", phiPlotInterval*i,val);
            }
            fprintf(fp_ntrts_txt, "\n");
        }
        fclose(fp_ntrts_txt);   


    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle)
            return 0.0f;

        return warp::squareToUniformSpherePdf();
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        if (!(bRec.typeMask & EDiffuseReflection))
            return Spectrum(0.0f);

        Float _pdf;
        return this->sample(bRec,_pdf,sample);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        if (!(bRec.typeMask & EDiffuseReflection))
            return Spectrum(0.0f);

        bRec.wo = warp::squareToUniformSphere(sample);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;
        pdf = warp::squareToUniformSpherePdf();
        return eval(bRec);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        BSDF::addChild(name, child);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        stream->writeFloat(m_kappa);
		stream->writeFloat(m_eta);
        stream->writeFloat(m_alpha);
        stream->writeFloat(m_beta_m);
        stream->writeFloat(m_beta_n);

		m_sigma_ca.serialize(stream);
        m_sigma_ms.serialize(stream);
        m_sigma_ma.serialize(stream);

        stream->writeFloat(m_g);
        stream->writeFloat(m_l);

        //manager->serialize(stream, m_diffuseReflectance.get());
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "Marschner Hair[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  sigma = \"" << m_sigma_ca.toString() << "\"," << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()

private:
    // Considering l layers together, the reflectance is then given by [Stokes 1860]
    // Eq.7 in Physically-Accurate Fur Reflectance: Modeling, Measurement and Rendering 
    inline Float dielectricReflectanceLayers(Float eta, Float cosThetaI, Float layer) const
    {
        Float fresnel = dielectricReflectance(eta, cosThetaI);
        return (fresnel * layer) / (1 + (layer-1) * fresnel);
    }

    inline Float dielectricReflectance(Float eta, Float cosThetaI) const
    {
        Float cosThetaT;
        if (cosThetaI < 0.0f) {
            eta = 1.0f/eta;
            cosThetaI = -cosThetaI;
        }
        Float sinThetaTSq = eta*eta*(1.0f - cosThetaI*cosThetaI);
        if (sinThetaTSq > 1.0f) {
            cosThetaT = 0.0f;
            return 1.0f;
        }
        cosThetaT = std::sqrt(std::max(1.0f - sinThetaTSq, 0.0f));

        Float Rs = (eta*cosThetaI - cosThetaT)/(eta*cosThetaI + cosThetaT);
        Float Rp = (eta*cosThetaT - cosThetaI)/(eta*cosThetaT + cosThetaI);

        return (Rs*Rs + Rp*Rp)*0.5f;
    }

    inline Float SchlickFresnel(Float eta, Float cosThetaI) const {
        Float F0 = (1-eta)/(1+eta);
        F0 = F0*F0;
        return F0 + (1-F0)*std::pow((1-cosThetaI), 5);
    }

    inline Float SchlickFresnelLayer(Float eta, Float cosThetaI, Float layer) const {
        Float fresnel = SchlickFresnel(eta, cosThetaI);
        return (fresnel * layer) / (1 + (layer-1) * fresnel);
    }

    std::vector<Spectrum> N(Float Phi, Float cosThetaD) const {
        std::vector<Spectrum> result(5);

        //dEon 14, epic 16
        auto cosPhi = std::cos(Phi);
        auto cosHalfPhi = std::sqrt(0.5f+0.5f*cosPhi);

        //assume eta = 1.55, epic 16
        auto n_prime = 1.19 / cosThetaD + 0.36 * cosThetaD;
        float a = 1 / n_prime;
        float h = cosHalfPhi * (1 + a * (0.6 - 0.8 * cosPhi));



        auto Nr = 0.25f * math::clamp(std::sqrt(0.5f+0.5f*cosPhi), 0.0f, 1.0f);
        auto Ar = SchlickFresnelLayer(m_eta, std::sqrt(0.5f+0.5f*cosThetaD), m_l);
        result[0] = Spectrum(Nr*Ar);


        auto sinGammaTSqr = h*h*a*a;
        auto sm = math::safe_sqrt(m_kappa*m_kappa - sinGammaTSqr);
        auto sc = math::safe_sqrt(1-sinGammaTSqr)-sm;
        auto Cc = (-m_sigma_ca*4).exp();
        auto Cma = (-m_sigma_ma*4).exp();
        auto Cms = (-m_sigma_ms*4).exp();
        auto Cm = Cma * Cms;
        Spectrum Tc = Cc.pow(sc/(2*cosThetaD));
        Spectrum Tm = Cm.pow(sm/(2*cosThetaD));
        
        //approximate of Dtt, epic 16
        //float Ntt = Logistic(Phi, 0.4, 3.14);
        float Ntt = exp(-2.6 * cosPhi - 3.1);
        float Ftt = SchlickFresnelLayer(m_eta, cosThetaD * std::sqrt(math::clamp(1.0f - h * h,0.0f,1.0f)), m_l);
        Ftt = (1-Ftt)*(1-Ftt);
        result[1] = Tc * Ntt * Ftt * Tm;

        //we assume h = sqrt(3)/2
        sm = math::safe_sqrt(m_kappa*m_kappa-3.0f/4.0f);
        sc = 1/2.0f - sm;
        //brute force fit when beta_n = 12
        float Ntrt = std::exp((6.3f*cosThetaD+0.7f)*cosPhi-(5*cosThetaD+2));
        float Fhackh = SchlickFresnelLayer(m_eta, cosThetaD * 0.5f, m_l);
        float Ftrt = Fhackh * (1-Fhackh)*(1-Fhackh);
        Tc = Cc.pow(sc/(cosThetaD));
        Tm = Cm.pow(sm/(cosThetaD));

        result[2] = Tc * Tm * Ntrt * Ftrt;

        Spectrum Atts = Cc.pow((sc+1-m_kappa)/(4*cosThetaD))*Cma.pow(m_kappa/(4*cosThetaD));
        float Ntts = 0.05*std::cos(2*Phi) + 0.16f;

        result[3] = Ntts*Atts*Fhackh;

        float Ntrts = 0.08*std::cos(1.5*Phi + 1.7) + 0.18f;
        Spectrum Atrts = Cc.pow((3*sc+1-m_kappa)/(4*cosThetaD))*Cma.pow((2*sm+m_kappa)/(4*cosThetaD))*Cms.pow(sm/(8*cosThetaD));

        result[4] = Atrts * Fhackh * (1-Fhackh) * Ntrts;
        
        return result;
    }

    // Standard Guassian
    Float gaussian(Float x, Float mu, Float sigma) const
    {	
        return (Float)(std::exp(-(x-mu)*(x-mu)/(2*sigma*sigma))
            / (std::sqrt(2*M_PI) * std::abs(sigma)));
    }

    Float Phi(float gammaI, float gammaT, int p)
    {
        return 2.0f*p*gammaT - 2.0f*gammaI + p*M_PI;
    }

    Float Phis(float gammaI, float gammaT, int p)
    {
        return gammaT - gammaI + (p-1)*(M_PI+2*gammaT);
    }

    std::vector<Float> M(Float thetaO, Float thetaI) const {
        std::vector<Float> result(4);
        for (int i = 0; i < 3; ++i) {
            result[i] = gaussian(thetaO, -thetaI + _alpha[i], _beta_m[i]);
        }

        Float thetaO0 = M_PI_2 + thetaO;
        result[3] = std::abs(std::sin(thetaO0))*0.5;
        return result;
    }

    // Standard normalized Gaussian
    float g(float beta, float theta)
    {
        return std::exp(-theta*theta/(2.0f*beta*beta))/(std::sqrt(2.0f*M_PI)*std::abs(beta));
    }

    float D(float beta, float phi)
    {
        float result = 0.0f;
        float delta;
        float shift = 0.0f;
        do {
            delta = g(beta, phi + shift) + g(beta, phi - shift - 2.0f*M_PI);
            result += delta;
            shift += 2.0f*M_PI;
        } while (delta > 1e-4f);
        return result;
    }


    void precomputeAzimuthalDistributions()
    {
        const int Resolution = PrecomputedAzimuthalLobe::AzimuthalResolution;
        std::unique_ptr<Spectrum[]> valuesR  (new Spectrum[Resolution*Resolution]);
        std::unique_ptr<Spectrum[]> valuesTT (new Spectrum[Resolution*Resolution]);
        std::unique_ptr<Spectrum[]> valuesTRT(new Spectrum[Resolution*Resolution]);

        std::unique_ptr<Spectrum[]> valuesTTs (new Spectrum[Resolution*Resolution]);
        std::unique_ptr<Spectrum[]> valuesTRTs(new Spectrum[Resolution*Resolution]);

        PrecomputedScatterLobe* _ns = new PrecomputedScatterLobe("data/fur/medulla_azimuthal.bin");

        // Ideally we could simply make this a constexpr, but MSVC does not support that yet (boo!)
        #define NumPoints 140

        // Integrate over h
        GaussLegendre<NumPoints> integrator;
        const auto points = integrator.points();
        const auto weights = integrator.weights();

        // Cache the gammaI across all integration points
        std::array<float, NumPoints> gammaIs;
        for (int i = 0; i < NumPoints; ++i)
            gammaIs[i] = std::asin(points[i]);

        const int NumGaussianSamples = 2048;
        std::unique_ptr<float[]> Ds[3];
        for (int p = 0; p < 3; ++p) {
            Ds[p].reset(new float[NumGaussianSamples]);
            for (int i = 0; i < NumGaussianSamples; ++i)
                Ds[p][i] = D(_beta_n[p], i/(NumGaussianSamples - 1.0f)*M_PI*2.0f);
        }

        // Simple wrapped linear interpolation of the precomputed table
        auto approxD = [&](int p, float phi) {
            float u = std::abs(phi*(0.5f*M_1_PI*(NumGaussianSamples - 1)));
            int x0 = int(u);
            int x1 = x0 + 1;
            u -= x0;
            return Ds[p][x0 % NumGaussianSamples]*(1.0f - u) + Ds[p][x1 % NumGaussianSamples]*u;
        };

        // Here follows the actual precomputation of the azimuthal scattering functions
        // The scattering functions are parametrized with the azimuthal angle phi,
        // and the cosine of the half angle, cos(thetaD).
        // This parametrization makes the azimuthal function relatively smooth and allows using
        // really low resolutions for the table (64x64 in this case) without any visual
        // deviation from ground truth, even at the lowest supported roughness setting
        for (int y = 0; y < Resolution; ++y) {
            float cosHalfAngle = y/(Resolution - 1.0f);

            // Precompute reflection Fresnel factor and reduced absorption coefficient
            float iorPrime = std::sqrt(m_eta*m_eta - (1.0f - cosHalfAngle*cosHalfAngle))/cosHalfAngle;
            float cosThetaT = std::sqrt(1.0f - (1.0f - cosHalfAngle*cosHalfAngle)*math::sqr(1.0f/m_eta));
            //Spectrum sigma_cPrime = m_sigma_ca/cosThetaT;
            //Spectrum sigma_mPrime = (m_sigma_ms + m_sigma_ma)/cosThetaT;

            // Precompute gammaT, f_t and internal absorption across all integration points
            std::array<float, NumPoints> fresnelTerms, gammaTs;
            std::array<Spectrum, NumPoints> absorptions_c; // Tc integrator
            std::array<Spectrum, NumPoints> absorptions_m; // Tm integrator
            std::array<Spectrum, NumPoints> absorptions_tts; // Tc integrator
            std::array<Spectrum, NumPoints> absorptions_trts; // Tm integrator
            for (int i = 0; i < NumPoints; ++i) {
                auto sinGammaT = math::clamp(points[i]/iorPrime, -1.0f, 1.0f);
                gammaTs[i] = std::asin(sinGammaT);
                //fresnelTerms[i] = dielectricReflectanceLayers(1.0f/iorPrime, std::cos(gammaIs[i]), m_l);
                fresnelTerms[i] = dielectricReflectanceLayers(1.0f/m_eta, cosHalfAngle*std::cos(gammaIs[i]), m_l);
                auto sm = math::safe_sqrt(m_kappa*m_kappa - sinGammaT*sinGammaT);
                auto sc = math::safe_sqrt(1 - sinGammaT*sinGammaT)-sm;
                absorptions_c[i] = (-m_sigma_ca*2.0f*sc/cosThetaT).exp();
                absorptions_m[i] = (-(m_sigma_ms + m_sigma_ma)*2.0f*sm/cosThetaT).exp();
                absorptions_tts[i] = (-((sc+1-m_kappa)*m_sigma_ca+m_kappa*m_sigma_ma)/cosThetaT).exp();
                absorptions_trts[i] = (-((3*sc+1-m_kappa)*m_sigma_ca+(2*sm+m_kappa)*m_sigma_ma+2*sm*m_sigma_ms)/cosThetaT).exp();
            }

            for (int phiI = 0; phiI < Resolution; ++phiI) {
                float phi = M_PI*2.0f*phiI/(Resolution - 1.0f);

                float integralR = 0.0f;
                Spectrum integralTT(0.0f);
                float integralTRT(0.0f);
                float integralTTs(0.0f);
                float integralTRTs(0.0f);

                // Here follows the integration across the fiber width, h.
                // Since we were able to precompute most of the factors that
                // are constant w.r.t phi for a given h,
                // we don't have to do much work here.
                for (int i = 0; i < integrator.numSamples(); ++i) {
                    float fR = fresnelTerms[i];
                    Spectrum Tc = absorptions_c[i];
                    Spectrum Tm = absorptions_m[i];
                    Spectrum Ttts = absorptions_tts[i];
                    Spectrum Ttrts = absorptions_trts[i];

                    float AR = fR;
                    Spectrum ATT = (1.0f - fR)*(1.0f - fR)*Tc*Tm;
                    Spectrum ATRT = ATT*fR*Tc*Tm;
                    Spectrum ATTs = fR*Ttts;
                    Spectrum ATRTs = (1.0f - fR)*fR*Ttrts;

                    integralR    += weights[i]*approxD(0, phi - Phi(gammaIs[i], gammaTs[i], 0))*AR;
                    integralTT   += weights[i]*approxD(1, phi - Phi(gammaIs[i], gammaTs[i], 1))*ATT;
                    integralTRT  += weights[i]*approxD(2, phi - Phi(gammaIs[i], gammaTs[i], 2));//*ATRT;
                    auto hPrime = math::clamp(points[i]/iorPrime, -1.0f, 1.0f); //sin gammaT
                    integralTTs  += weights[i]*_ns->evalN(_sigma_ms/m_kappa, hPrime/m_kappa, m_g, phi - Phis(gammaIs[i], gammaTs[i], 1));//*ATTs;
                    integralTRTs += weights[i]*_ns->evalN(_sigma_ms/m_kappa, hPrime/m_kappa, m_g, phi - Phis(gammaIs[i], gammaTs[i], 2));//*ATRTs;
                }
                //std::cout<<"IntegralTT: "<<cosHalfAngle<<","<<phi<<":\t"<<integralTT.toString()<<std::endl;
                valuesR   [phiI + y*Resolution] = Spectrum(0.5f*integralR);
                valuesTT  [phiI + y*Resolution] = 0.5f*integralTT;
                valuesTRT [phiI + y*Resolution] = Spectrum(0.5f*integralTRT);
                valuesTTs [phiI + y*Resolution] = Spectrum(0.5f*integralTTs);
                valuesTRTs[phiI + y*Resolution] = Spectrum(0.5f*integralTRTs);
            }
        }

        // Hand the values off to the helper class to construct sampling CDFs and so forth
        _nR    = new PrecomputedAzimuthalLobe(std::move(valuesR));
        _nTT   = new PrecomputedAzimuthalLobe(std::move(valuesTT));
        _nTRT  = new PrecomputedAzimuthalLobe(std::move(valuesTRT));
        _nTTs  = new PrecomputedAzimuthalLobe(std::move(valuesTTs));
        _nTRTs = new PrecomputedAzimuthalLobe(std::move(valuesTRTs));
    }

    Float m_eta;			// refractive index of cortex and medulla
    Float m_kappa;          // medullary index (rel. radius length)
	Spectrum m_sigma_ca;	// absorption coefficient in cortex
    Spectrum m_sigma_ms;	// scattering coefficient in medulla
    Spectrum m_sigma_ma;	// absorption coefficient in medulla
    Float m_alpha;		    // scale tilt for cuticle
    Float m_g;              //anisotropy factor of scattering in medulla
    Float m_l;              // layers of cuticle

    Float m_beta_m;         // longitudinal roughness of cuticle
    Float m_beta_n;         // azimuthal roughness of cuticle 

    Float _alpha[3], _beta_m[3], _beta_n[3];
    Float _sigma_ca, _sigma_ms, _sigma_ma;

    PrecomputedAzimuthalLobe* _nR;
    PrecomputedAzimuthalLobe* _nTT;
    PrecomputedAzimuthalLobe* _nTRT;
    PrecomputedAzimuthalLobe* _nTTs;
    PrecomputedAzimuthalLobe* _nTRTs;
    PrecomputedScatterLobe* _ms;

    //ref<Texture> m_diffuseReflectance;
};

MTS_IMPLEMENT_CLASS_S(YanFurRT, false, BSDF)
MTS_EXPORT_PLUGIN(YanFurRT, "Yan Fur RT")
MTS_NAMESPACE_END
