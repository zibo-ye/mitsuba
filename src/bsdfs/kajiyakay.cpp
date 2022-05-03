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
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

class KajiyaKay : public BSDF {
public:
    KajiyaKay(const Properties &props)
        : BSDF(props) {
        m_diffuseReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("diffuseReflectance", Spectrum(0.5f)));
        m_specularReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("specularReflectance", Spectrum(0.2f)));
        m_exponent = new ConstantFloatTexture(
            props.getFloat("exponent", 30.0f));
        m_specularSamplingWeight = 0.0f;
    }

    KajiyaKay(Stream *stream, InstanceManager *manager)
     : BSDF(stream, manager) {
        m_diffuseReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_exponent = static_cast<Texture *>(manager->getInstance(stream));

        configure();
    }

    void configure() {
        m_components.clear();
        m_components.push_back(EGlossyReflection | EFrontSide |
            ((!m_specularReflectance->isConstant()
              || !m_exponent->isConstant()) ? ESpatiallyVarying : 0));
        m_components.push_back(EDiffuseReflection | EFrontSide
            | (m_diffuseReflectance->isConstant() ? 0 : ESpatiallyVarying));

        /* Verify the input parameters and fix them if necessary */
        std::pair<Texture *, Texture *> result = ensureEnergyConservation(
            m_specularReflectance, m_diffuseReflectance,
            "specularReflectance", "diffuseReflectance", 1.0f);
        m_specularReflectance = result.first;
        m_diffuseReflectance = result.second;

        /* Compute weights that steer samples towards
           the specular or diffuse components */
        Float dAvg = m_diffuseReflectance->getAverage().getLuminance(),
              sAvg = m_specularReflectance->getAverage().getLuminance();
        m_specularSamplingWeight = sAvg / (dAvg + sAvg);

        m_usesRayDifferentials =
            m_diffuseReflectance->usesRayDifferentials() ||
            m_specularReflectance->usesRayDifferentials() ||
            m_exponent->usesRayDifferentials();

        BSDF::configure();
    }

    Spectrum getDiffuseReflectance(const Intersection &its) const {
        return m_diffuseReflectance->eval(its);
    }

    Spectrum getSpecularReflectance(const Intersection &its) const {
        return m_specularReflectance->eval(its);
    }

    /// Reflection in local coordinates
    inline Vector reflect(const Vector &wi) const {
        return Vector(-wi.x, -wi.y, wi.z);
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 || measure != ESolidAngle)
            return Spectrum(0.0f);

        bool hasSpecular = (bRec.typeMask & EGlossyReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
                && (bRec.component == -1 || bRec.component == 1);

        Spectrum result(0.0f);
        if (hasSpecular) {
            Float alpha    = dot(bRec.wo, reflect(bRec.wi)),
                  exponent = m_exponent->eval(bRec.its).average();

            if (alpha > 0.0f) {
                result += m_specularReflectance->eval(bRec.its) *
                    ((exponent + 2) * INV_TWOPI * std::pow(alpha, exponent));
            }
        }

        if (hasDiffuse)
            result += m_diffuseReflectance->eval(bRec.its) * INV_PI;

        return result * Frame::cosTheta(bRec.wo);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 || measure != ESolidAngle)
            return 0.0f;

        bool hasSpecular = (bRec.typeMask & EGlossyReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
                && (bRec.component == -1 || bRec.component == 1);

        Float diffuseProb = 0.0f, specProb = 0.0f;

        if (hasDiffuse)
            diffuseProb = warp::squareToCosineHemispherePdf(bRec.wo);

        if (hasSpecular) {
            Float alpha    = dot(bRec.wo, reflect(bRec.wi)),
                  exponent = m_exponent->eval(bRec.its).average();
            if (alpha > 0)
                specProb = std::pow(alpha, exponent) *
                    (exponent + 1.0f) / (2.0f * M_PI);
        }

        if (hasDiffuse && hasSpecular)
            return m_specularSamplingWeight * specProb +
                   (1-m_specularSamplingWeight) * diffuseProb;
        else if (hasDiffuse)
            return diffuseProb;
        else if (hasSpecular)
            return specProb;
        else
            return 0.0f;
    }

    inline Spectrum sample(BSDFSamplingRecord &bRec, Float &_pdf, const Point2 &_sample) const {
        Point2 sample(_sample);

        bool hasSpecular = (bRec.typeMask & EGlossyReflection)
                && (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse  = (bRec.typeMask & EDiffuseReflection)
                && (bRec.component == -1 || bRec.component == 1);

        if (!hasSpecular && !hasDiffuse)
            return Spectrum(0.0f);

        bool choseSpecular = hasSpecular;

        if (hasDiffuse && hasSpecular) {
            if (sample.x <= m_specularSamplingWeight) {
                sample.x /= m_specularSamplingWeight;
            } else {
                sample.x = (sample.x - m_specularSamplingWeight)
                    / (1-m_specularSamplingWeight);
                choseSpecular = false;
            }
        }

        if (choseSpecular) {
            Vector R = reflect(bRec.wi);
            Float exponent = m_exponent->eval(bRec.its).average();

            /* Sample from a KajiyaKay lobe centered around (0, 0, 1) */
            Float sinAlpha = std::sqrt(1-std::pow(sample.y, 2/(exponent + 1)));
            Float cosAlpha = std::pow(sample.y, 1/(exponent + 1));
            Float phi = (2.0f * M_PI) * sample.x;
            Vector localDir = Vector(
                sinAlpha * std::cos(phi),
                sinAlpha * std::sin(phi),
                cosAlpha
            );

            /* Rotate into the correct coordinate system */
            bRec.wo = Frame(R).toWorld(localDir);
            bRec.sampledComponent = 1;
            bRec.sampledType = EGlossyReflection;

            if (Frame::cosTheta(bRec.wo) <= 0)
                return Spectrum(0.0f);
        } else {
            bRec.wo = warp::squareToCosineHemisphere(sample);
            bRec.sampledComponent = 0;
            bRec.sampledType = EDiffuseReflection;
        }
        bRec.eta = 1.0f;

        _pdf = pdf(bRec, ESolidAngle);

        if (_pdf == 0)
            return Spectrum(0.0f);
        else
            return eval(bRec, ESolidAngle) / _pdf;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        Float pdf;
        return KajiyaKay::sample(bRec, pdf, sample);
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
            if (name == "exponent")
                m_exponent = static_cast<Texture *>(child);
            else if (name == "specularReflectance")
                m_specularReflectance = static_cast<Texture *>(child);
            else if (name == "diffuseReflectance")
                m_diffuseReflectance = static_cast<Texture *>(child);
            else
                BSDF::addChild(name, child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        manager->serialize(stream, m_diffuseReflectance.get());
        manager->serialize(stream, m_specularReflectance.get());
        manager->serialize(stream, m_exponent.get());
    }

    Float getRoughness(const Intersection &its, int component) const {
        Assert(component == 0 || component == 1);
        /* Find the Beckmann-equivalent roughness */
        if (component == 0)
            return std::sqrt(2 / (2+m_exponent->eval(its).average()));
        else
            return std::numeric_limits<Float>::infinity();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "KajiyaKay[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  diffuseReflectance = " << indent(m_diffuseReflectance->toString()) << "," << endl
            << "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
            << "  specularSamplingWeight = " << m_specularSamplingWeight << "," << endl
            << "  diffuseSamplingWeight = " << (1-m_specularSamplingWeight) << "," << endl
            << "  exponent = " << indent(m_exponent->toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    ref<Texture> m_diffuseReflectance;
    ref<Texture> m_specularReflectance;
    ref<Texture> m_exponent;
    Float m_specularSamplingWeight;
};

MTS_IMPLEMENT_CLASS_S(KajiyaKay, false, BSDF)
MTS_EXPORT_PLUGIN(KajiyaKay, "KajiyaKay Hair BRDF");
MTS_NAMESPACE_END
