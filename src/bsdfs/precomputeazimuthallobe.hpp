#ifndef PRECOMPUTEAZIMU_HPP_
#define PRECOMPUTEAZIMU_HPP_

#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

class InterpolatedDistribution1D
{
    int _size;
    int _numDistributions;

    std::vector<float> _pdfs;
    std::vector<float> _cdfs;
    std::vector<float> _sums;

    float cdfs(int x, int distribution) const
    {
        return _cdfs[x + distribution * (_size + 1)];
    }

    float pdfs(int x, int distribution) const
    {
        return _pdfs[x + distribution * _size];
    }

    float &cdfs(int x, int distribution)
    {
        return _cdfs[x + distribution * (_size + 1)];
    }

    float &pdfs(int x, int distribution)
    {
        return _pdfs[x + distribution * _size];
    }

public:
    InterpolatedDistribution1D(std::vector<float> weights, int size, int numDistributions)
        : _size(size),
          _numDistributions(numDistributions),
          _pdfs(std::move(weights)),
          _cdfs((size + 1) * numDistributions),
          _sums(numDistributions)
    {
        for (int dist = 0; dist < _numDistributions; ++dist)
        {
            cdfs(0, dist) = 0.0f;
            for (int x = 0; x < _size; ++x)
                cdfs(x + 1, dist) = pdfs(x, dist) + cdfs(x, dist);

            _sums[dist] = cdfs(_size, dist);

            if (_sums[dist] < 1e-4f)
            {
                // Revert to uniform sampling for near-degenerate distributions
                float ratio = 1.0f / _size;
                for (int x = 0; x < _size; ++x)
                {
                    pdfs(x, dist) = ratio;
                    cdfs(x, dist) = x * ratio;
                }
            }
            else
            {
                float scale = 1.0f / _sums[dist];
                for (int x = 0; x < _size; ++x)
                {
                    pdfs(x, dist) *= scale;
                    cdfs(x, dist) *= scale;
                }
            }
            cdfs(_size, dist) = 1.0f;
        }
    }

    void warp(float distribution, float &u, int &x) const
    {
        int d0 = math::clamp(int(distribution), 0, _numDistributions - 1);
        int d1 = std::min(d0 + 1, _numDistributions - 1);
        float v = math::clamp(distribution - d0, 0.0f, 1.0f);

        int lower = 0, upper = _size;
        float lowerU = 0.0f, upperU = 1.0f;

        while (upper - lower != 1)
        {
            int midpoint = (upper + lower) / 2;
            float midpointU = cdfs(midpoint, d0) * (1.0f - v) + cdfs(midpoint, d1) * v;
            if (midpointU < u)
            {
                lower = midpoint;
                lowerU = midpointU;
            }
            else
            {
                upper = midpoint;
                upperU = midpointU;
            }
        }

        x = lower;
        u = math::clamp((u - lowerU) / (upperU - lowerU), 0.0f, 1.0f);
    }

    float pdf(float distribution, int x) const
    {
        int d0 = math::clamp(int(distribution), 0, _numDistributions - 1);
        int d1 = std::min(d0 + 1, _numDistributions - 1);
        float v = math::clamp(distribution - d0, 0.0f, 1.0f);

        return pdfs(x, d0) * (1.0f - v) + pdfs(x, d1) * v;
    }

    float sum(float distribution) const
    {
        int d0 = math::clamp(int(distribution), 0, _numDistributions - 1);
        int d1 = std::min(d0 + 1, _numDistributions - 1);
        float v = math::clamp(distribution - d0, 0.0f, 1.0f);

        return _sums[d0] * (1.0f - v) + _sums[d1] * v;
    }
};

class PrecomputedAzimuthalLobe
{
public:
    static const int AzimuthalResolution = 64;

private:
    std::unique_ptr<Spectrum[]> _table;
    std::unique_ptr<InterpolatedDistribution1D> _sampler;

public:
    PrecomputedAzimuthalLobe(std::unique_ptr<Spectrum[]> table)
        : _table(std::move(table))
    {
        const int Size = AzimuthalResolution;

        std::vector<float> weights(Size * Size);
        for (int i = 0; i < Size * Size; ++i)
            weights[i] = _table[i].max();

        // Dilate weights slightly to stay conservative
        for (int y = 0; y < Size; ++y)
        {
            for (int x = 0; x < Size - 1; ++x)
                weights[x + y * Size] = std::max(weights[x + y * Size], weights[x + 1 + y * Size]);
            for (int x = Size - 1; x > 0; --x)
                weights[x + y * Size] = std::max(weights[x + y * Size], weights[x - 1 + y * Size]);
        }
        for (int x = 0; x < Size; ++x)
        {
            for (int y = 0; y < Size - 1; ++y)
                weights[x + y * Size] = std::max(weights[x + y * Size], weights[x + (y + 1) * Size]);
            for (int y = Size - 1; y > 0; --y)
                weights[x + y * Size] = std::max(weights[x + y * Size], weights[x + (y - 1) * Size]);
        }
    }

    void sample(float cosThetaD, float xi, float &phi, float &pdf) const
    {
        float v = (AzimuthalResolution - 1) * cosThetaD;

        int x;
        _sampler->warp(v, xi, x);

        phi = M_PI * 2.0f * (x + xi) * (1.0f / AzimuthalResolution);
        pdf = _sampler->pdf(v, x) * float(AzimuthalResolution * M_1_PI * 0.5f);
    }

    Spectrum eval(float phi, float cosThetaD) const
    {
        float u = (AzimuthalResolution - 1) * phi * M_1_PI * 0.5f;
        float v = (AzimuthalResolution - 1) * cosThetaD;
        int x0 = math::clamp(int(u), 0, AzimuthalResolution - 2);
        int y0 = math::clamp(int(v), 0, AzimuthalResolution - 2);
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        u = math::clamp(u - x0, 0.0f, 1.0f);
        v = math::clamp(v - y0, 0.0f, 1.0f);

        return (_table[x0 + y0 * AzimuthalResolution] * (1.0f - u) + _table[x1 + y0 * AzimuthalResolution] * u) * (1.0f - v) +
               (_table[x0 + y1 * AzimuthalResolution] * (1.0f - u) + _table[x1 + y1 * AzimuthalResolution] * u) * v;
    }

    float pdf(float phi, float cosThetaD) const
    {
        float u = (AzimuthalResolution - 1) * phi * M_1_PI * 0.5f;
        float v = (AzimuthalResolution - 1) * cosThetaD;
        return _sampler->pdf(v, int(u)) * float(AzimuthalResolution * M_1_PI * 0.5f);
    }

    float weight(float cosThetaD) const
    {
        float v = (AzimuthalResolution - 1) * cosThetaD;
        return _sampler->sum(v) * (M_PI * 2.0f / AzimuthalResolution);
    }
};
MTS_NAMESPACE_END
#endif