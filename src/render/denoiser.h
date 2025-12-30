#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "core/color.h"
#include "render/film.h"
#include "render/gbuffer.h"

struct DenoiseParams {
    int iterations = 5;
    float sigma_depth = 0.05f;   // relative depth
    float sigma_normal = 0.25f;  // 1 - dot(n0, n1)
    float sigma_albedo = 0.25f;  // albedo L2 distance
};

namespace denoise_detail {
inline float safe_exp_neg(float x) {
    x = std::max(0.0f, x);
    if (x > 80.0f) {
        return 0.0f;
    }
    return std::exp(-x);
}

inline float gaussian_weight(float delta, float sigma) {
    if (!(sigma > 0.0f)) {
        return 1.0f;
    }
    const float x = delta / sigma;
    return safe_exp_neg(0.5f * x * x);
}

inline float guidance_weight(std::size_t p,
                             std::size_t q,
                             const GBuffer& g,
                             const DenoiseParams& params) {
    const bool vp = g.valid_at(p) != 0;
    const bool vq = g.valid_at(q) != 0;
    if (vp != vq) {
        return 0.0f;
    }
    if (!vp) {
        return 1.0f;
    }

    const Vec3 np = g.normal_at(p);
    const Vec3 nq = g.normal_at(q);
    const float nd = std::max(0.0f, 1.0f - dot(np, nq));

    const float zp = g.depth_at(p);
    const float zq = g.depth_at(q);
    const float zden = std::max(1e-4f, std::min(zp, zq));
    const float zd = std::fabs(zp - zq) / zden;

    const Color ap = g.albedo_at(p);
    const Color aq = g.albedo_at(q);
    const Color da = ap - aq;
    const float ad = da.length();

    float w = 1.0f;
    w *= gaussian_weight(nd, params.sigma_normal);
    w *= gaussian_weight(zd, params.sigma_depth);
    w *= gaussian_weight(ad, params.sigma_albedo);
    return w;
}
}  // namespace denoise_detail

inline Film denoise_atrous(const Film& noisy, const GBuffer& gbuffer, const DenoiseParams& params) {
    const int width = noisy.width();
    const int height = noisy.height();
    if (width != gbuffer.width() || height != gbuffer.height()) {
        return noisy;
    }

    const int iterations = std::max(0, params.iterations);
    if (iterations == 0) {
        return noisy;
    }

    const std::vector<Color>& src = noisy.pixels();
    std::vector<Color> ping = src;
    std::vector<Color> pong(ping.size(), Color(0.0f));

    constexpr int kRadius = 2;
    constexpr float k1d[5] = {1.0f, 4.0f, 6.0f, 4.0f, 1.0f};

    for (int it = 0; it < iterations; ++it) {
        const int step = 1 << it;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                const std::size_t p =
                    static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                    static_cast<std::size_t>(x);

                Color sum(0.0f);
                float sum_w = 0.0f;

                for (int ky = -kRadius; ky <= kRadius; ++ky) {
                    const int sy = y + ky * step;
                    if (sy < 0 || sy >= height) continue;
                    const float wy = k1d[ky + kRadius];

                    for (int kx = -kRadius; kx <= kRadius; ++kx) {
                        const int sx = x + kx * step;
                        if (sx < 0 || sx >= width) continue;

                        const float wx = k1d[kx + kRadius];
                        const float wk = wx * wy;

                        const std::size_t q =
                            static_cast<std::size_t>(sy) * static_cast<std::size_t>(width) +
                            static_cast<std::size_t>(sx);

                        const float wg = denoise_detail::guidance_weight(p, q, gbuffer, params);
                        const float w = wk * wg;

                        sum += ping[q] * w;
                        sum_w += w;
                    }
                }

                pong[p] = (sum_w > 0.0f) ? (sum / sum_w) : ping[p];
            }
        }
        ping.swap(pong);
    }

    Film out(width, height);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const std::size_t i =
                static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                static_cast<std::size_t>(x);
            out.set_pixel(x, y, ping[i]);
        }
    }

    return out;
}

