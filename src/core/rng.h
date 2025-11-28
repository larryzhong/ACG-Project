#pragma once

#include <cstdint>

// Simple linear congruential generator for reproducible randomness.
struct RNG {
    std::uint64_t state;

    RNG() : state(1u) {}
    explicit RNG(std::uint64_t seed) : state(seed ? seed : 1u) {}

    // Returns a float in [0, 1).
    float uniform() {
        // Parameters from Numerical Recipes.
        state = state * 6364136223846793005ull + 1u;
        const std::uint32_t bits = static_cast<std::uint32_t>(state >> 32);
        return static_cast<float>(bits) / static_cast<float>(0xFFFFFFFFu);
    }
};

