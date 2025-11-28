#include <cstdlib>
#include <iostream>
#include <string>

#include "core/color.h"
#include "core/vec3.h"
#include "io/image_io.h"
#include "render/film.h"

int main(int argc, char** argv) {
    int width = 400;
    int height = 225;
    std::string output = "gradient.ppm";

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--output" && i + 1 < argc) {
            output = argv[++i];
        } else if (arg == "--width" && i + 1 < argc) {
            width = std::atoi(argv[++i]);
        } else if (arg == "--height" && i + 1 < argc) {
            height = std::atoi(argv[++i]);
        }
    }

    if (width <= 0 || height <= 0) {
        std::cerr << "Invalid image dimensions.\n";
        return 1;
    }

    Film film(width, height);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const float u = static_cast<float>(x) / static_cast<float>(width - 1);
            const float v = static_cast<float>(y) / static_cast<float>(height - 1);

            Color color(u, v, 0.2f);
            film.add_sample(x, y, color);
        }
    }

    write_ppm(output, film);

    std::cout << "Wrote gradient image to " << output << "\n";

    return 0;
}

