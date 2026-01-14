
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "core/rng.h"
#include "core/ray.h"
#include "core/vec3.h"
#include "scene/bvh.h"
#include "scene/hittable.h"
#include "scene/sphere.h"

namespace {

// -----------------------------
// Median-split BVH (baseline)
// -----------------------------

class MedianBVHNode : public Hittable {
public:
	MedianBVHNode() = default;

	MedianBVHNode(std::vector<HittablePtr>& objects, std::size_t start, std::size_t end) {
		const std::size_t object_span = end - start;
		if (object_span == 0) {
			return;
		}

		constexpr std::size_t kLeafSize = 4;
		if (object_span <= kLeafSize) {
			is_leaf_ = true;
			leaf_objects_.reserve(object_span);
			bool first = true;
			for (std::size_t i = start; i < end; ++i) {
				if (!objects[i]) continue;
				leaf_objects_.push_back(objects[i]);
				const AABB b = objects[i]->bounding_box();
				if (first) {
					box_ = b;
					first = false;
				} else {
					box_ = surrounding_box(box_, b);
				}
			}
			return;
		}

		const int axis = choose_split_axis(objects, start, end);
		const std::size_t mid = start + object_span / 2;

		auto comparator = [axis](const HittablePtr& a, const HittablePtr& b) {
			if (!a) return true;
			if (!b) return false;
			return centroid_axis(a->bounding_box(), axis) < centroid_axis(b->bounding_box(), axis);
		};

		std::nth_element(objects.begin() + static_cast<std::ptrdiff_t>(start),
						 objects.begin() + static_cast<std::ptrdiff_t>(mid),
						 objects.begin() + static_cast<std::ptrdiff_t>(end),
						 comparator);

		left_ = std::make_shared<MedianBVHNode>(objects, start, mid);
		right_ = std::make_shared<MedianBVHNode>(objects, mid, end);

		const AABB box_left = left_ ? left_->bounding_box() : AABB();
		const AABB box_right = right_ ? right_->bounding_box() : AABB();
		box_ = surrounding_box(box_left, box_right);
	}

	bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const override {
		if (!box_.hit(r, t_min, t_max)) {
			return false;
		}

		if (is_leaf_) {
			bool hit_anything = false;
			float closest_so_far = t_max;
			HitRecord temp;
			for (const auto& obj : leaf_objects_) {
				if (obj && obj->hit(r, t_min, closest_so_far, temp)) {
					hit_anything = true;
					closest_so_far = temp.t;
					rec = temp;
				}
			}
			return hit_anything;
		}

		HitRecord left_rec;
		HitRecord right_rec;
		const bool hit_left = left_ && left_->hit(r, t_min, t_max, left_rec);
		const bool hit_right = right_ && right_->hit(r, t_min, t_max, right_rec);

		if (hit_left && hit_right) {
			rec = (left_rec.t < right_rec.t) ? left_rec : right_rec;
			return true;
		}
		if (hit_left) {
			rec = left_rec;
			return true;
		}
		if (hit_right) {
			rec = right_rec;
			return true;
		}
		return false;
	}

	AABB bounding_box() const override {
		return box_;
	}

private:
	static float centroid_axis(const AABB& b, int axis) {
		if (axis == 0) return 0.5f * (b.min.x + b.max.x);
		if (axis == 1) return 0.5f * (b.min.y + b.max.y);
		return 0.5f * (b.min.z + b.max.z);
	}

	static int choose_split_axis(const std::vector<HittablePtr>& objects,
								 std::size_t start,
								 std::size_t end) {
		AABB bounds;
		bool first_box = true;
		for (std::size_t i = start; i < end; ++i) {
			if (!objects[i]) continue;
			const AABB box = objects[i]->bounding_box();
			if (first_box) {
				bounds = box;
				first_box = false;
			} else {
				bounds = surrounding_box(bounds, box);
			}
		}

		const float dx = bounds.max.x - bounds.min.x;
		const float dy = bounds.max.y - bounds.min.y;
		const float dz = bounds.max.z - bounds.min.z;
		if (dx > dy && dx > dz) return 0;
		if (dy > dz) return 1;
		return 2;
	}

	HittablePtr left_;
	HittablePtr right_;
	AABB box_;
	bool is_leaf_ = false;
	std::vector<HittablePtr> leaf_objects_;
};

float lerp(float a, float b, float t) {
	return a + (b - a) * t;
}

float rand_range(RNG& rng, float a, float b) {
	return lerp(a, b, rng.uniform());
}

Vec3 test_random_in_unit_sphere(RNG& rng) {
	// Rejection sampling.
	for (int i = 0; i < 64; ++i) {
		const Vec3 p(
			rand_range(rng, -1.0f, 1.0f),
			rand_range(rng, -1.0f, 1.0f),
			rand_range(rng, -1.0f, 1.0f));
		if (p.length_squared() < 1.0f) {
			return p;
		}
	}
	return Vec3(1.0f, 0.0f, 0.0f);
}

Vec3 test_random_unit_vector(RNG& rng) {
	const Vec3 p = test_random_in_unit_sphere(rng);
	const float len2 = p.length_squared();
	if (len2 <= 1e-20f) {
		return Vec3(0.0f, 0.0f, 1.0f);
	}
	return p / std::sqrt(len2);
}

struct HitResult {
	bool hit = false;
	HitRecord rec{};
};

HitResult hit_bruteforce(const std::vector<HittablePtr>& objects,
						 const Ray& r,
						 float t_min,
						 float t_max) {
	HitResult out;
	float closest_so_far = t_max;
	HitRecord tmp;

	for (const auto& obj : objects) {
		if (!obj) {
			continue;
		}
		if (obj->hit(r, t_min, closest_so_far, tmp)) {
			out.hit = true;
			closest_so_far = tmp.t;
			out.rec = tmp;
		}
	}

	return out;
}

std::shared_ptr<BVHNode> build_sah_bvh_once(const std::vector<HittablePtr>& objects) {
	if (objects.empty()) {
		return {};
	}
	std::vector<HittablePtr> objects_copy = objects;
	return std::make_shared<BVHNode>(objects_copy, 0, objects_copy.size());
}

std::shared_ptr<MedianBVHNode> build_median_bvh_once(const std::vector<HittablePtr>& objects) {
	if (objects.empty()) {
		return {};
	}
	std::vector<HittablePtr> objects_copy = objects;
	return std::make_shared<MedianBVHNode>(objects_copy, 0, objects_copy.size());
}

HitResult hit_accel(const Hittable& accel, const Ray& r, float t_min, float t_max) {
	HitResult out;
	out.hit = accel.hit(r, t_min, t_max, out.rec);
	return out;
}

bool nearly_equal(float a, float b, float rel = 1e-4f, float abs = 1e-5f) {
	const float diff = std::fabs(a - b);
	if (diff <= abs) {
		return true;
	}
	const float scale = std::max(1.0f, std::max(std::fabs(a), std::fabs(b)));
	return diff <= rel * scale;
}

bool vec_nearly_equal(const Vec3& a, const Vec3& b, float eps = 1e-3f) {
	return nearly_equal(a.x, b.x, 1e-4f, eps) &&
		   nearly_equal(a.y, b.y, 1e-4f, eps) &&
		   nearly_equal(a.z, b.z, 1e-4f, eps);
}

std::vector<HittablePtr> make_uniform_spheres(RNG& rng, std::size_t count) {
	std::vector<HittablePtr> objects;
	objects.reserve(count);

	for (std::size_t i = 0; i < count; ++i) {
		const Vec3 center(
			rand_range(rng, -5.0f, 5.0f),
			rand_range(rng, -2.0f, 2.0f),
			rand_range(rng, -5.0f, 5.0f));
		const float radius = rand_range(rng, 0.05f, 0.6f);
		objects.push_back(std::make_shared<Sphere>(center, radius, MaterialPtr()));
	}

	return objects;
}

std::vector<HittablePtr> make_clustered_spheres(RNG& rng, std::size_t count) {
	std::vector<HittablePtr> objects;
	objects.reserve(count);

	const Vec3 cluster_center(
		rand_range(rng, -2.0f, 2.0f),
		rand_range(rng, -1.0f, 1.0f),
		rand_range(rng, -2.0f, 2.0f));

	for (std::size_t i = 0; i < count; ++i) {
		const bool in_cluster = (rng.uniform() < 0.9f);
		const Vec3 center = in_cluster
			? cluster_center + Vec3(rand_range(rng, -0.5f, 0.5f),
									rand_range(rng, -0.5f, 0.5f),
									rand_range(rng, -0.5f, 0.5f))
			: Vec3(rand_range(rng, -25.0f, 25.0f),
				   rand_range(rng, -10.0f, 10.0f),
				   rand_range(rng, -25.0f, 25.0f));

		const float radius = in_cluster ? rand_range(rng, 0.05f, 0.3f)
										: rand_range(rng, 0.05f, 1.0f);
		objects.push_back(std::make_shared<Sphere>(center, radius, MaterialPtr()));
	}

	return objects;
}

Ray make_random_ray(RNG& rng) {
	const Vec3 origin(
		rand_range(rng, -15.0f, 15.0f),
		rand_range(rng, -10.0f, 10.0f),
		rand_range(rng, -15.0f, 15.0f));
	const Vec3 dir = test_random_unit_vector(rng);
	return Ray(origin, dir, 0.0f);
}

void add_axis_rays(std::vector<Ray>& rays, const Vec3& origin) {
	// Stress AABB slab code for zero components.
	rays.emplace_back(origin, Vec3(1.0f, 0.0f, 0.0f));
	rays.emplace_back(origin, Vec3(-1.0f, 0.0f, 0.0f));
	rays.emplace_back(origin, Vec3(0.0f, 1.0f, 0.0f));
	rays.emplace_back(origin, Vec3(0.0f, -1.0f, 0.0f));
	rays.emplace_back(origin, Vec3(0.0f, 0.0f, 1.0f));
	rays.emplace_back(origin, Vec3(0.0f, 0.0f, -1.0f));
	rays.emplace_back(origin, Vec3(1.0f, 1.0f, 0.0f));
	rays.emplace_back(origin, Vec3(1.0f, 0.0f, 1.0f));
	rays.emplace_back(origin, Vec3(0.0f, 1.0f, 1.0f));
}

bool check_one_scene(const char* name,
					 const std::vector<HittablePtr>& objects,
					 std::uint64_t seed,
					 std::size_t ray_count) {
	const float t_min = 1e-4f;
	const float t_max = 1e30f;
	RNG rng(seed);

	std::vector<Ray> rays;
	rays.reserve(ray_count + 32);
	for (std::size_t i = 0; i < ray_count; ++i) {
		rays.push_back(make_random_ray(rng));
	}
	add_axis_rays(rays, Vec3(0.0f, 0.0f, 0.0f));

	const auto sah_bvh_ptr = build_sah_bvh_once(objects);
	const auto median_bvh_ptr = build_median_bvh_once(objects);
	if (!sah_bvh_ptr || !median_bvh_ptr) {
		std::cerr << "[BVH TEST] Scene '" << name << "' has no objects.\n";
		return true;
	}

	std::size_t mismatches = 0;
	for (std::size_t i = 0; i < rays.size(); ++i) {
		const Ray& r = rays[i];
		const HitResult brute = hit_bruteforce(objects, r, t_min, t_max);
		const HitResult bvh_median = hit_accel(*median_bvh_ptr, r, t_min, t_max);
		const HitResult bvh_sah = hit_accel(*sah_bvh_ptr, r, t_min, t_max);

		auto mismatch = [&](const HitResult& candidate) {
			if (brute.hit != candidate.hit) return true;
			if (!brute.hit) return false;
			if (!std::isfinite(brute.rec.t) || !std::isfinite(candidate.rec.t)) return true;
			if (!nearly_equal(brute.rec.t, candidate.rec.t)) return true;
			const Vec3 p_brute = r.at(brute.rec.t);
			const Vec3 p_cand = r.at(candidate.rec.t);
			return !vec_nearly_equal(p_brute, p_cand);
		};

		if (mismatch(bvh_median) || mismatch(bvh_sah)) {
			++mismatches;
		}

		if (mismatches > 0) {
			std::cerr << "[BVH TEST] Mismatch in scene '" << name << "' at ray " << i << "\n";
			std::cerr << "  origin = (" << r.origin.x << ", " << r.origin.y << ", " << r.origin.z << ")\n";
			std::cerr << "  dir    = (" << r.direction.x << ", " << r.direction.y << ", " << r.direction.z << ")\n";
			std::cerr << "  brute.hit = " << brute.hit;
			if (brute.hit) std::cerr << ", t = " << brute.rec.t;
			std::cerr << "\n";
			std::cerr << "  bvh(median).hit = " << bvh_median.hit;
			if (bvh_median.hit) std::cerr << ", t = " << bvh_median.rec.t;
			std::cerr << "\n";
			std::cerr << "  bvh(sah).hit    = " << bvh_sah.hit;
			if (bvh_sah.hit) std::cerr << ", t = " << bvh_sah.rec.t;
			std::cerr << "\n";
			return false;
		}
	}

	std::cout << "[BVH TEST] OK: " << name << " (objects=" << objects.size()
			  << ", rays=" << rays.size() << ")\n";
	return true;
}

struct BenchResult {
	double build_ms = 0.0;
	double brute_ms = 0.0;
	double bvh_ms = 0.0;
	double sah_build_ms = 0.0;
	double sah_ms = 0.0;
	std::size_t brute_hits = 0;
	std::size_t bvh_hits = 0;
	std::size_t sah_hits = 0;
	std::uint64_t checksum = 0;
};

std::uint32_t float_bits(float v) {
	std::uint32_t u = 0;
	static_assert(sizeof(float) == sizeof(std::uint32_t), "unexpected float size");
	std::memcpy(&u, &v, sizeof(u));
	return u;
}

BenchResult benchmark_one_scene(const char* name,
								const std::vector<HittablePtr>& objects,
								std::uint64_t seed,
								std::size_t ray_count,
								bool include_bruteforce) {
	using clock = std::chrono::steady_clock;

	const float t_min = 1e-4f;
	const float t_max = 1e30f;

	RNG rng(seed);
	std::vector<Ray> rays;
	rays.reserve(ray_count);
	for (std::size_t i = 0; i < ray_count; ++i) {
		rays.push_back(make_random_ray(rng));
	}

	BenchResult res;

	const auto t0 = clock::now();
	const auto median_bvh_ptr = build_median_bvh_once(objects);
	const auto t1 = clock::now();
	res.build_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

	const auto t1b = clock::now();
	const auto sah_bvh_ptr = build_sah_bvh_once(objects);
	const auto t1c = clock::now();
	res.sah_build_ms = std::chrono::duration<double, std::milli>(t1c - t1b).count();

	if (!median_bvh_ptr || !sah_bvh_ptr) {
		std::cout << "[BVH BENCH] " << name << ": empty scene\n";
		return res;
	}

	std::uint64_t checksum = 0;
	std::size_t hits = 0;
	if (include_bruteforce) {
		// Brute-force trace
		const auto t2 = clock::now();
		for (const auto& r : rays) {
			const HitResult brute = hit_bruteforce(objects, r, t_min, t_max);
			if (brute.hit) {
				++hits;
				checksum ^= (static_cast<std::uint64_t>(float_bits(brute.rec.t)) << 1);
			}
		}
		const auto t3 = clock::now();
		res.brute_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();
		res.brute_hits = hits;
		res.checksum ^= checksum;
	} else {
		res.brute_ms = 0.0;
		res.brute_hits = 0;
	}

	// BVH trace (median baseline)
	checksum = 0;
	hits = 0;
	const auto t4 = clock::now();
	for (const auto& r : rays) {
		const HitResult bvh = hit_accel(*median_bvh_ptr, r, t_min, t_max);
		if (bvh.hit) {
			++hits;
			checksum ^= (static_cast<std::uint64_t>(float_bits(bvh.rec.t)) << 1);
		}
	}
	const auto t5 = clock::now();
	res.bvh_ms = std::chrono::duration<double, std::milli>(t5 - t4).count();
	res.bvh_hits = hits;
	res.checksum ^= (checksum << 32);

	// BVH trace (SAH)
	checksum = 0;
	hits = 0;
	const auto t6 = clock::now();
	for (const auto& r : rays) {
		const HitResult sah = hit_accel(*sah_bvh_ptr, r, t_min, t_max);
		if (sah.hit) {
			++hits;
			checksum ^= (static_cast<std::uint64_t>(float_bits(sah.rec.t)) << 1);
		}
	}
	const auto t7 = clock::now();
	res.sah_ms = std::chrono::duration<double, std::milli>(t7 - t6).count();
	res.sah_hits = hits;
	res.checksum ^= (checksum << 48);

	std::cout << "[BVH BENCH] " << name
			  << " | objects=" << objects.size()
			  << " | rays=" << rays.size()
			  << " | build(bvh)=" << res.build_ms << "ms"
			  << " | build(sah)=" << res.sah_build_ms << "ms";

	if (include_bruteforce) {
		std::cout << " | brute=" << res.brute_ms << "ms"
				  << " | bvh=" << res.bvh_ms << "ms"
				  << " | sah=" << res.sah_ms << "ms"
				  << " | speedup(brute/bvh)=" << (res.bvh_ms > 0.0 ? (res.brute_ms / res.bvh_ms) : 0.0)
				  << " | speedup(brute/sah)=" << (res.sah_ms > 0.0 ? (res.brute_ms / res.sah_ms) : 0.0)
				  << " | speedup(bvh/sah)=" << (res.sah_ms > 0.0 ? (res.bvh_ms / res.sah_ms) : 0.0)
				  << " | hits(brute/bvh/sah)=" << res.brute_hits << "/" << res.bvh_hits << "/" << res.sah_hits;
	} else {
		std::cout << " | brute=SKIP(no-brute)"
				  << " | bvh=" << res.bvh_ms << "ms"
				  << " | sah=" << res.sah_ms << "ms"
				  << " | speedup(bvh/sah)=" << (res.sah_ms > 0.0 ? (res.bvh_ms / res.sah_ms) : 0.0)
				  << " | hits(bvh/sah)=" << res.bvh_hits << "/" << res.sah_hits;
	}

	std::cout << " | checksum=" << res.checksum << "\n";

	return res;
}

bool arg_is_flag(const char* s) {
	return s && s[0] == '-' && s[1] == '-';
}

} // namespace

int main(int argc, char** argv) {
	std::uint64_t seed = 1337u;
	std::size_t ray_count = 200000;
	std::size_t object_count = 5000;

	bool run_correctness = true;
	bool run_speed = false;
	bool speed_include_bruteforce = true;

	int numeric_index = 0;
	for (int i = 1; i < argc; ++i) {
		const std::string arg = argv[i] ? argv[i] : "";
		if (arg == "--speed") {
			run_speed = true;
		} else if (arg == "--perf" || arg == "--bench-accel") {
			// Convenience: only compare BVH vs SAH performance.
			run_speed = true;
			run_correctness = false;
			speed_include_bruteforce = false;
		} else if (arg == "--no-brute" || arg == "--accel-only") {
			speed_include_bruteforce = false;
		} else if (arg == "--correctness") {
			run_correctness = true;
		} else if (arg == "--no-correctness") {
			run_correctness = false;
		} else {
			// Positional numerics: seed, ray_count, object_count
			if (numeric_index == 0) {
				seed = static_cast<std::uint64_t>(std::strtoull(argv[i], nullptr, 10));
			} else if (numeric_index == 1) {
				ray_count = static_cast<std::size_t>(std::strtoull(argv[i], nullptr, 10));
			} else if (numeric_index == 2) {
				object_count = static_cast<std::size_t>(std::strtoull(argv[i], nullptr, 10));
			}
			++numeric_index;
		}
	}

	RNG rng(seed);

	const auto uniform = make_uniform_spheres(rng, object_count);
	const auto clustered = make_clustered_spheres(rng, object_count);

	if (run_correctness) {
		if (!check_one_scene("uniform_spheres", uniform, seed ^ 0x1234u, ray_count)) {
			return 1;
		}

		if (!check_one_scene("clustered_spheres", clustered, seed ^ 0xBEEF'CAFEu, ray_count)) {
			return 1;
		}

		std::cout << "[BVH TEST] All tests passed.\n";
	}

	if (run_speed) {
		(void)benchmark_one_scene("uniform_spheres", uniform, seed ^ 0x1234u, ray_count, speed_include_bruteforce);
		(void)benchmark_one_scene("clustered_spheres", clustered, seed ^ 0xBEEF'CAFEu, ray_count, speed_include_bruteforce);
	}
	return 0;
}
