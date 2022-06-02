#include "Pathtracer.h"
#include <memory>
#include <iostream>
#include <map>
#include <algorithm>
#include "material.h"
#include "embree.h"
#include "sampling.h"

using namespace std;
using namespace glm;

namespace pathtracer
{
	///////////////////////////////////////////////////////////////////////////////
	// Global variables
	///////////////////////////////////////////////////////////////////////////////
	Settings settings;
	Environment environment;
	Image rendered_image;
	PointLight point_light;

	///////////////////////////////////////////////////////////////////////////
	// Restart rendering of image
	///////////////////////////////////////////////////////////////////////////
	void restart()
	{
		// No need to clear image,
		rendered_image.number_of_samples = 0;
	}

	///////////////////////////////////////////////////////////////////////////
	// On window resize, window size is passed in, actual size of pathtraced
	// image may be smaller (if we're subsampling for speed)
	///////////////////////////////////////////////////////////////////////////
	void resize(int w, int h)
	{
		rendered_image.width = w / settings.subsampling;
		rendered_image.height = h / settings.subsampling;
		rendered_image.data.resize(rendered_image.width * rendered_image.height);
		restart();
	}

	///////////////////////////////////////////////////////////////////////////
	// Return the radiance from a certain direction wi from the environment
	// map.
	///////////////////////////////////////////////////////////////////////////
	vec3 Lenvironment(const vec3& wi)
	{
		const float theta = acos(std::max(-1.0f, std::min(1.0f, wi.y)));
		float phi = atan(wi.z, wi.x);
		if (phi < 0.0f)
			phi = phi + 2.0f * M_PI;
		vec2 lookup = vec2(phi / (2.0 * M_PI), theta / M_PI);
		return environment.multiplier * environment.map.sample(lookup.x, lookup.y);
	}

	///////////////////////////////////////////////////////////////////////////
	// Calculate the radiance going from one point (r.hitPosition()) in one
	// direction (-r.d), through path tracing.
	///////////////////////////////////////////////////////////////////////////


	vec3 Li(Ray& primary_ray)
	{
		vec3 L = vec3(0.0f);
		vec3 path_throughput = vec3(1.0);
		Ray current_ray = primary_ray;
		Ray shadow_ray;


		//if(occl){
		//shadow ray cache.

		///////////////////////////////////////////////////////////////////
		// Get the intersection information from the ray
		///////////////////////////////////////////////////////////////////
		Intersection hit = getIntersection(current_ray);
		//
		vec3 wi = normalize(point_light.position - hit.position);

		shadow_ray.o = hit.position + (EPSILON * hit.shading_normal);
		shadow_ray.d = wi;
		bool occl = occluded(shadow_ray);
		if (!occl) {
			///////////////////////////////////////////////////////////////////
			// Create a Material tree for evaluating brdfs and calculating
			// sample directions.
			///////////////////////////////////////////////////////////////////

			Diffuse diffuse(hit.material->m_color);
			BlinnPhong dielectric(hit.material->m_shininess, hit.material->m_fresnel, &diffuse);
			BlinnPhongMetal metal(hit.material->m_color, hit.material->m_shininess, hit.material->m_fresnel);
			LinearBlend metal_blend(hit.material->m_metalness, &metal, &dielectric);
			LinearBlend reflectivity_blend(hit.material->m_reflectivity, &metal_blend, &diffuse);
			BRDF& mat = dielectric;
			///////////////////////////////////////////////////////////////////
			// Calculate Direct Illumination from light.
			///////////////////////////////////////////////////////////////////
			{
				const float distance_to_light = length(point_light.position - hit.position);
				const float falloff_factor = 1.0f / (distance_to_light * distance_to_light);
				vec3 Li = point_light.intensity_multiplier * point_light.color * falloff_factor;
				L = mat.f(wi, hit.wo, hit.shading_normal) * Li * std::max(0.0f, dot(wi, hit.shading_normal));

			}
		}
		// Return the final outgoing radiance for the primary ray
		//}
		return L;
	}
	vec2 ConcentricSampleDisk(const vec2& u) {
		vec2 uOffset = 2.f * u - vec2(1, 1);
		if (uOffset.x == 0 && uOffset.y == 0)
			return vec2(0, 0);

		float
			theta, r;
		if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
			r = uOffset.x;
			theta = M_PI * (uOffset.y / uOffset.x);
		}
		else {
			r = uOffset.y;
			theta = M_PI / 2 - M_PI / 4 * (uOffset.x / uOffset.y);
		}
		return r * vec2(std::cos(theta), std::sin(theta));
	}
	struct Interaction {
		vec3 n;
		vec3 p;
	};
	class LightDisk {
	public:
		float radius;
		float height;
		vec3 position;
		vec3 normal;
		vec3 Lemit;
		LightDisk(float r, float h, vec3 p, vec3 n, vec3 lemit) : radius(r), height(h), position(p), normal(n), Lemit(lemit)
		{
		}

		float LightDisk::Area() const {
			return M_PI * radius * radius;
		}

		Interaction LightDisk::Sample(const vec2& u) const {
			vec2 pd = ConcentricSampleDisk(u);
			vec3 pObj(pd.x * radius, pd.y * radius, height);
			Interaction it;
			//(*ObjectToWorld)
			it.n = normalize(normal);
			//if (reverseOrientation) it.n *= -1;
			it.p = position + pObj;
			return it;
		}

		float LightDisk::Pdf(const Interaction& ref,
			const vec3& wi) const {

			Ray ray = Ray(ref.p, wi);
			//Float tHit;
			//SurfaceInteraction isectLight;
			if (!intersect(ray)) {
				return 0;
			}
			float a = LightDisk::Area();
			Intersection hit = getIntersection(ray);
			float pdf = pow(length(ref.p - hit.position), 2) / (a * (dot(hit.geometry_normal, -wi)));
			return pdf;
		}

		float LightDisk::Pdf_Li(const Interaction& ref,
			const vec3& wi) const {
			return Pdf(ref, wi);
		}
		//,VisibilityTester* vis
		vec3 LightDisk::Sample_Li(const Interaction& ref,
			const vec2& u, vec3* wi, float* pdf, Ray* intersect) const {
			Interaction pShape = Sample(u);
			//pShape.mediumInterface = mediumInterface;
			*wi = normalize(pShape.p - ref.p);
			*pdf = std::max(0.01f, Pdf(ref, *wi));
			//*vis = VisibilityTester(ref, pShape);

			intersect->o = ref.p;
			intersect->d = normalize(-pShape.p+ref.p);

			return L(pShape, -*wi);
		}
		vec3 L(const Interaction& intr, const vec3& w) const {
			return dot(intr.n, w) > 0.f ? Lemit : vec3(0);
		}


	};

	LightDisk disk = LightDisk(0.1, 0, vec3(3, 0, 0), vec3(0, -1, 0), vec3(0.9, 0, 0));
	LightDisk disk2 = LightDisk(0.2, 0, vec3(0, 10, 0), vec3(0, -1, 0), vec3(0.6, 0.1, 0.3));

	std::vector<LightDisk> lightData = { disk2,disk };

	vec3 LiPointLights(Ray& primary_ray) {
		vec3 L = vec3(0.0f);
		vec3 path_throughput = vec3(1.0);
		Ray current_ray = primary_ray;
		Ray shadow_ray;
		int j = 0;
		LightDisk diskit = lightData[j];
		for (int bounces = 0; bounces < settings.max_bounces; bounces++)
		{
			// Get the intersection information from the ray
			Intersection hit = getIntersection(current_ray);
			// Create a Material tree
			Diffuse diffuse(hit.material->m_color);
			BlinnPhong dielectric(hit.material->m_shininess, hit.material->m_fresnel, &diffuse);
			BlinnPhongMetal metal(hit.material->m_color, hit.material->m_shininess, hit.material->m_fresnel);
			LinearBlend metal_blend(hit.material->m_metalness, &metal, &dielectric);
			LinearBlend reflectivity_blend(hit.material->m_reflectivity, &metal_blend, &diffuse);

			BRDF& mat = reflectivity_blend;
			vec3 wi = normalize(diskit.position - hit.position);




			const float distance_to_light = length(diskit.position - hit.position);
			const float falloff_factor = 1.0f / (distance_to_light * distance_to_light);
			//sample
			Interaction i;
			i.p = hit.position + (EPSILON * hit.geometry_normal);
			i.n = hit.geometry_normal;
			Ray raySample;
			vec2 u = vec2(randf(), randf());
			vec3 sample_wi;
			float sample_pdf;
			vec3 sample_li = diskit.Sample_Li(i, u, &sample_wi, &sample_pdf, &raySample);

			vec3 Li = point_light.intensity_multiplier * diskit.Lemit * falloff_factor;

			shadow_ray.o = hit.position + (EPSILON * hit.geometry_normal);
			shadow_ray.d = sample_wi;
			bool occl = occluded(shadow_ray);
			//direct
			if (!occl) {
				L += path_throughput * mat.f(wi, hit.wo, hit.shading_normal) * Li * std::max(0.0f, dot(wi, hit.shading_normal))*10000.f;
			}
			//L += hit.material->m_emission * hit.material->m_color;

			vec3 sample_wi_brdf;
			vec3 sample_brdf;
			float sample_pdf_brdf;
			sample_brdf = mat.sample_wi(sample_wi_brdf, hit.wo, hit.shading_normal, sample_pdf_brdf);

			if (sample_pdf_brdf < 0.0001) {
				return L;
			}
			float cosineterm = abs(dot(sample_wi_brdf, hit.shading_normal));
			path_throughput = path_throughput * (sample_brdf * cosineterm) / sample_pdf_brdf;
			// If pathThroughput is zero there is no need to continue
			if (path_throughput == vec3(0)) return L;
			// Create next ray on path (existing instance can't be reused)
			//currentRay < -Create new ray instance from intersection point in outgoing direction
			//vec3 dir = normalize(hit.shading_normal - hit.position);
			current_ray = Ray(hit.position + EPSILON * hit.geometry_normal, sample_wi_brdf);
			//current_ray = raySample;

			// Bias the ray slightly to avoid self-intersection 
			if (!intersect(current_ray)) {
				return L + path_throughput * Lenvironment(current_ray.d);
			}
		}
		return L;
	}
	vec3 NewLi(Ray& primary_ray)
	{
		vec3 L = vec3(0.0f);
		vec3 path_throughput = vec3(1.0);
		Ray current_ray = primary_ray;
		Ray shadow_ray;
		for (int bounces = 0; bounces < settings.max_bounces; bounces++)
		{
			// Get the intersection information from the ray
			Intersection hit = getIntersection(current_ray);
			// Create a Material tree
			Diffuse diffuse(hit.material->m_color);
			BlinnPhong dielectric(hit.material->m_shininess, hit.material->m_fresnel, &diffuse);
			BlinnPhongMetal metal(hit.material->m_color, hit.material->m_shininess, hit.material->m_fresnel);
			LinearBlend metal_blend(hit.material->m_metalness, &metal, &dielectric);
			LinearBlend reflectivity_blend(hit.material->m_reflectivity, &metal_blend, &diffuse);

			BRDF& mat = reflectivity_blend;
			vec3 wi = normalize(point_light.position - hit.position);

			shadow_ray.o = hit.position + (EPSILON * hit.geometry_normal);
			shadow_ray.d = wi;
			bool occl = occluded(shadow_ray);
			//if (!occl) {
				// Direct illumination
			const float distance_to_light = length(point_light.position - hit.position);
			const float falloff_factor = 1.0f / (distance_to_light * distance_to_light);
			vec3 Li = point_light.intensity_multiplier * point_light.color * falloff_factor;
			if (!occl) {
				L += path_throughput * mat.f(wi, hit.wo, hit.shading_normal) * Li * std::max(0.0f, dot(wi, hit.shading_normal));
			}

			L += hit.material->m_emission * hit.material->m_color;
			//L+= pdf2 *0.01f;

			// Sample an incoming direction (and the brdf and pdf for that direction) random?
			vec3 sample_wi;
			vec3 sample_brdf;
			float sample_pdf;
			sample_brdf = mat.sample_wi(sample_wi, hit.wo, hit.shading_normal, sample_pdf);
			// Add emitted radiance from intersection
			if (sample_pdf < 0.0001) {
				return L;
			}
			float cosineterm = abs(dot(sample_wi, hit.shading_normal));
			path_throughput = path_throughput * (sample_brdf * cosineterm) / sample_pdf;

			// If pathThroughput is zero there is no need to continue
			if (path_throughput == vec3(0)) return L;
			// Create next ray on path (existing instance can't be reused)
			//currentRay < -Create new ray instance from intersection point in outgoing direction
			//vec3 dir = normalize(hit.shading_normal - hit.position);
			current_ray = Ray(hit.position + EPSILON * hit.geometry_normal, sample_wi);

			// Bias the ray slightly to avoid self-intersection 
			//+ path_throughput * Lenvironment(current_ray.d)
			if (!intersect(current_ray)) {
				return L + path_throughput * Lenvironment(current_ray.d);

			}
			//}

			// Add emitted radiance from intersection
			//L += pathThroughput * emitted light
		}
		return L;
	}

	///////////////////////////////////////////////////////////////////////////
	// Used to homogenize points transformed with projection matrices
	///////////////////////////////////////////////////////////////////////////
	inline static glm::vec3 homogenize(const glm::vec4& p)
	{
		return glm::vec3(p * (1.f / p.w));
	}

	///////////////////////////////////////////////////////////////////////////
	// Trace one path per pixel and accumulate the result in an image
	///////////////////////////////////////////////////////////////////////////
	void tracePaths(const glm::mat4& V, const glm::mat4& P)
	{
		// Stop here if we have as many samples as we want
		if ((int(rendered_image.number_of_samples) > settings.max_paths_per_pixel)
			&& (settings.max_paths_per_pixel != 0))
		{
			return;
		}
		vec3 camera_pos = vec3(glm::inverse(V) * vec4(0.0f, 0.0f, 0.0f, 1.0f));
		// Trace one path per pixel (the omp parallel stuf magically distributes the
		// pathtracing on all cores of your CPU).
		int num_rays = 0;
		vector<vec4> local_image(rendered_image.width * rendered_image.height, vec4(0.0f));

#pragma omp parallel for
		for (int y = 0; y < rendered_image.height; y++)
		{
			for (int x = 0; x < rendered_image.width; x++)
			{
				vec3 color;
				Ray primaryRay;
				primaryRay.o = camera_pos;
				// Create a ray that starts in the camera position and points toward
				// the current pixel on a virtual screen. randf() 
				vec2 screenCoord = vec2(float(x + randf()) / float(rendered_image.width),
					float(y + randf()) / float(rendered_image.height));
				// Calculate direction
				vec4 viewCoord = vec4(screenCoord.x * 2.0f - 1.0f, screenCoord.y * 2.0f - 1.0f, 1.0f, 1.0f);
				vec3 p = homogenize(inverse(P * V) * viewCoord);
				primaryRay.d = normalize(p - camera_pos);
				// Intersect ray with scene
				if (intersect(primaryRay))
				{
					// If it hit something, evaluate the radiance from that point
					color = NewLi(primaryRay);
				}
				else
				{
					// Otherwise evaluate environment
					color = Lenvironment(primaryRay.d);
				}
				// Accumulate the obtained radiance to the pixels color
				float n = float(rendered_image.number_of_samples);
				rendered_image.data[y * rendered_image.width + x] =
					rendered_image.data[y * rendered_image.width + x] * (n / (n + 1.0f))
					+ (1.0f / (n + 1.0f)) * color;
			}
		}
		rendered_image.number_of_samples += 1;
	}
}; // namespace pathtracer
