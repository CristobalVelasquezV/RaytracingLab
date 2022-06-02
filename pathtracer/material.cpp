#include "material.h"
#include "sampling.h"

namespace pathtracer
{
	///////////////////////////////////////////////////////////////////////////
	// A Lambertian (diffuse) material
	///////////////////////////////////////////////////////////////////////////
	vec3 Diffuse::f(const vec3& wi, const vec3& wo, const vec3& n)
	{
		if (dot(wi, n) <= 0.0f)
			return vec3(0.0f);
		if (!sameHemisphere(wi, wo, n))
			return vec3(0.0f);
		return (1.0f / M_PI) * color;
	}

	vec3 Diffuse::sample_wi(vec3& wi, const vec3& wo, const vec3& n, float& p)
	{
		vec3 tangent = normalize(perpendicular(n));
		vec3 bitangent = normalize(cross(tangent, n));
		vec3 sample = cosineSampleHemisphere();
		wi = normalize(sample.x * tangent + sample.y * bitangent + sample.z * n);
		if (dot(wi, n) <= 0.0f)
			p = 0.0f;
		else
			p = max(0.0f, dot(n, wi)) / M_PI;
		return f(wi, wo, n);
	}

	///////////////////////////////////////////////////////////////////////////
	// A Blinn Phong Dielectric Microfacet BRFD
	///////////////////////////////////////////////////////////////////////////
	vec3 BlinnPhong::refraction_brdf(const vec3& wi, const vec3& wo, const vec3& n)
	{
		float dotnwi = dot(n, wi);
		if (dotnwi <= 0) {
			return vec3(0);
		}
		vec3 wh = normalize(wi + wo);
		float dotnwh = max(0.01f, dot(n, wh));
		float fwi = R0 + ((1.0f - R0) * pow((1.0f - dot(wh, wi)), 5.0f));
		float refraction_brdf = 1 - fwi;
		if (refraction_layer != NULL) {
			return refraction_brdf * refraction_layer->f(wi, wo, n);
		}

		return vec3(0);
	}
	vec3 BlinnPhong::reflection_brdf(const vec3& wi, const vec3& wo, const vec3& n)
	{
		float dotnwi = dot(n, wi);
		if (dotnwi <= 0) {
			return vec3(0, 0, 0);
		}
		vec3 wh = normalize(wi + wo);
		float dotnwh = max(0.01f, dot(n, wh));
		float fwi = R0 + ((1.0f - R0) * pow((1.0f - dot(wh, wi)), 5.0f));
		float dwh = ((shininess + 2.0f) / (M_PI * 2.0f)) * pow(dotnwh, shininess);
		float dotnwo = max(0.01f, dot(n, wo));
		float dotwowh = max(0.01f, dot(wo, wh));
		float constant = (2 * dotnwh) / dotwowh;
		float gwiwo = min(1.f, min((constant * dotnwo), (constant * dotnwi)));
		float reflection_brdf = (fwi * dwh * gwiwo) / (4 * dotnwi * dotnwo);

		return vec3(reflection_brdf);
	}

	vec3 BlinnPhong::f(const vec3& wi, const vec3& wo, const vec3& n)
	{
		return reflection_brdf(wi, wo, n) + refraction_brdf(wi, wo, n);
	}

	vec3 BlinnPhong::sample_wi(vec3& wi, const vec3& wo, const vec3& n, float& p)
	{
		//vec3 tangent = normalize(perpendicular(n));
		//vec3 bitangent = normalize(cross(tangent, n));
		//vec3 sample = cosineSampleHemisphere();
		//wi = normalize(sample.x * tangent + sample.y * bitangent + sample.z * n);
		//if(dot(wi, n) <= 0.0f)
		//	p = 0.0f;
		//else
		//	p = max(0.0f, dot(n, wi)) / M_PI;
		//return f(wi, wo, n);
		vec3 tangent = normalize(perpendicular(n));
		vec3 bitangent = normalize(cross(tangent, n));
		float phi = 2.0f * M_PI * randf();
		float cos_theta = pow(randf(), 1.0f / (shininess + 1));
		float sin_theta = sqrt(max(0.0f, 1.0f - cos_theta * cos_theta));
		vec3 wh = normalize(sin_theta * cos(phi) * tangent +
			sin_theta * sin(phi) * bitangent +
			cos_theta * n);


		if (randf() < 0.5) {
			// Sample a direction based on the Microfacet brdf

			float pwh = ((shininess + 1) * pow(dot(n, wh), shininess)) / (2 * M_PI);
			p = pwh / (4 * dot(wo, wh));
			wi = -reflect(wo, wh);
			if (dot(wo, n) <= 0.0f) return vec3(0.0f);
			p = p * 0.5;
			return reflection_brdf(wi, wo, n);
		}
		else {
			// Sample a direction for the underlying layer
			if (refraction_layer != NULL) {
				vec3 brdf = refraction_layer->sample_wi(wi, wo, n, p);
				p = p * 0.5;
				// We need to attenuate the refracted brdf with (1 - F)
				float F = R0 + (1.0f - R0) * pow(1.0f - abs(dot(wh, wi)), 5.0f);
				return (1 - F) * brdf;
			}
			else
			{
				return vec3(0);
			}

		}

	}

	///////////////////////////////////////////////////////////////////////////
	// A Blinn Phong Metal Microfacet BRFD (extends the BlinnPhong class)
	///////////////////////////////////////////////////////////////////////////

	vec3 BlinnPhongMetal::refraction_brdf(const vec3& wi, const vec3& wo, const vec3& n)
	{
		return vec3(0.0f);
	}
	vec3 BlinnPhongMetal::reflection_brdf(const vec3& wi, const vec3& wo, const vec3& n)
	{
		return BlinnPhong::reflection_brdf(wi, wo, n) * color;
	};

	///////////////////////////////////////////////////////////////////////////
	// A Linear Blend between two BRDFs
	///////////////////////////////////////////////////////////////////////////

	vec3 LinearBlend::f(const vec3& wi, const vec3& wo, const vec3& n)
	{
		vec3 f0 = bsdf0->f(wi, wo, n);
		vec3 f1 = bsdf1->f(wi, wo, n);

		return  f0 * w + (f1 * (1 - w));
	}

	vec3 LinearBlend::sample_wi(vec3& wi, const vec3& wo, const vec3& n, float& p)
	{
		if (randf() < w) {
			return bsdf0->sample_wi(wi, wo, n, p);
		}
		else {
			return bsdf1->sample_wi(wi, wo, n, p);
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// A perfect specular refraction.
	///////////////////////////////////////////////////////////////////////////
} // namespace pathtracer