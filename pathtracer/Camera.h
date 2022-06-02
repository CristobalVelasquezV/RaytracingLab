#pragma once
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include "embree.h"
#include <iostream> 
// Camera Declarations
namespace pathtracer {
	class Camera {
	public:
		// Camera Interface
		Camera(const glm::mat4& CameraToWorld, float shutterOpen,
			float shutterClose, Film* film, const Medium* medium);
		virtual ~Camera();
		virtual float GenerateRay(const CameraSample& sample, Ray* ray) const = 0;
		virtual float GenerateRayDifferential(const CameraSample& sample,
			RayDifferential* rd) const;
		virtual glm::vec3 We(const Ray& ray, glm::vec2* pRaster2 = nullptr) const;
		virtual void Pdf_We(const Ray& ray, float* pdfPos, float* pdfDir) const;
		//virtual Spectrum Sample_Wi(const Interaction& ref, const glm::vec2& u,
		//	Vector3f* wi, float* pdf, glm::vec2* pRaster,
		//	VisibilityTester* vis) const;

		// Camera Public Data
		glm::mat4 CameraToWorld;
		const float shutterOpen, shutterClose;
		//Film* film;

	};

	struct CameraSample {
		glm::vec2 pFilm;
		glm::vec2 pLens;
		float time;
	};



	class ProjectiveCamera : public Camera {
	public:
		// ProjectiveCamera Public Methods
		ProjectiveCamera(const glm::mat4& CameraToWorld,
			const glm::mat4& CameraToScreen,
			const Bounds2f& screenWindow, float shutterOpen,
			float shutterClose, float lensr, float focald, Film* film,
			const Medium* medium)
			: Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
			CameraToScreen(CameraToScreen) {
			// Initialize depth of field parameters
			lensRadius = lensr;
			focalDistance = focald;

			// Compute projective camera glm::mat4ations

			// Compute projective camera screen glm::mat4ations
			ScreenToRaster =
				scale(glm::vec3(film->fullResolution.x, film->fullResolution.y, 1)) *
				scale(glm::vec3(1 / (screenWindow.pMax.x - screenWindow.pMin.x),
					1 / (screenWindow.pMin.y - screenWindow.pMax.y), 1)) *
				translate(glm::vec3(-screenWindow.pMin.x, -screenWindow.pMax.y, 0));
			RasterToScreen = inverse(ScreenToRaster);
			RasterToCamera = inverse(CameraToScreen) * RasterToScreen;
		}

	protected:
		// ProjectiveCamera Protected Data
		glm::mat4 CameraToScreen, RasterToCamera;
		glm::mat4 ScreenToRaster, RasterToScreen;
		float lensRadius, focalDistance;
	};

}
