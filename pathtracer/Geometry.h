#pragma once
#include <glm/glm.hpp>
namespace pathtracer {
	// Bounds Declarations
	template <typename T>
	class Bounds2 {
	public:
		// Bounds2 Public Methods
		Bounds2() {
			T minNum = std::numeric_limits<T>::lowest();
			T maxNum = std::numeric_limits<T>::max();
			pMin = glm::vec<2,T,highp>(maxNum, maxNum);
			pMax = glm::vec<2,T,highp>(minNum, minNum);
		}
		explicit Bounds2(const glm::vec<2,T,highp>& p) : pMin(p), pMax(p) {}
		Bounds2(const glm::vec<2,T,highp>& p1, const glm::vec<2,T,highp>& p2) {
			pMin = glm::vec<2,T,highp>(std::min(p1.x, p2.x), std::min(p1.y, p2.y));
			pMax = glm::vec<2,T,highp>(std::max(p1.x, p2.x), std::max(p1.y, p2.y));
		}
		template <typename U>
		explicit operator Bounds2<U>() const {
			return Bounds2<U>((glm::vec<2,U,highp>)pMin, (glm::vec<2,U,highp>)pMax);
		}

		glm::vec<2,U,highp> Diagonal() const { return pMax - pMin; }
		T Area() const {
			glm::vec<2,U,highp> d = pMax - pMin;
			return (d.x * d.y);
		}
		int MaximumExtent() const {
			glm::vec<2,U,highp> diag = Diagonal();
			if (diag.x > diag.y)
				return 0;
			else
				return 1;
		}
		inline const glm::vec<2,T,highp>& operator[](int i) const {
			DCHECK(i == 0 || i == 1);
			return (i == 0) ? pMin : pMax;
		}
		inline glm::vec<2,T,highp>& operator[](int i) {
			DCHECK(i == 0 || i == 1);
			return (i == 0) ? pMin : pMax;
		}
		bool operator==(const Bounds2<T>& b) const {
			return b.pMin == pMin && b.pMax == pMax;
		}
		bool operator!=(const Bounds2<T>& b) const {
			return b.pMin != pMin || b.pMax != pMax;
		}
		glm::vec<2,T,highp> Lerp(const Point2f& t) const {
			return glm::vec<2,T,highp>(pbrt::Lerp(t.x, pMin.x, pMax.x),
				pbrt::Lerp(t.y, pMin.y, pMax.y));
		}
		glm::vec<2,U,highp> Offset(const glm::vec<2,T,highp>& p) const {
			glm::vec<2,U,highp> o = p - pMin;
			if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
			if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
			return o;
		}
		void BoundingSphere(glm::vec<2,T,highp>* c, Float* rad) const {
			*c = (pMin + pMax) / 2;
			*rad = Inside(*c, *this) ? Distance(*c, pMax) : 0;
		}
		friend std::ostream& operator<<(std::ostream& os, const Bounds2<T>& b) {
			os << "[ " << b.pMin << " - " << b.pMax << " ]";
			return os;
		}

		// Bounds2 Public Data
		glm::vec<2,T,highp> pMin, pMax;
	};

	template <typename T>
	class Bounds3 {
	public:
		// Bounds3 Public Methods
		Bounds3() {
			T minNum = std::numeric_limits<T>::lowest();
			T maxNum = std::numeric_limits<T>::max();
			pMin = glm::vec<3,U,highp>(maxNum, maxNum, maxNum);
			pMax = glm::vec<3,U,highp>(minNum, minNum, minNum);
		}
		explicit Bounds3(const glm::vec<3,U,highp>& p) : pMin(p), pMax(p) {}
		Bounds3(const glm::vec<3,U,highp>& p1, const glm::vec<3,U,highp>& p2)
			: pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
				std::min(p1.z, p2.z)),
			pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
				std::max(p1.z, p2.z)) {}
		const glm::vec<3,U,highp>& operator[](int i) const;
		glm::vec<3,U,highp>& operator[](int i);
		bool operator==(const Bounds3<T>& b) const {
			return b.pMin == pMin && b.pMax == pMax;
		}
		bool operator!=(const Bounds3<T>& b) const {
			return b.pMin != pMin || b.pMax != pMax;
		}
		glm::vec<3,U,highp> Corner(int corner) const {
			DCHECK(corner >= 0 && corner < 8);
			return glm::vec<3,U,highp>((*this)[(corner & 1)].x,
				(*this)[(corner & 2) ? 1 : 0].y,
				(*this)[(corner & 4) ? 1 : 0].z);
		}
		glm::vec<3,U,highp> Diagonal() const { return pMax - pMin; }
		T SurfaceArea() const {
			glm::vec<3,U,highp> d = Diagonal();
			return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
		}
		T Volume() const {
			glm::vec<3,U,highp> d = Diagonal();
			return d.x * d.y * d.z;
		}
		int MaximumExtent() const {
			glm::vec<3,U,highp> d = Diagonal();
			if (d.x > d.y && d.x > d.z)
				return 0;
			else if (d.y > d.z)
				return 1;
			else
				return 2;
		}
		glm::vec<3,U,highp> Lerp(const Point3f& t) const {
			return glm::vec<3,U,highp>(pbrt::Lerp(t.x, pMin.x, pMax.x),
				pbrt::Lerp(t.y, pMin.y, pMax.y),
				pbrt::Lerp(t.z, pMin.z, pMax.z));
		}
		glm::vec<3,U,highp> Offset(const glm::vec<3,U,highp>& p) const {
			glm::vec<3,U,highp> o = p - pMin;
			if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
			if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
			if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
			return o;
		}
		void BoundingSphere(glm::vec<3,U,highp>* center, Float* radius) const {
			*center = (pMin + pMax) / 2;
			*radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
		}
		template <typename U>
		explicit operator Bounds3<U>() const {
			return Bounds3<U>((Point3<U>)pMin, (Point3<U>)pMax);
		}
		//bool IntersectP(const Ray& ray, Float* hitt0 = nullptr,
		//	Float* hitt1 = nullptr) const;
		//inline bool IntersectP(const Ray& ray, const Vector3f& invDir,
		//	const int dirIsNeg[3]) const;
		friend std::ostream& operator<<(std::ostream& os, const Bounds3<T>& b) {
			os << "[ " << b.pMin << " - " << b.pMax << " ]";
			return os;
		}

		// Bounds3 Public Data
		glm::vec<3,U,highp> pMin, pMax;
	};

	typedef Bounds2<float> Bounds2f;
	typedef Bounds2<int> Bounds2i;
	typedef Bounds3<float> Bounds3f;
	typedef Bounds3<int> Bounds3i;
}