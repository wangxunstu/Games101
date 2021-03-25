//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"
#include "global.hpp"
enum MaterialType { DIFFUSE,REFR,SPEC};

// Reflection Declarations
float FrDielectric(float cosThetaI, float etaI, float etaT);
Vector3f FrConductor(float cosThetaI, const Vector3f &etaI,
	const Vector3f &etaT, const Vector3f &k);

// BSDF Inline Functions
inline float CosTheta(const Vector3f &w) { return w.z; }
inline float Cos2Theta(const Vector3f &w) { return w.z * w.z; }
inline float AbsCosTheta(const Vector3f &w) { return std::abs(w.z); }
inline float Sin2Theta(const Vector3f &w) {
	return std::max((float)0, (float)1 - Cos2Theta(w));
}

inline float SinTheta(const Vector3f &w) { return std::sqrt(Sin2Theta(w)); }

inline float TanTheta(const Vector3f &w) { return SinTheta(w) / CosTheta(w); }

inline float Tan2Theta(const Vector3f &w) {
	return Sin2Theta(w) / Cos2Theta(w);
}

inline float CosPhi(const Vector3f &w) {
	float sinTheta = SinTheta(w);
	return (sinTheta == 0) ? 1 : clamp( -1, 1, w.x / sinTheta );
}

inline float SinPhi(const Vector3f &w) {
	float sinTheta = SinTheta(w);
	return (sinTheta == 0) ? 0 : clamp( -1, 1, w.y / sinTheta);
}

inline float Cos2Phi(const Vector3f &w) { return CosPhi(w) * CosPhi(w); }

inline float Sin2Phi(const Vector3f &w) { return SinPhi(w) * SinPhi(w); }

inline float CosDPhi(const Vector3f &wa, const Vector3f &wb) {
	return clamp(
		(wa.x * wb.x + wa.y * wb.y) / std::sqrt((wa.x * wa.x + wa.y * wa.y) *
		(wb.x * wb.x + wb.y * wb.y)),
		-1, 1);
}

inline Vector3f Reflect(const Vector3f &wo, const Vector3f &n) {
	return -wo + 2 * dotProduct(wo, n) * n;
}

inline bool Refract(const Vector3f &wi, const Vector3f &n, float eta,
	Vector3f *wt) {
	// Compute $\cos \theta_\roman{t}$ using Snell's law
	float cosThetaI = dotProduct(n, wi);
	float sin2ThetaI = std::max(float(0), float(1 - cosThetaI * cosThetaI));
	float sin2ThetaT = eta * eta * sin2ThetaI;

	// Handle total internal reflection for transmission
	if (sin2ThetaT >= 1) return false;
	float cosThetaT = std::sqrt(1 - sin2ThetaT);
	*wt = eta * -wi + (eta * cosThetaI - cosThetaT) * Vector3f(n);
	return true;
}

inline bool SameHemisphere(const Vector3f &w, const Vector3f &wp) {
	return w.z * wp.z > 0;
}


// BSDF Declarations
enum BxDFType {
	BSDF_REFLECTION = 1 << 0,
	BSDF_TRANSMISSION = 1 << 1,
	BSDF_DIFFUSE = 1 << 2,
	BSDF_GLOSSY = 1 << 3,
	BSDF_SPECULAR = 1 << 4,
	BSDF_ALL = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR | BSDF_REFLECTION |
	BSDF_TRANSMISSION,
};

class Fresnel {
public:
	// Fresnel Interface
	virtual ~Fresnel();
	virtual Vector3f Evaluate(float cosI) const = 0;
	virtual std::string ToString() const = 0;
};


class FresnelConductor : public Fresnel {
public:
	// FresnelConductor Public Methods
	Vector3f Evaluate(float cosThetaI) const;
	FresnelConductor(const Vector3f &etaI, const Vector3f &etaT,
		const Vector3f &k)
		: etaI(etaI), etaT(etaT), k(k) {}
	std::string ToString() const;

private:
	Vector3f etaI, etaT, k;
};

class FresnelDielectric : public Fresnel {
public:
	// FresnelDielectric Public Methods
	Vector3f Evaluate(float cosThetaI) const;
	FresnelDielectric(float etaI, float etaT) : etaI(etaI), etaT(etaT) {}
	std::string ToString() const;

private:
	float etaI, etaT;
};

class FresnelNoOp : public Fresnel {
public:
	Vector3f Evaluate(float) const { return Vector3f(1.); }
	std::string ToString() const { return "[ FresnelNoOp ]"; }
};

class Material{
public:

    // Compute reflection direction
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I - 2 * dotProduct(I, N) * N;
    }

    // Compute refraction direction using Snell's law
    //
    // We need to handle with care the two possible situations:
    //
    //    - When the ray is inside the object
    //
    //    - When the ray is outside.
    //
    // If the ray is outside, you need to make cosi positive cosi = -N.I
    //
    // If the ray is inside, you need to invert the refractive indices and negate the normal N
    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior,float& etai_etat) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
        float eta = etai / etat;

		etai_etat = eta;

        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected
    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {  std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }

	
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }

    Vector3f toWorld(const Vector3f &a, const Vector3f &N){
        Vector3f B, C;
        if (std::fabs(N.x) > std::fabs(N.y)){
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vector3f(N.z * invLen, 0.0f, -N.x *invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y *invLen);
        }
        B = crossProduct(C, N);
        return a.x * B + a.y * C + a.z * N;
    }

public:
    MaterialType m_type;
    //Vector3f m_color;
    Vector3f m_emission;
    float ior;
    Vector3f Kd, Ks;
    float specularExponent;
    //Texture tex;

    inline Material(MaterialType t=DIFFUSE, Vector3f e=Vector3f(0,0,0),float _ior=1.0f,Vector3f _kd=Vector3f(0,0,0),Vector3f _ks=Vector3f(0,0,0),float _spec=1.0f);
    inline MaterialType getType();
    //inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    // sample a ray by Material properties
    inline Vector3f sample(const Vector3f &wi, const Vector3f &N,Vector3f& wo, float & pdf);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);

};

Material::Material(MaterialType t, Vector3f e , float _ior , Vector3f _kd , Vector3f _ks , float _spec ){
    m_type = t;
    //m_color = c;
    m_emission = e;

	ior = _ior;
	Kd = _kd;
	Ks = _ks;
	specularExponent = _spec;
}

MaterialType Material::getType(){return m_type;}
///Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}
bool Material::hasEmission() {
    if (m_emission.norm() > EPSILON) return true;
    else return false;
}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}


Vector3f Material::sample(const Vector3f &wo, const Vector3f &N, Vector3f& wi, float & _pdf){
    switch(m_type)
	{
        case DIFFUSE:
        {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);
            wi=toWorld(localRay, N);
            
			wi.normalized();

			_pdf = pdf(wi, wo, N);


			return eval(wi, wo, N);


            break;
        }
		case REFR:
		{

			float etai_etat = 0;

	        Vector3f refractDir = refract(-wo, N, 1.5f,etai_etat);
			Vector3f reflectDir = reflect(-wo, N);

			float kr = 0.0f;
			fresnel(-wo, N, 1.5f, kr);

			//float prob = 0.5 + 0.25*kr;


			if (get_random_float() > kr)
			{

				wi = refractDir.normalized();
				_pdf = 1.0f;
				return 1.0/ abs(dotProduct(wi,N));
			}
			else 
			{
				_pdf = 1.0f;
				wi = reflectDir.normalized();
				return (etai_etat*etai_etat )/ abs(dotProduct(wi, N));
			}

		}
		case SPEC:
		{
			return reflect(-wi, N);
		}
    }
}

float Material::pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample probability 1 / (2 * PI)
            if (dotProduct(wi, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
		case REFR:
		{
			return 1.0f;
		}
		case SPEC:
		{
			return 1.0f;
		}
    }
}

Vector3f Material::eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type)
	{
        case DIFFUSE:
        {
            // calculate the contribution of diffuse   model
            float cosalpha = dotProduct(N, wi);
            if (cosalpha > 0.0f) {
                Vector3f diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return Vector3f(0.0f);
            break;
        }
		case REFR:
		{

		}
		case SPEC:
		{

		}
    }
}

#endif //RAYTRACING_MATERIAL_H
