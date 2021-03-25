//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}


Vector3f Scene::Li(const Ray &ray, int depth) const
{
	Vector3f L(0.0);
	Vector3f pathThrought(1.0f);
	Ray pathRay(ray.origin,ray.direction);

	for (auto bounce=0;;bounce++)
	{
		Intersection intersection = intersect(pathRay);
		if (pathThrought.isZero())
			break;
		if(!intersection.happened || bounce>maxDepth)
			break;
		// Explicitly sample light sources
		L += pathThrought*estimateDirectLight(-pathRay.direction,intersection);

		Vector3f wo = normalize(-pathRay.direction);
		Vector3f p = intersection.coords;
		Vector3f n = normalize(intersection.normal);
		Vector3f wi;
		float pdf = 0.0f;
		Vector3f f = intersection.m->sample(wo, n, wi, pdf);

		pathThrought = pathThrought* f*abs(dotProduct(wi, n)) / pdf;
		pathRay = Ray(p + Vector3f(0.001,0.001,0.001)*wi, wi);

		if (bounce > 5)
		{
			float RR =std::min(1.0f, Vector3f::Luminance(pathThrought));
			if (get_random_float() > RR)
				break;
			pathThrought =pathThrought/RR;
		}
	}	
	return L;
}

Vector3f Scene::estimateDirectLight(const Vector3f& wo, const Intersection& intersection)const
{
	Vector3f L(0.0);
	Vector3f p = intersection.coords;

	float lightPdf=0.0f;
	Intersection inter;
	sampleLight(inter, lightPdf);
	Vector3f x = inter.coords;
	Vector3f wi = normalize(x - p);
	Vector3f ni = normalize(inter.normal);
	Vector3f Li = inter.emit;

	lightPdf *= (dotProduct(x-p,x-p)/dotProduct(ni,-wi));
	Vector3f f = intersection.m->eval(wi, wo, intersection.normal)*
		                             abs(dotProduct(wi,intersection.normal));

	Intersection visibility = intersect(Ray(p, wi));
	if (!f.isZero() && (visibility.coords-x).norm()<0.0001f )
	{
		L += f*Li / lightPdf;
	}

	return L;
}
