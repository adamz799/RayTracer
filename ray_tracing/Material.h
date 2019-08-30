#pragma once

#ifndef MATERIAL
#define MATERIAL

#include "head.h"
#include "vec.h"
#include "texture.h"
#include "pdf.h"



class Material {
public:
	virtual bool scatter(const ray &r_in, const hit_record &rec, scatter_record &srec)const {
		return false;
	}
	virtual float scattering_pdf(const ray &r_in, const hit_record &rec, const ray &scattered) const {
		return false;
	}
	virtual vec4 emitted(const ray &r, hit_record &rec, float u, float v, const vec4 &p)const {
		return vec4(0.f);
	}
};

class Lambertian : public Material {
public:
	Texture  *albedo;
	Lambertian(Texture *a) : albedo(a) {};
	virtual bool scatter(const ray &r_in, const hit_record &rec, scatter_record &srec)const {
		srec.is_specular = false;
		srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
		srec.pdf_ptr = new CosinePDF(rec.normal);
		return true;
	}

	virtual float scattering_pdf(const ray &r_in, const hit_record &rec, const ray &scattered)const{
		float cosine = dot(rec.normal, scattered.direction().unit());
		if (cosine < 0) { return 0; }
		return cosine / M_PI;
	}
};

class Metal : public Material {
public:
	vec4 albedo;
	float fuzz;
	Metal(const vec4 &a, float f) : albedo(a) { fuzz = (f < 1.) ? f : 1.; };
	virtual bool scatter(const ray &r_in, const hit_record &rec, scatter_record &srec)const {
		{
			vec4 reflected = reflect(r_in.direction().unit(), rec.normal);
			srec.specular_ray = ray(rec.p, reflected + fuzz * random_in_unit_sphere());
			srec.attenuation = albedo;
			srec.is_specular = true;
			srec.pdf_ptr = NULL;
			return true;
		}
	};
};


class Dielectric : public Material {
public:
	float ref_idx;
	Dielectric(float ri) :ref_idx(ri) {};
	virtual bool scatter(const ray &r_in, const hit_record &rec, scatter_record &srec)const {
		{
			srec.attenuation = vec4(0.95f);
			srec.is_specular = true;
			vec4 outward_normal;
			vec4 reflected;
			float ni_over_nt;
			vec4 refracted;
			float reflect_prob;
			float cosine;

			float dot_val = dot(r_in.direction(), rec.normal);
			if (dot_val > 0) {
				outward_normal = -rec.normal;
				ni_over_nt = ref_idx;
				cosine = ref_idx * dot_val / r_in.direction().length();
			}
			else
			{
				outward_normal = rec.normal;
				ni_over_nt = 1. / ref_idx;
				cosine = -dot_val / r_in.direction().length();
			}
			reflected = reflect(r_in.direction(), outward_normal);

			if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted))
			{
				reflect_prob = schlick(cosine, ref_idx);
				/*if (dot(r_in.direction(), rec.normal) > 0) {
					reflect_prob = 0;
				}*/
			}
			else
			{
				srec.specular_ray = ray(rec.p, reflected);
				reflect_prob = 1.;
				//return true;
			}
			if (randf() < reflect_prob)
			{
				srec.specular_ray = ray(rec.p, reflected);
			}
			else
			{
				srec.specular_ray = ray(rec.p, refracted);
			}
			return true;
		}
	};
};

class DiffuseLight : public Material {
public:
	Texture * emit;

	DiffuseLight(Texture *a): emit(a){}
	virtual vec4 emitted(const ray &r, hit_record &rec, float u, float v, const vec4 &p)const {
		return emit->value(u, v, p);
	}


};

class Isotropic : public Material {
public:
	Texture *albedo;

	Isotropic(Texture *a):albedo(a){}

	virtual bool scatter(const ray &r_in, const hit_record &rec, scatter_record &srec)const {
		srec.is_specular = true;
		srec.specular_ray = ray(rec.p, random_in_unit_sphere());
		srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
		return true;
	}

};


#endif // !MATERIAL

