#pragma once

#include "vec.h"
#include "ray.h"
#include "utils.h"
#include "head.h"
#include "texture.h"
#include "pdf.h"

#ifndef _MATERIAL_
#define _MATERIAL_


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
	Lambertian(Texture *a);
	virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec)const;
	virtual float scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered)const;
};

class Metal : public Material {
public:
	vec4 albedo;
	float fuzz;
	Metal(const vec4 &a, float f);
	virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec)const;
};


class Dielectric : public Material {
public:
	float ref_idx;
	Dielectric(float ri);
	virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec)const;
};

class DiffuseLight : public Material {
public:
	Texture * emit;
	DiffuseLight(Texture* a);
	virtual vec4 emitted(const ray& r, hit_record& rec, float u, float v, const vec4& p)const;
};

class Isotropic : public Material {
public:
	Texture *albedo;
	Isotropic(Texture* a);
	virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec)const;
};


#endif // !_MATERIAL_

