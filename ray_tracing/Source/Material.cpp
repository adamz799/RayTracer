#pragma once

#include "Material.h"
#include "vec.h"
#include "ray.h"
#include "utils.h"
#include "head.h"
#include "texture.h"
#include "pdf.h"

Lambertian::Lambertian(Texture* a) : albedo(a) {};

bool Lambertian::scatter(const ray& r_in, const hit_record& rec, scatter_record& srec)const {
	srec.is_specular = false;
	srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
	srec.pdf_ptr = new CosinePDF(rec.normal);
	return true;
}

float Lambertian::scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered)const {
	float cosine = dot(rec.normal, scattered.direction().unit());
	if (cosine < 0) { return 0; }
	return cosine / M_PI;
}



Metal::Metal(const vec4& a, float f) : albedo(a) { fuzz = (f < 1.f) ? f : 1.f; }

bool Metal::scatter(const ray& r_in, const hit_record& rec, scatter_record& srec)const {
	vec4 reflected = reflect(r_in.direction().unit(), rec.normal);
	srec.specular_ray = ray(rec.p, reflected + fuzz * random_in_unit_sphere());
	srec.attenuation = albedo;
	srec.is_specular = true;
	srec.pdf_ptr = NULL;
	return true;

}


Dielectric::Dielectric(float ri) :ref_idx(ri) {}

bool Dielectric::scatter(const ray & r_in, const hit_record & rec, scatter_record & srec)const {
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


DiffuseLight::DiffuseLight(Texture* a) : emit(a) {}

vec4 DiffuseLight::emitted(const ray& r, hit_record& rec, float u, float v, const vec4& p)const {
	return emit->value(u, v, p);
}




Isotropic::Isotropic(Texture* a) :albedo(a) {};

bool Isotropic::scatter(const ray& r_in, const hit_record& rec, scatter_record& srec)const {
	srec.is_specular = true;
	srec.specular_ray = ray(rec.p, random_in_unit_sphere());
	srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
	return true;
}
