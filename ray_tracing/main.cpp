#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <windows.h>
#include <WinUser.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <io.h>
#include <fcntl.h>
#include <time.h>
#include "head.h"
#include "vec.h"
#include "device.h"
#include "buffer.h"
#include "Object.h"
#include "Camera.h"
#include "Material.h"
#include "texture.h"
#include "pdf.h"


int width = 800, height = 600;

//int width = 1040, height = 650;
int ns = 1600;
bool quit = false;

vec4 color(const ray& r, const BVHNode *world, HitableObj *light_shape, int depth)
{
	hit_record rec;
	if (world->hit(r, 1e-5, 1e6, rec))
	{
		scatter_record srec;
		vec4 emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
		if (depth < 50 && rec.mat_ptr->scatter(r, rec, srec))
		{
			if (srec.is_specular) {
				return emitted + srec.attenuation*color(srec.specular_ray, world, light_shape, depth + 1);
			}
			else {
				/*HitablePDF plight(light_shape, rec.p);
				MixPDF p(&plight, srec.pdf_ptr);*/
				PDF *p = srec.pdf_ptr;
				ray scattered(rec.p, p->generate(), r.time());
				float pdf = p->value(scattered.direction());
				delete srec.pdf_ptr;
				return emitted + srec.attenuation*rec.mat_ptr->scattering_pdf(r, rec, scattered) * color(scattered, world, light_shape, depth + 1) / pdf;
			}
		}
		else
		{
			return emitted;
		}
	}
	else
	{
		vec4 dir = r.direction().unit();
		float t = 0.5*(dir.y() + 1.0);
		return 0.25*((1. - t)*vec4(1., 1., 1.) + t * vec4(0.5, 0.7, 1.));
	}
}

HitableList *random_scene()
{
	int n = 500;
	HitableObj **list = new HitableObj*[n + 1];

	Texture *checker = new CheckerTexture(new ConstantTexture(vec3(0.2, 0.3, 0.5)), new ConstantTexture(vec3(0.9f)));
	list[0] = new Sphere(vec3(0., -1000, 0.), 1000, new Lambertian(checker));
	list[1] = new Parallelogram(vec4(-4.f, 2.8f, -5.f), vec4(6.f, 6.f, 0.f), vec4(0.0f, 0.0f, 6.0f), new DiffuseLight(new ConstantTexture(vec4(6.f))));
	int i = 2;
	for (int a = -7; a < 7; ++a)
	{
		for (int b = -7; b < 7; ++b)
		{
			float choose_mat = randf();
			vec3 center(a + 0.9*randf(), 0.2, b + randf()*0.9);
			if ((center - vec3(4., 0., 2.)).length() > 0.9)
			{
				float radius = 0.15 + 0.05*randf();
				if (choose_mat < 0.65)
				{
					/*if (choose_mat < 0.2) {
						list[i++] = new MovingSphere(center, center + vec3(0, 0.3*randf(), 0), 0, 1, radius, new Lambertian(new NoiseTexture(4.5+randf())));
					}
					else {
						list[i++] = new MovingSphere(center, center + vec3(0, 0.3*randf(), 0), 0, 1, radius, new Lambertian(new ConstantTexture(vec3(randf(), randf(), randf()))));
					}*/
					list[i++] = new Sphere(center, radius, new Lambertian(new ConstantTexture(vec3(randf(), randf(), randf()))));

				}
				else if (choose_mat < 0.85)
				{
					list[i++] = new Sphere(center, radius, new Metal(0.85*(vec3(randf(), randf(), randf())), 0.4*randf()));
				}
				else
				{
					list[i++] = new Sphere(center, radius, new Dielectric(randf()*0.8 + 1.2));
				}
			}
		}

	}
	float height = 1.0f;
	list[i++] = new Sphere(vec4(0.f, height, 0.f), height, new Dielectric(1.5));
	//list[i++] = new Sphere(vec4(-2.0f, height, 0.0f), height, new Lambertian(new ImageTexture("earth.jpg")));
	HitableObj *obj = new Sphere(vec4(-2.0f, height, 0.0f), height, NULL);
	list[i++] = new ConstantMedium(obj, 0.4f, new ConstantTexture(vec4(1.)));
	list[i++] = new Sphere(vec3(0, height + 2, 0), 0.5*height, new DiffuseLight(new NoiseTexture(5.0)));
	list[i++] = new Sphere(vec4(2.0f, height, 0.0f), height, new Metal(vec4(0.9f, 0.76f, 0.8f), 0.05f));

	return new HitableList(list, i);
}


void writeToPPM(const Buffer &buffer, const char* name) {
	FILE *f = fopen(name, "w");
	fprintf(f, "P3\n%d %d\n255\n", buffer.width, buffer.height);

	for (int j = buffer.height - 1; j >= 0; --j)
	{
		for (int i = 0; i < buffer.width; ++i)
		{
			Color *p = buffer.ptr + buffer.width * j + i;
			int ir = int(p->R * 255);
			int ig = int(p->G * 255);
			int ib = int(p->B * 255);
			//printf("x:%d y:%d \n", i, j);
			fprintf(f, "%d %d %d\n", ir, ig, ib);
		}
	}
	fclose(f);
}

int WINAPI WinMain(HINSTANCE hinstance, HINSTANCE hprevinstance, LPSTR lpcmdline, int ncmdshow)
{
	AllocConsole();
	freopen("CONOUT$", "w+t", stdout);
	freopen("CONIN$", "r+t", stdin);

	MSG msg;
	Device device(hinstance, width, height);
	ShowWindow(device.hwnd, SW_SHOW);
	UpdateWindow(device.hwnd);

	srand(static_cast <unsigned> (time(0)));

	vec4 t1(1, 0.5, 0, 0), t2(0.2, 1, 0, 0);
	float t = dot(t1, t2);



	HitableList *world = random_scene();

	BVHNode *scene = new BVHNode(world->list, world->list_size, 0, 0);

	//vec3 ori(15.,1.8,3.), look_at(2., 0.75, 0.35) , up(0., 1., 0.) ;
	vec3 ori(3.5, 3.8, 25.), look_at(0.f, 0.65f, 0.35f), up(0.f, 1.f, 0.f);
	float focus_distance = (look_at - ori).length();
	//Camera init_camera(ori, look_at,up, 145.f, float(width) / float(height));
	FocusCamera init_camera(ori, look_at, up, 145.f, float(width) / float(height), 0.10f, focus_distance);

	Buffer buffer(width, height);
	for (int k = 0; k < ns; ++k)
	{
		for (int j = height - 1; j >= 0; --j)
		{
			/*time_t start, stop;
			start = time(NULL);*/
			for (int i = 0; i < width; ++i)
			{
				int offset = buffer.width*j + i;
				vec4 col((*(buffer.ptr + offset)).R, (*(buffer.ptr + offset)).G, (*(buffer.ptr + offset)).B);
				float delta_i = randf();
				float delta_j = randf();
				float u = float(i + delta_i) / float(width) * 2 - 1.;
				float v = float(j + delta_j) / float(height) * 2 - 1.;

				ray r = init_camera.get_ray(u, v);
				col += de_nan(color(r, scene, world->list[1], 0) - col) / (k + 1);

				for (int i = 0; i < 3; ++i) {
					if (col[i] > 1) { col[i] = 1; }
				}
				(*(buffer.ptr + offset)).R = col.r();
				(*(buffer.ptr + offset)).G = col.g();
				(*(buffer.ptr + offset)).B = col.b();

				//stop = time(NULL);
				//printf("Use Time:%ld\n", (stop - start));
				//std::cout << ir << " " << ig << " " << ib << "\n";
			}

		}
		device.draw(buffer);
	}
	writeToPPM(buffer, "light_t.ppm");
	device.draw(buffer);
	while (!quit)
	{
		if (PeekMessage(&msg, NULL, NULL, NULL, PM_REMOVE))
		{
			if (msg.message == WM_QUIT)
				quit = true;
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
		//buffer.clear(0.f);
		device.draw(buffer);
	}

	device.release();
	FreeConsole();

	return 0;
}