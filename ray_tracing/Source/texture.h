#pragma once

#include <math.h>
#include <memory>
#include "vec.h"
#include "Noise.h"
#include "stb_image.h"

#ifndef _TEXTURE_
#define _TEXTURE_


class Texture {
public:
	virtual vec4 value(float u, float v, const vec4 &p) const = 0;
	virtual ~Texture() {};
};

class ConstantTexture :public Texture {
public:
	vec4 color;
	ConstantTexture(){}
	ConstantTexture(const vec4 &c) :color(c) {}
	virtual vec4 value(float u, float v, const vec4 &p)const {
		return color;
	}
	
};

class CheckerTexture : public Texture {
public:
	std::shared_ptr<Texture> even, odd;
	CheckerTexture() :even(nullptr), odd(nullptr) {}
	CheckerTexture(std::shared_ptr<Texture> t0, std::shared_ptr<Texture> t1):even(t0), odd(t1) {}
	
	virtual vec4 value(float u, float v, const vec4 &p)const {
		float sines = sin(10 * p.x())*sin(10 * p.y())*sin(10*p.z());
		if (sines < 0.f) { return odd->value(u, v, p); }
		else { return even->value(u, v, p); }
	}
};

class NoiseTexture : public Texture {
public:
	Perlin noise;
	float scale;

	NoiseTexture(){
		scale = 1.f;
	}
	NoiseTexture(float s):scale(s){}
	virtual vec4 value(float u, float v, const vec4 &p) const {
		return (vec4(232,93,25)/255.f)*(2.5f+sin(scale*p.z()+10*noise.turb(p)));
	}
};


class ImageTexture : public Texture {
public:
	unsigned char *image;
	int width, height, nrChannels;

	ImageTexture() :image(nullptr){}
	ImageTexture(const char* file){
		stbi_set_flip_vertically_on_load(true);

		image = NULL;
		image = stbi_load(file, &width, &height, &nrChannels, 0);

		if (image == NULL) {
			std::cout << file << "\nError: Texture load failed!" << std::endl;
			return;
		}		
	}

	virtual vec4 value(float u, float v, const vec4 &P)const {
		int i = u * width;
		int j = v * height;
		if (i < 0) { i = 0; }
		if (j < 0) { j = 0; }
		if (i > width - 1) { i = width - 1; }
		if (j > height - 1) { j = height - 1; }
		int offset = nrChannels * i + nrChannels * width*j;
		float r = int(image[offset]) / 255.0f;
		float g = int(image[offset +1]) / 255.0f;
		float b = int(image[offset +2]) / 255.0f;
		return vec4(r, g, b);
	}

	~ImageTexture() {
		stbi_image_free(image);
	}
};

#endif // !_TEXTURE_

