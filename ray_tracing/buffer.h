#pragma once

//#include <vector>
#include "vec.h"

#ifndef _BUFFER_
#define _BUFFER_

class Buffer
{
public:
	vec4* ptr;//RGBA in float
	int width, height;
	vec4 operator[](int i) const { return ptr[i]; }
	vec4& operator[](int i) { return ptr[i]; }

	Buffer() :width(0), height(0), ptr(nullptr) {}
	Buffer(int _width, int _height);
	~Buffer();

	void resize(int _width, int _height);
	void setValue(float v);
};

class DepthBuffer
{
public:
	float* ptr;//Depth in float
	int width, height;
	float operator[](int i) const { return ptr[i]; }
	float& operator[](int i) { return ptr[i]; }

	DepthBuffer() :ptr(nullptr), width(-1), height(-1) {}
	DepthBuffer(int _width, int _height);
	~DepthBuffer();

	void resize(int _width, int _height);
	void setValue(float v);
};

#endif // !_BUFFER_

