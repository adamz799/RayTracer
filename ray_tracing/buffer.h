#pragma once

#include <vector>
#include "head.h"
#include "utils.h"

#ifndef _BUFFER_
#define _BUFFER_

class Buffer
{
public:
	vec4* ptr;//RGBA in float
	int width, height;
	vec4 operator[](int i) const { return ptr[i]; }
	vec4& operator[](int i) { return ptr[i]; }

	Buffer();
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

	DepthBuffer();
	DepthBuffer(int _width, int _height);
	~DepthBuffer();

	void resize(int _width, int _height);
	void setValue(float v);
};

#endif // !_BUFFER_

