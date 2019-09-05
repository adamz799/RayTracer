#pragma once

#include "Buffer.h"

Buffer::Buffer():width(0),height(0),ptr(nullptr)
{
}

Buffer::Buffer(int _width, int _height)
{
	ptr = nullptr;
	try {
		ptr = new vec4[_width * _height];
	}
	catch (std::bad_alloc) {
		std::cout << "Bad alloc: memory alloc failed!" << std::endl;
		ptr = nullptr;
	}
	if (ptr) {
		memset(ptr, 0, sizeof(vec4) * _width * _height);
		width = _width;
		height = _height;
	}
}

Buffer::~Buffer()
{
	memset(ptr, 0, sizeof(vec4) * width * height);
	delete[] ptr;
	ptr = nullptr;
}

void Buffer::resize(int _width, int _height)
{
	vec4* temp = ptr;
	try {
		ptr = new vec4[_width * _height];
	}
	catch (std::bad_alloc) {
		std::cout << "Bad alloc: resize failed!" << std::endl;
		ptr = temp;
		return;
	}
	if (!ptr) {
		ptr = temp;
	}
	else {
		memset(ptr, 0, sizeof(vec4) * _width * _height);
		width = _width;
		height = _height;
		delete[] temp;
	}
}

void Buffer::setValue(float v)
{
	int* t = (int*)(&v);
	memset(ptr, *t, sizeof(vec4)*width*height);
}




DepthBuffer::DepthBuffer():ptr(nullptr),width(-1),height(-1){}

DepthBuffer::DepthBuffer(int _width, int _height)
{
	ptr = NULL;
	try {
		ptr = new float[_width * _height];
	}
	catch (std::bad_alloc) {
		std::cout << "Bad alloc: memory alloc failed!" << std::endl;
		ptr = NULL;

	}
	if (ptr) {
		std::fill(ptr, ptr + _width * _height, -2.f);
		//memset(ptr, 0, sizeof(float)*_width*_height);
		width = _width;
		height = _height;
	}
}

DepthBuffer::~DepthBuffer()
{
	memset(ptr, 0, sizeof(float) * width * height);
	delete[] ptr;
	ptr = nullptr;
}

void DepthBuffer::resize(int _width, int _height)
{
	float* temp = ptr;
	try {
		ptr = new float[_width * _height];
	}
	catch (std::bad_alloc) {
		std::cout << "Bad alloc: resize failed!" << std::endl;
		ptr = temp;
		return;
	}
	if (!ptr) {
		ptr = temp;
	}
	else {
		std::fill(ptr, ptr + _width * _height, -2.f);
		width = _width;
		height = _height;
		delete[] temp;
	}
}

void DepthBuffer::setValue(float v)
{
	std::fill(ptr, ptr + width * height, v);
}