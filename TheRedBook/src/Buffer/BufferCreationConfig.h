#pragma once

struct VertObj {
	int PointerID;
	int dataType;
	int size;	
	int stride;
};

struct BufferCreationConfig {

	std::vector<VertObj> vertObj;	
};