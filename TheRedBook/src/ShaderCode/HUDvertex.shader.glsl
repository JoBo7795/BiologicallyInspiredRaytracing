#version 460 core

layout(location = 0) in vec2 aPos;
layout(location = 2) in vec2 aTexCoords;

out vec2 fTexCoords;
out vec3 camPos;

uniform vec3 viewPos;

void main()
{
	gl_Position = vec4(aPos,0.0f, 1.0f);	

	fTexCoords = aTexCoords;
	camPos = viewPos;
}