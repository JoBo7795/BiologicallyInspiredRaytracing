#version 460 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aColor;
layout(location = 2) in vec2 aTexCoords;
layout(location = 3) in vec3 aNormal;
layout(location = 4) in mat4 aInstanceMatrix;
layout(location = 8) in vec2 aTexModOffset;


out vec3 fPos;
out vec4 fVecColor;
out vec2 fTexCoords;
out vec3 fNormal;

uniform mat4 translation = mat4(1.0f), rotation = mat4(1.0f), scale = mat4(1.0f);
uniform mat4 view;
uniform mat4 projection;

void main()
{
	float xOffset, yOffset;

	gl_Position = projection * view * aInstanceMatrix * vec4(aPos, 1.0f);

	fVecColor = vec4(aColor, 1.0f);

	fPos = aPos;
	fTexCoords = aTexCoords / 5 + aTexModOffset;
	fNormal = aNormal;
}