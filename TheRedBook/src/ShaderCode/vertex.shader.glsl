#version 460 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec2 aTexCoords;
layout (location = 3) in vec3 aNormal;


out vec3 fPos;
out vec3 camPos;
out vec4 fVecColor;
out vec2 fTexCoords;
out vec3 fNormal;
out mat4 fMvp;

uniform mat4 translation = mat4(1.0f),rotation = mat4(1.0f),scale = mat4(1.0f);
uniform mat4 view;
uniform vec3 viewPos;
uniform mat4 projection;

void main()
{

	gl_Position = vec4(aPos, 1.0f);//projection * view * (translation * rotation * scale) * vec4(aPos, 1.0f);

	fVecColor = vec4(aColor,1.0f);

	fPos = aPos;
	fTexCoords = aTexCoords;
	fNormal = aNormal;
	fMvp =  view;// *vec4(aPos, 1.0f);
	camPos = viewPos;
}