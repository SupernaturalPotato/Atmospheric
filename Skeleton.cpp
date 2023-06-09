#include "framework.h"

const char * const vertexSource = R"(
#version 330 core

out vec2 delta;

void main(){

	const vec2 quadVertices[6] = vec2[6] (
		vec2(-1.f, 1.f), vec2(-1.f, -1.f), vec2(1.f, -1.f),
		vec2(-1.f, 1.f), vec2(1.f, -1.f), vec2(1.f, 1.f)
	);
	
	delta = quadVertices[gl_VertexID];
	gl_Position = vec4(quadVertices[gl_VertexID], 0.f, 1.f);	
}
)";

const char* const fragmentSource = R"(
#version 330 core
precision highp float;

in vec2 delta;
out vec3 color;

uniform vec3 camDir;
uniform vec3 towardsTheSun;

const float EARTH_RADIUS = 6371.f;
const vec3 EARTH_CENTER = vec3(0.f, 0.f, 0.f );
const vec3 CAM_POS = vec3(0.f, 6371.010f, 0.f);

const vec3 betaRayleigh = vec3(5.80e-6f, 1.65e-5f, 3.31e-5f);
const vec3 betaExtinctionRayleigh = betaRayleigh;
const vec3 betaMie = vec3(21e-6f, 21e-6f, 21e-6f);
const vec3 betaExtinctionMie = (1.f / 0.9f) * betaMie;
const vec3 betaExtinctionOzone = vec3(1e-12f, 1e-10f, 1e-9f);
vec3 SUN_STRENGTH = 24.0f * vec3(1.f, 1.f, 1.f);

const float PI = 3.14159265f;
const int INTEGRATION_STEPS = 20;
const float ATMOSPHERE_DENSITY_SCALING = 1.2f;

float phaseRayleigh(float theta) { 
	const float normalizationFactor = 3.f / (16.f * PI);
	float ct = cos(theta); 
	return normalizationFactor * (1.f + ct * ct);
}

float phaseHenyeyGreenstein(float theta, float g) {
	
	const float normalizationFactor = 3.f / (8.f * PI);
	float cosTheta = cos(theta);
	float secondDenom = 1.f + g * g - 2.f * g * cosTheta;

	return normalizationFactor *
	3.f * (1.f - g * g) * (1.f + cosTheta * cosTheta) 
	/ ( 2.f * ( 2.f + g * g) * sqrt(secondDenom * secondDenom * secondDenom) );
}

float densityRayleigh(float height) { return ATMOSPHERE_DENSITY_SCALING * exp(height / -8.f); }
float densityMie(float height) { return ATMOSPHERE_DENSITY_SCALING * exp(height / -1.2f); }
float densityOzone(float height) { return 6e-7f * densityRayleigh(height); }
float height(vec3 point) { return length(point) - EARTH_RADIUS; }

bool solveQuadratic(float a, float b, float c, out float outX1, out float outX2) {

	float D = b * b - 4.f * a * c;
	if (D < 0.0f) return false;

	D = sqrt(D);
	float denom = 2.f * a;

	outX1 = (-b + D) / (denom + 0.00001f);
	outX2 = (-b - D) / (denom + 0.00001f);
	return true;
}

float forwardRaySphere(vec3 rayStart, vec3 rayDir, float radius) {

	float t1, t2;

	float a = dot(rayDir, rayDir);
	float b = 2.f * dot(rayStart, rayDir);
	float c = dot(rayStart, rayStart) - radius * radius;

	if (solveQuadratic(a, b, c, t1, t2) == false) return -1.f;
	if (t1 > 0.f && t2 > 0.f) return t1 < t2 ? t1 : t2;
	if (t1 > 0.f) return t1;
	if (t2 > 0.f) return t2;
	return -1.f;
}

vec3 transmittanceBetween(vec3 pa, vec3 pb) {

	float L = distance(pa, pb);
	if (L < 0.02f) return vec3(0.f, 0.f, 0.f);

	float ds = L / INTEGRATION_STEPS * 1000.f;
	vec3 p = pa;
	float h = height(p);
	vec3 step = (pb - pa) * (1.f / INTEGRATION_STEPS);
	float totalDensityRay = 0.f, totalDensityMie = 0.0f, totalDensityOzone = 0.0f;

	p += 0.5f * step;

	for (int i = 0; i < INTEGRATION_STEPS; ++i) {

		h = height(p);
		if(h < 0.f) return vec3(10000.f, 10000.f, 10000.f);
		totalDensityRay += densityRayleigh(h);
		totalDensityMie += densityMie(h);
		totalDensityOzone += densityOzone(h);

		p += step;
	}

	totalDensityRay *= ds;
	totalDensityMie *= ds;
	totalDensityOzone *= ds;

	return totalDensityRay * betaRayleigh + totalDensityMie * betaExtinctionMie + totalDensityOzone * betaExtinctionOzone;
}

vec3 intensitySingleScattering(vec3 p0, vec3 v, vec3 l, vec3 incidentLight) {
	
	float hitScale = forwardRaySphere(p0, -v, EARTH_RADIUS);
	if (hitScale > -0.5f) {
		return vec3(0.07f, 0.07f, 0.07f);
	}

	float theta = acos(dot(v, -l) * 0.99f);

	//the edge of the atmosphere
	hitScale = forwardRaySphere(p0, -v, EARTH_RADIUS + 100.f);
	vec3 pb = p0 + hitScale * (-v);

	//initialize the start point to the ray start point
	vec3 pa = p0;
	
	float L = length(pb - pa);
	float ds = L / INTEGRATION_STEPS * 1000.f;
	vec3 step = (pb - pa) * (1.f / INTEGRATION_STEPS);

	vec3 totalRayleigh = vec3(0.f, 0.f, 0.f);
	vec3 totalMie = vec3(0.f, 0.f, 0.f);

	vec3 p = pa;
	p += step * 0.5f;

	for (int i = 0; i < INTEGRATION_STEPS; ++i) {
		
		float h = height(p);

		float scalePc = forwardRaySphere(p, -l, EARTH_RADIUS + 100.f);	
		vec3 pc = p + -l * scalePc;

		vec3 curTransPcP = transmittanceBetween(pc, p);
		vec3 curTransPPa = transmittanceBetween(p, p0);

		vec3 expTransmittance = exp(-curTransPcP - curTransPPa);
		totalRayleigh += densityRayleigh(h) * expTransmittance;
		totalMie += densityMie(h) * expTransmittance;

		p += step;
	}

	totalRayleigh *= ds * SUN_STRENGTH;
	totalMie *= ds * SUN_STRENGTH;

	vec3 ret = SUN_STRENGTH * (
		totalRayleigh * phaseRayleigh(theta) * betaRayleigh +
		totalMie * phaseHenyeyGreenstein(theta, -0.8f) * betaMie);

	if (dot(v, l) > 0.9993f) {
		vec3 directTransmittance = transmittanceBetween(pb, p0);
		ret += SUN_STRENGTH * SUN_STRENGTH * PI * exp(-directTransmittance);
	}
	return ret;
}

void main() {

	float tanHalfFov = tan(0.5f * 1.5f);
	vec3 camRight = normalize(cross( camDir, vec3(0.f, 1.f, 0.f)));
	vec3 camUp = cross(camRight, camDir);

	camRight *= tanHalfFov * 16.f / 9.f;
	camUp *= tanHalfFov;
	
	vec3 rayDir = normalize(normalize(camDir) + camRight * delta.x + camUp * delta.y);

	color = intensitySingleScattering(CAM_POS, -rayDir, -normalize(towardsTheSun), SUN_STRENGTH);
	color = 1.f - exp(-0.02f * color);
}

)";

struct Camera {
	vec3 forward, right, up;

	float pitch = 0.f;
	float yaw = 0.f;
	float fovY = 1.f;

	void update() {

		forward = vec3(cosf(pitch) * -sinf(yaw), sinf(pitch), -cosf(pitch) * cosf(yaw));
		right = normalize(cross(forward, vec3(0.f, 1.f, 0.f)));
		up = cross(right, forward);
	}
};

GPUProgram gpuProgram;
GLuint uniCamDir, uniTowardsTheSun;
Camera camera;

unsigned HVBO;

bool keysDown[256];

void tryBuffers() {

	std::vector<vec3> vertices = {
		vec3(-1.f, -1.f, 0.f),
		vec3(1.f, -1.f, 0.f),
		vec3(0.f, 1.f, 0.f)
	};
	std::vector<GLuint> indices = {
		0, 1, 2
	};

	glGenBuffers(1, &HVBO);

	glBindBuffer(GL_ARRAY_BUFFER, HVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vec3) * 3 + sizeof(GLuint) * 3, nullptr, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, 3 * sizeof(vec3), vertices.data());
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, HVBO);
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 3 * sizeof(vec3), 3 * sizeof(GLuint), indices.data());	

	glBindBuffer(GL_ARRAY_BUFFER, HVBO);
	
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), 0);

}

void drawTriangle() {

	glBindBuffer(GL_ARRAY_BUFFER, HVBO);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, HVBO);
	glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, (void*)9);
}

void onInitialization() {

	for (bool& k : keysDown) k = false;

	camera.pitch = 0.4f;
	camera.fovY = 1.5f;
	camera.yaw = 0.f;
	camera.update();

	gpuProgram.create(vertexSource, fragmentSource, "outColor");
	gpuProgram.Use();
	uniCamDir = glGetUniformLocation(gpuProgram.getId(), "camDir");
	uniTowardsTheSun = glGetUniformLocation(gpuProgram.getId(), "towardsTheSun");

	glViewport(0, 0, windowWidth, windowHeight);
	glEnable(GL_FRAMEBUFFER_SRGB);
	glClearColor(0, 0, 0, 1);

	tryBuffers();
}

void onDisplay() {

	glClear(GL_COLOR_BUFFER_BIT);
	glDrawArrays(GL_TRIANGLES, 0, 6);
	drawTriangle();
	glutSwapBuffers();
}

void onKeyboard(unsigned char key, int pX, int pY) { 
	keysDown[key] = true; 
}

void onKeyboardUp(unsigned char key, int pX, int pY) {
	keysDown[key] = false;
}

float lastFrame = 0.f;
bool leftPressed = false;
int lastX = 0, lastY = 0;
float angle = 0.0f;
float tTime = 0.f;

void onMouseMotion(int pX, int pY) {

	float dX = (float)(pX - lastX);
	float dY = (float)(pY - lastY);

	if (leftPressed) {
		camera.pitch += dY * 0.0025f;
		camera.yaw += dX * 0.0025f;
		camera.update();
	}

	lastX = pX;
	lastY = pY;
}

void onMouse(int button, int state, int pX, int pY) {

	leftPressed = state == GLUT_DOWN;
	lastX = pX; lastY = pY;
}

void onIdle() {

	tTime = (float)glutGet(GLUT_ELAPSED_TIME);
	float dt = (tTime - lastFrame) / 1000.f;
	
	if (keysDown['r']) angle += 0.0625f * dt;
	if (keysDown['f']) angle -= 0.0625f * dt;

	glUniform3fv(uniCamDir, 1, &camera.forward.x);
	glUniform3f(uniTowardsTheSun, 0.f, sinf(angle), -cosf(angle));

	glutPostRedisplay();
	
	lastFrame = tTime;
}
