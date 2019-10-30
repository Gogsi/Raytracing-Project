#include "flyscene.hpp"
#include <GLFW/glfw3.h>
#include <queue>
#include <thread>
#include <ctime>
#include "ShaderLib.h"


using namespace std;

#define MAX_BOUNCES 2

void Flyscene::initialize(int width, int height) {
	// initiliaze the Phong Shading effect for the Opengl Previewer
	phong.initialize();

	// set the camera's projection matrix
	flycamera.setPerspectiveMatrix(60.0, width / (float)height, 0.1f, 100.0f);
	flycamera.setViewport(Eigen::Vector2f((float)width, (float)height));

  // load the OBJ file and materials
  Tucano::MeshImporter::loadObjFile(mesh, materials,
                                    "resources/models/planeobj.obj");

	// normalize the model (scale to unit cube and center at origin)
	mesh.normalizeModelMatrix();

	// pass all the materials to the Phong Shader
	for (int i = 0; i < materials.size(); ++i)
		phong.addMaterial(materials[i]);


	for (int i = 0; i < mesh.getNumberOfFaces(); ++i) {
		triangles.push_back(mesh.getFace(i));
	}


	// set the color and size of the sphere to represent the light sources
	// same sphere is used for all sources
	lightrep.setColor(Eigen::Vector4f(1.0, 1.0, 0.0, 1.0));
	lightrep.setSize(0.15);

	// create a first ray-tracing light source at some random position
	lights.push_back(Eigen::Vector3f(-1.0, 1.0, 1.0));
	//sphericalLights.push_back(make_pair(Eigen::Vector3f(-1.0, 1.0, 1.0), 0.1));

	// scale the camera representation (frustum) for the ray debug
	camerarep.shapeMatrix()->scale(0.2);

	// the debug ray is a cylinder, set the radius and length of the cylinder
	ray.setSize(0.005, 10.0);

	// craete a first debug ray pointing at the center of the screen
	createDebugRay(Eigen::Vector2f(width / 2.0, height / 2.0));

	glEnable(GL_DEPTH_TEST);


	// for (int i = 0; i<mesh.getNumberOfFaces(); ++i){
	//   Tucano::Face face = mesh.getFace(i);    
	//   for (int j =0; j<face.vertex_ids.size(); ++j){
	//     std::cout<<"vid "<<j<<" "<<face.vertex_ids[j]<<std::endl;
	//     std::cout<<"vertex "<<mesh.getVertex(face.vertex_ids[j]).transpose()<<std::endl;
	//     std::cout<<"normal "<<mesh.getNormal(face.vertex_ids[j]).transpose()<<std::endl;
	//   }
	//   std::cout<<"mat id "<<face.material_id<<std::endl<<std::endl;
	//   std::cout<<"face   normal "<<face.normal.transpose() << std::endl << std::endl;
	// }


  std::cout << "..." << std::endl;
  // create the array of boxes
  this->root_box = Box(mesh);
  
  // KD Trees :
 // divideBox_KD(1000);
 //#define show_KD

  // Flat structure:
  this->boxes = divideBox(root_box, 100);
  #define show_flat
  
  // if u want to visualize the bounding boxes 
  //#define show_bounding

  initBoundingBoxes();

}

void Flyscene::paintGL(void) {


	// update the camera view matrix with the last mouse interactions
	flycamera.updateViewMatrix();
	Eigen::Vector4f viewport = flycamera.getViewport();

	// clear the screen and set background color
	glClearColor(0.9, 0.9, 0.9, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	scene_light.resetViewMatrix();
	scene_light.viewMatrix()->translate(-lights.back());
	//scene_light.viewMatrix()->translate(-sphericalLights.back().first);

  // render the scene using OpenGL and one light source
  phong.render(mesh, flycamera, scene_light);

  // render the bounding boxes:
#ifdef show_bounding
  renderBoundingBoxes();	
#endif
  // render the ray and camera representation for ray debug

  for (auto i = 0; i < rays.size(); i++)
  {
	  rays.at(i).render(flycamera, scene_light);
  }

  for (auto i = 0; i < ray_hitbox.size(); i++) {
	  ray_hitbox.at(i).render(flycamera, scene_light);
  }


  //ray.render(flycamera, scene_light);
	camerarep.render(flycamera, scene_light);

	// render the ray and camera representation for ray debug
	ray.render(flycamera, scene_light);
	camerarep.render(flycamera, scene_light);

	// render ray tracing light sources as yellow spheres
	for (int i = 0; i < lights.size(); ++i) {
		lightrep.resetModelMatrix();
		lightrep.modelMatrix()->translate(lights[i]);
		lightrep.render(flycamera, scene_light);
	}

	for (int i = 0; i < sphericalLights.size(); ++i) {
		lightrep.resetModelMatrix();
		lightrep.modelMatrix()->translate(sphericalLights[i].first);
		lightrep.render(flycamera, scene_light);
	}

	// render coordinate system at lower right corner
	flycamera.renderAtCorner();
}

void Flyscene::simulate(GLFWwindow* window) {
	// Update the camera.
	// NOTE(mickvangelderen): GLFW 3.2 has a problem on ubuntu where some key
	// events are repeated: https://github.com/glfw/glfw/issues/747. Sucks.
	float dx = (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS ? 0.1 : 0.0) -
		(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS ? 0.1 : 0.0);
	float dy = (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS ||
		glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS
		? 0.1
		: 0.0) -
		(glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS ||
			glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS
			? 0.1
			: 0.0);
	float dz = (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS ? 0.1 : 0.0) -
		(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS ? 0.1 : 0.0);
	flycamera.translate(dx, dy, dz);
}

void Flyscene::createDebugRay(const Eigen::Vector2f& mouse_pos) {
	rays.clear();
	ray_hitbox.clear();

	//ray.resetModelMatrix();
	// from pixel position to world coordinates
	Eigen::Vector3f screen_pos = flycamera.screenToWorld(mouse_pos);

	// direction from camera center to click position
	Eigen::Vector3f dir = (screen_pos - flycamera.getCenter()).normalized();

	ReflectDebugRay(flycamera.getCenter(), dir, 0);
  
	// position and orient the cylinder representing the ray
	//ray.setOriginOrientation(flycamera.getCenter(), dir);

	// place the camera representation (frustum) on current camera location,
	camerarep.resetModelMatrix();
	camerarep.setModelMatrix(flycamera.getViewMatrix().inverse());
}

// Debug Ray
void Flyscene::ReflectDebugRay(Eigen::Vector3f origin, Eigen::Vector3f dir, int bounce) {

	Tucano::Shapes::Cylinder reflectionRay = Tucano::Shapes::Cylinder(0.1, 1.0, 16, 64);

	reflectionRay.setOriginOrientation(origin, dir);

	int max_bounce = 10;

	if (bounce > max_bounce) {
		return;
	}

	Ray r = Ray(origin, dir);
	Tucano::Face closest_triangle;
	HitInfo smallestHit;
	float smallestT = INFINITY;

	for (auto i = 0; i < boxes.size(); i++)
	{
		Box curr_box = boxes.at(i);
		if (curr_box.isInBox(origin)) {
			HitInfo result_triangle = intersectTriangle(curr_box.triangles, origin, dir);
			if (result_triangle.t != INFINITY && smallestT > result_triangle.t) {
				smallestT = result_triangle.t;
				smallestHit = result_triangle;
				closest_triangle = curr_box.triangles.at(result_triangle.faceId);
			}
			else {
				continue;
			}
		}
		else {
			HitInfo result_box = intersectBox(curr_box, origin, dir);


			if (result_box.t != INFINITY && result_box.t > 0)
			{

				ray_hitbox.push_back(bounding_boxes.at(i));

				HitInfo result_triangle = intersectTriangle(curr_box.triangles, origin, dir);
				if (result_triangle.t != INFINITY && smallestT > result_triangle.t) {
					smallestT = result_triangle.t;
					smallestHit = result_triangle;
					closest_triangle = curr_box.triangles.at(result_triangle.faceId);
				}
			}
		}
	}

	if (smallestT == INFINITY) {
		reflectionRay.setSize(0.0025, 10.0);
		reflectionRay.setColor(Eigen::Vector4f(1.0, 0.5, 0.5, 0.0));
		rays.push_back(reflectionRay);

		return;
	}
	else {
		reflectionRay.setSize(0.0025, smallestT);
		reflectionRay.setColor(Eigen::Vector4f(1.0, 1.0, 1.0, 0.0));
		rays.push_back(reflectionRay);


		for (auto i = 0; i < lights.size(); i++)
		{
			if (canSeeLight(lights.at(i), smallestHit.point)) {
				Tucano::Shapes::Cylinder rayToLight = Tucano::Shapes::Cylinder(0.1, 1.0, 16, 64);
				Vector3f direction = (lights.at(i) - smallestHit.point).normalized();

				rayToLight.setOriginOrientation(smallestHit.point, direction);
				rayToLight.setColor(Eigen::Vector4f(0.0, 0.0, 0.0, 0.0));

				rayToLight.setSize(0.0025, (lights.at(i) - smallestHit.point).norm());
				rays.push_back(rayToLight);
			}
		}

		Ray newRay = r.reflectRay(smallestHit.normal, smallestHit.point);
		Eigen::Vector3f origin2 = newRay.getOrigin();
		Eigen::Vector3f dir2 = newRay.getDirection();
		return ReflectDebugRay(origin2, dir2, bounce + 1);
	}
}


void Flyscene::raytraceScene(int width, int height) {
  std::cout << "ray tracing ..." << std::endl;
  std::clock_t c_start = std::clock();
  auto timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
  std::cout << "Time at the start : " << ctime(&timenow) << std::endl;

  // if no width or height passed, use dimensions of current viewport
  Eigen::Vector2i image_size(width, height);
  if (width == 0 || height == 0) {
    image_size = flycamera.getViewportSize();
  }

  // create 2d vector to hold pixel colors and resize to match image size
  vector<vector<Eigen::Vector3f>> pixel_data;
  pixel_data.resize(image_size[1]);
  for (int i = 0; i < image_size[1]; ++i)
    pixel_data[i].resize(image_size[0]);

  // origin of the ray is always the camera center
  Eigen::Vector3f origin = flycamera.getCenter();
  Eigen::Vector3f screen_coords;

  // Threads
  int number_threads = 16;
  vector<thread> threads;

  for (auto i = 0; i < number_threads; i++) {
	  thread curr_thread(&Flyscene::updating_pixels, this, std::ref(pixel_data), std::ref(origin), std::ref(image_size), number_threads, i);
	  threads.push_back(std::move(curr_thread));
  }

  for (auto i = 0; i < threads.size(); i++) {
	  //std::cout << i << std::endl;
	  threads.at(i).join();
	  //std::cout << i << std::endl;
  }

  // write the ray tracing result to a PPM image
  Tucano::ImageImporter::writePPMImage("result.ppm", pixel_data);
  std::clock_t c_end = std::clock();
  auto timeend = chrono::system_clock::to_time_t(chrono::system_clock::now());
  std::cout << "Time at the end : " << ctime(&timeend) << std::endl;
  std::cout << "Time elapsed: " << 1000*(c_end - c_start) / (CLOCKS_PER_SEC) << " ms"<< std::endl;
  std::cout << "ray tracing done! " << std::endl;}

void Flyscene::updating_pixels(vector<vector<Eigen::Vector3f>>& pixel_data, Eigen::Vector3f& origin, Eigen::Vector2i& image_size, int number_threads, int thread_id) {

	sphereFace = Tucano::Face();
#define SCALE 0.5f

	std::vector<Eigen::Vector4f> spherePositions{
		Eigen::Vector4f(-0.5f, 0.25f, 0.25f, 0.75f) * SCALE,
		Eigen::Vector4f(0.75f, 0.25f, -0.1f, 0.75f) * SCALE,
		Eigen::Vector4f(0.0f, 0.33f, -0.5f, 1.0f) * SCALE,
		Eigen::Vector4f(-0.75f, 0.33f, -0.2f, 1.0f) * SCALE,
		Eigen::Vector4f(0.6f, 0.167f, 0.4f, 0.5f) * SCALE
	};
	std::vector<Eigen::Vector3f> sphereColors{
		Eigen::Vector3f(0.8f, 0.2401f, 0.0f), // red sphere
		Eigen::Vector3f(0.5f, 0.5f, 0), // green sphere
		Eigen::Vector3f(0.98824, 0.72941, 0.01176), // orange sphere
		Eigen::Vector3f(0.0f, 0.5f, 0.5f), // blue sphere
		Eigen::Vector3f(1.0f, 0.0f, 1.0f) // pink sphere
	};

	// HAS TO BE SAME SIZE
	if (spherePositions.size() != sphereColors.size()) {
		abort();
	}

	for (size_t i = 0; i < spherePositions.size(); i++)
	{
		auto tempMat = currentScene.addSphere(spherePositions[i].head<3>(), spherePositions[i].w(), sphereColors[i]);
		int mat_id = materials.size();
		materials.push_back(tempMat);
		currentScene.getSphere(i).setMaterialID(mat_id);
	}

	for (int j = thread_id; j < image_size[1]; j += number_threads) {
		for (int i = 0; i < image_size[0]; i++) {
			// create a ray from the camera passing through the pixel (i,j)
			Eigen::Vector3f screen_coords = flycamera.screenToWorld(Eigen::Vector2f(i, j));
			// launch raytracing for the given ray and write result to pixel data
			pixel_data[i][j] = traceRay(0, Ray(origin, screen_coords - origin));
		}
	}
}

void Flyscene::renderBoundingBoxes() {
#ifdef show_flat 
	for (auto i = 0; i < bounding_boxes.size(); i++)
	{
		bounding_boxes.at(i).render(flycamera, scene_light);
	}
#endif // show_flat

#ifdef show_KD 
	for (auto i = 0; i < bounding_boxes.size(); i++)
	{
		bounding_boxes.at(i).render(flycamera, scene_light);
	}
#endif
}


void Flyscene::initBoundingBoxes() {

	for (auto i = 0; i < boxes.size(); i++)
	{
		Box curr_box = this->boxes.at(i);
		Tucano::Shapes::Box helper = Tucano::Shapes::Box(curr_box.tmax.x() - curr_box.tmin.x(), curr_box.tmax.y() - curr_box.tmin.y(),
			curr_box.tmax.z() - curr_box.tmin.z());

		Affine3f matrix = Affine3f::Identity();

		mesh.resetModelMatrix();
		mesh.setModelMatrix(matrix);

		Eigen::Vector3f translation_vector = Eigen::Vector3f((curr_box.tmax.x() + curr_box.tmin.x()) / 2, (curr_box.tmax.y() + curr_box.tmin.y()) / 2,
			(curr_box.tmax.z() + curr_box.tmin.z()) / 2);


		matrix.translate(translation_vector);
		helper.setModelMatrix(matrix);


		float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float g = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float b = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

		Eigen::Vector4f color = Vector4f(r, g, b, 0.1);
		helper.setColor(color);
		this->bounding_boxes.push_back(helper);
	}	
}

Eigen::Vector3f Flyscene::traceRay(int bounce, Ray ray) {

#ifndef KD
	//std::cout << boxes.size() << std::endl;
	Tucano::Face closest_triangle;
	HitInfo smallestHit;
	float smallestT = INFINITY;

	for (size_t i = 0; i < currentScene.getNumberOfSpheres(); i++)
	{
		HitInfo resultSphere = currentScene.getSphere(i).intersects(ray);
		if (resultSphere.t != INFINITY && resultSphere.t < smallestT) {
			smallestT = resultSphere.t;
			smallestHit = resultSphere;
			closest_triangle = sphereFace;
			closest_triangle.material_id = currentScene.getSphere(i).getMaterialID();
		}
	}

	for (auto i = 0; i < boxes.size(); i++)
	{
		Box curr_box = boxes.at(i);
		HitInfo result_box = intersectBox(curr_box, ray.getOrigin(), ray.getDirection());

		if (result_box.t != INFINITY)
		{
			HitInfo result_triangle = intersectTriangle(curr_box.triangles, ray.getOrigin(), ray.getDirection());
			if (result_triangle.t != INFINITY && smallestT > result_triangle.t) {
				smallestT = result_triangle.t;
				smallestHit = result_triangle;
				closest_triangle = curr_box.triangles.at(result_triangle.faceId);
			}
		}
	}

	if (smallestT != INFINITY) {
		return Shader(bounce, closest_triangle, smallestHit, ray);
	}

	if (bounce == 0) return Eigen::Vector3f(1.0, 1.0, 1.0);
	return Eigen::Vector3f(0.0, 0.0, 0.0);
#else
	Tucano::Face closest_triangle;
	HitInfo smallestHit;
	float smallestT = INFINITY;
	std::queue<Box> queue_boxes;
	queue_boxes.push(root_box);

	while (queue_boxes.size() > 0) {
		Box box = queue_boxes.front();

		HitInfo result_box = intersectBox(box, ray.getOrigin(), ray.getDirection());

		if (result_box.t != INFINITY) {
			if (box.children.size() != 0)
			{
				queue_boxes.push(box.children.at(0));
				queue_boxes.push(box.children.at(1));
			}
			else {
				HitInfo result_triangle = intersectTriangle(box.triangles, ray.getOrigin(), ray.getDirection());
				if (result_triangle.t != INFINITY && smallestT > result_triangle.t) {
					smallestT = result_triangle.t;
					smallestHit = result_triangle;
					closest_triangle = box.triangles.at(result_triangle.faceId);
				}
			}
		}
		queue_boxes.pop();
	}

	//std::cout << "smallestT: " << smallestT << std::endl;

	//std::cout << "Box hit: " << i << std::endl;
	if (smallestT != INFINITY) {
		return Shader(bounce, closest_triangle, smallestHit, ray);
	}

	if (bounce == 0) return Eigen::Vector3f(1.0, 1.0, 1.0);

	return Eigen::Vector3f(0.0, 0.0, 0.0);
#endif
}

Eigen::Vector3f Flyscene::Shader(int bounce, Tucano::Face face, HitInfo hit, Ray ray) {

	//Tucano::Face face = mesh.getFace(hit.faceId);
	auto mat = materials[face.material_id];

	Eigen::Vector3f normalN = hit.normal;
	Eigen::Vector3f totalColor = Eigen::Vector3f(0.0, 0.0, 0.0);

	for (auto lightPosition : lights)
	{
		totalColor += calculateColor(bounce, lightPosition, hit, normalN, ray, mat); // Not sure what the reflection factor is. So any bugs could be caused by this
	}

	for (auto spherical : sphericalLights)
	{
		Eigen::Vector3f normal = hit.point - spherical.first;
		vector<Eigen::Vector3f> points = Shader::getNPointsOnCircle(spherical.first, spherical.second, normal.normalized(), 5);
		points.push_back(spherical.first);

		for (auto lightPosition : points)
		{
			totalColor += calculateColor(bounce, lightPosition, hit, normalN, ray, mat) / points.size(); // Not sure what the reflection factor is. So any bugs could be caused by this

		}
	}

	return totalColor;
}

Eigen::Vector3f Flyscene::calculateColor(int bounce, Eigen::Vector3f lightPosition, HitInfo hit, Eigen::Vector3f normalN, Ray ray, Tucano::Material::Mtl mat)
{
	// LIGHT
	int numLights = lights.size() + sphericalLights.size();
	Eigen::Vector3f lightIntensity = Eigen::Vector3f(1.0, 1.0, 1.0) / numLights;
	Eigen::Vector3f lightDirection = (lightPosition - hit.point).normalized();
	Eigen::Vector3f reflectedLight = Shader::reflect(-lightDirection, normalN);

	// EYE
	Eigen::Vector3f eyeDirection = (ray.getOrigin() - hit.point).normalized();
	float dotted = eyeDirection.dot(reflectedLight.normalized());

	Eigen::Vector3f reflectedColor = Eigen::Vector3f(0, 0, 0);
	Eigen::Vector3f refractedColor = Eigen::Vector3f(0, 0, 0);
	Eigen::Vector3f globalIllum = Eigen::Vector3f(0, 0, 0);

#define GLOBAL_RESOLUTION 4

	if (bounce < MAX_BOUNCES) {
		// calc reflectedRay
		Ray reflectedRay = ray.reflectRay(hit.normal, hit.point);
		reflectedColor = traceRay(bounce + 1, reflectedRay);

		Ray refractedRay = ray.refractRay(hit.normal, hit.point, -eyeDirection, mat.getOpticalDensity());
		refractedColor = traceRay(bounce + 1, refractedRay);

		/*std::vector<Ray> illum_rays = ray.resendRay(hit.normal, hit.point, GLOBAL_RESOLUTION);
		std::vector<Eigen::Vector3f> global_colors;

		for (size_t i = 0; i < illum_rays.size(); i++)
		{
			global_colors.push_back(traceRay(bounce + 1, illum_rays[i]));
		}

		for (size_t i = 0; i < global_colors.size(); i++)
		{
			globalIllum += global_colors[i];
		}
		globalIllum /= GLOBAL_RESOLUTION;*/
	}



	if (!canSeeLight(lightPosition, hit.point)) {
		return (mat.getAmbient().cwiseProduct(lightIntensity)) +mat.getSpecular().cwiseProduct(reflectedColor)) + mat.getDissolveFactor() * refractedColor;

	}

	Eigen::Vector3f ambient = mat.getAmbient().cwiseProduct(lightIntensity);

	Eigen::Vector3f diffuse = std::max(normalN.dot(lightDirection), 0.0f) * mat.getDiffuse().cwiseProduct(lightIntensity);

	Eigen::Vector3f specular = lightIntensity.cwiseProduct(mat.getSpecular()) * std::pow(std::max(dotted, 0.0f), 33);

	Eigen::Vector3f color = ambient + diffuse + specular;

	return color  +mat.getSpecular().cwiseProduct(reflectedColor) + mat.getDissolveFactor() * refractedColor; // Not sure what the reflection factor is. So any bugs could be caused by this
}

bool Flyscene::canSeeLight(Eigen::Vector3f lightPos, Eigen::Vector3f position)
{
	Eigen::Vector3f toLight = lightPos - position;
	Ray ray = Ray(position, toLight.normalized());
	HitInfo sphereHit;
	float toLightLength = toLight.norm();
	bool behind = false;
	bool sphereBlock = false;
	
	for (size_t i = 0; i < currentScene.getNumberOfSpheres(); i++)
	{
		sphereHit = currentScene.getSphere(i).intersects(ray);
		if (sphereHit.t != INFINITY && sphereHit.t > toLightLength) {
			behind = true;
		}
		else if (sphereHit.t != INFINITY){
			sphereBlock = true;
		}
	}

	HitInfo trianglHit = intersectTriangle(triangles, position, toLight);

	if (trianglHit.t != INFINITY && trianglHit.t > toLightLength)
	{
		//eventualHit = trianglHit;
		behind = true;
	}
	return ((!sphereBlock && trianglHit.t == INFINITY) || behind == true);
}

HitInfo Flyscene::intersectPlane(Eigen::Vector3f& origin, Eigen::Vector3f& dir) {

	int max = mesh.getNumberOfFaces();

	HitInfo smallesInfo = HitInfo{ INFINITY, -1 };

	// for all faces....
	for (size_t i = 0; i < max; i++)
	{
		Tucano::Face curFace = mesh.getFace(i);

		Eigen::Vector3f v0 = (mesh.getModelMatrix() * mesh.getVertex(curFace.vertex_ids.at(0))).head<3>(); // You can pick any vertex

		Eigen::Vector3f normalizeNormal = curFace.normal.normalized();

		std::vector<Eigen::Vector3f> temp{ v0 };

		Intersect::Face* plane = new Intersect::Plane(temp, normalizeNormal);

		Ray ray = Ray(origin, dir);

		HitInfo newInfo = plane->intersects(ray);

		if (newInfo.t < smallesInfo.t) {
			smallesInfo.t = newInfo.t;
			smallesInfo.faceId = i;
		}

		delete plane;
	}


	return smallesInfo;
}

HitInfo Flyscene::intersectTriangle(vector<Tucano::Face>& faces, Eigen::Vector3f origin, Eigen::Vector3f dir) {

	int max = faces.size();

	/*float smallestT = INFINITY;
	int smallestFace = -1;
	Eigen::Vector3f smallestHitPoint = Eigen::Vector3f(0, 0, 0);
	Eigen::Vector3f smallestNormal = Eigen::Vector3f(0, 0, 0);*/

	HitInfo smallestInfo = HitInfo { INFINITY, -1, Eigen::Vector3f(0, 0, 0),Eigen::Vector3f(0, 0, 0) };

	// for all faces....
	for (size_t i = 0; i < max; i++)
	{
		Tucano::Face curFace = faces.at(i);

		Eigen::Vector3f v0 = (mesh.getShapeModelMatrix() * mesh.getVertex(curFace.vertex_ids[0])).head<3>();
		Eigen::Vector3f v1 = (mesh.getShapeModelMatrix() * mesh.getVertex(curFace.vertex_ids[1])).head<3>();
		Eigen::Vector3f v2 = (mesh.getShapeModelMatrix() * mesh.getVertex(curFace.vertex_ids[2])).head<3>();

		Eigen::Vector3f normalizeNormal = (mesh.getShapeModelMatrix().inverse().matrix().transpose() * Eigen::Vector4f(curFace.normal.x(), curFace.normal.y(), curFace.normal.z(), 0)).head<3>().normalized();

		std::vector<Eigen::Vector3f> temp{ v0,v1,v2 };

		Intersect::Face * triangle = new Intersect::Triangle(temp, normalizeNormal);

		Ray ray = Ray(origin, dir);

		HitInfo newInfo = triangle-> intersects(ray);

		if (newInfo.t < smallestInfo.t) {
			smallestInfo.faceId = i;
			smallestInfo.t = newInfo.t;
			smallestInfo.normal = newInfo.normal;
			smallestInfo.point = newInfo.point;
		}
		delete( triangle);
		
		/*float D = normalizeNormal.dot(v0);

		float denom = normalizeNormal.dot(dir);

		if (denom != 0.0) {
			float t = (D - origin.dot(normalizeNormal)) / denom;

			Eigen::Vector3f hitPoint = origin + t * dir;

			if (t >= 0.001 && isInTriangle(hitPoint, v0, v1, v2)) {
				if (t < smallestT) {
					smallestT = t;
					smallestFace = i;
					smallestNormal = normalizeNormal;
					smallestHitPoint = hitPoint;
				}
			}
		}*/
		
	}

	return smallestInfo;//HitInfo { smallestT, smallestFace, smallestNormal, smallestHitPoint };
}

//
//bool Flyscene::isInTriangle(Eigen::Vector3f& hit, Eigen::Vector3f& v0, Eigen::Vector3f& v1, Eigen::Vector3f& v2)
//{
//	Eigen::Vector3f u = v1 - v0;
//	Eigen::Vector3f v = v2 - v0;
//	Eigen::Vector3f w = hit - v0;
//
//	float uDotU = u.dot(u);
//	float uDotV = u.dot(v);
//	float vDotV = v.dot(v);
//
//	float wDotU = w.dot(u);
//	float wDotV = w.dot(v);
//
//	float denom = uDotU * vDotV - uDotV * uDotV;
//	float s1 = vDotV * wDotU - uDotV * wDotV;
//	float t1 = uDotU * wDotV - uDotV * wDotU;
//
//	float s = s1 / denom;
//	float t = t1 / denom;
//
//	if (s < 0 || t < 0 || s + t > 1) return false;
//	return true;
//}

HitInfo Flyscene::intersectBox(Box& box, Eigen::Vector3f origin, Eigen::Vector3f dest) {

	/*Eigen::Vector3f dir = dest;
	Eigen::Vector3f invDir = Eigen::Vector3f(1 / dir.x(), 1 / dir.y(), 1 / dir.z());

	float tmin, tmax, tymin, tymax, tzmin, tzmax;

	if (invDir.x() >= 0) {
		tmin = (box.tmin.x() - origin.x()) * invDir.x();
		tmax = (box.tmax.x() - origin.x()) * invDir.x();
	}
	else {
		tmin = (box.tmax.x() - origin.x()) * invDir.x();
		tmax = (box.tmin.x() - origin.x()) * invDir.x();
	}

	if (invDir.y() >= 0) {
		tymin = (box.tmin.y() - origin.y()) * invDir.y();
		tymax = (box.tmax.y() - origin.y()) * invDir.y();
	}
	else {
		tymin = (box.tmax.y() - origin.y()) * invDir.y();
		tymax = (box.tmin.y() - origin.y()) * invDir.y();
	}
	if ((tmin > tymax) || (tymin > tmax)) {
		return HitInfo{ INFINITY, -1 };
	}
	if (tymin > tmin) {
		tmin = tymin;
	}
	if (tymax < tmax) {
		tmax = tymax;
	}

	if (invDir.z() >= 0) {
		tzmin = (box.tmin.z() - origin.z()) * invDir.z();
		tzmax = (box.tmax.z() - origin.z()) * invDir.z();
	}
	else {
		tzmin = (box.tmax.z() - origin.z()) * invDir.z();
		tzmax = (box.tmin.z() - origin.z()) * invDir.z();
	}

	if ((tmin > tzmax) || (tzmin > tmax)) {
		return HitInfo{ INFINITY, -1 };
	}
	if (tzmin > tmin) {
		tmin = tzmin;
	}
	if (tzmax < tmax) {
		tmax = tzmax;
	}

	return HitInfo{ tmin, -1 };*/

	

	Intersect::Face * boxie = new Intersect::BoxObject(box);

	Ray ray = Ray(origin, dest);

	HitInfo res = boxie->intersects(ray);

	delete( boxie);

	return res;
}

Eigen::Vector3f Flyscene::averagePoint(Box& box) {

	Eigen::Vector3f average_point;
	float sum_x = 0;
	float sum_y = 0;
	float sum_z = 0;

	for (auto i = 0; i < box.triangles.size(); i++)
	{
		for (auto n = 0; n < 3; n++)
		{
			Eigen::Vector3f v = (mesh.getShapeModelMatrix() * mesh.getVertex(box.triangles.at(i).vertex_ids[n])).head<3>();
			sum_x += v.x();
			sum_y += v.y();
			sum_z += v.z();
		}
	}
	float size = box.triangles.size() * 3;
	return Eigen::Vector3f(sum_x / size, sum_y / size, sum_z / size);
}

vector<Box> Flyscene::divideBox(Box& bigBox, int max_numberFaces) {


	vector<Box> result;
	std::queue<Box> list_box;
	list_box.push(bigBox);
	while (list_box.size() > 0)
	{
		Box box = list_box.front();

		if (box.triangles.size() <= max_numberFaces && box.triangles.size() > 0) {
			result.push_back(box);
			list_box.pop();

		}
		else
		{
			int axis = axisToDivide(box.tmax, box.tmin);
			//std::cout << "Number of boxes: " << result.size() << std::endl;

			Eigen::Vector3f average_point = averagePoint(box);

			Eigen::Vector3f midMax;
			Eigen::Vector3f midMin;

			if (axis == 0) {
				float x = average_point.x();

				midMax = Eigen::Vector3f(x, box.tmax.y(), box.tmax.z());
				midMin = Eigen::Vector3f(x, box.tmin.y(), box.tmin.z());

			}
			else if (axis == 1) {
				float y = average_point.y();

				midMax = Eigen::Vector3f(box.tmax.x(), y, box.tmax.z());
				midMin = Eigen::Vector3f(box.tmin.x(), y, box.tmin.z());
			}
			else {
				float z = average_point.z();

				midMax = Eigen::Vector3f(box.tmax.x(), box.tmax.y(), z);
				midMin = Eigen::Vector3f(box.tmin.x(), box.tmin.y(), z);

			}

			Box box1 = Box(box.tmin, midMax);
			Box box2 = Box(midMin, box.tmax);

			for (auto i = 0; i < box.triangles.size(); i++)
			{
				Tucano::Face face = box.triangles.at(i);
				bool isInBox1 = isInBox(box1, face);
				bool isInBox2 = isInBox(box2, face);

				// adding the triangle to the list of trinagles if it lies in box1 or box2
				if (isInBox1) {
					box1.triangles.push_back(face);
				}
				if (isInBox2) {
					box2.triangles.push_back(face);
				}
			}

			list_box.push(box1);
			list_box.push(box2);

			list_box.pop();

		}
	}
	std::cout << "Number of boxes: " << result.size() << std::endl;
	return result;
}

int Flyscene::axisToDivide(Eigen::Vector3f& tmax, Eigen::Vector3f& tmin) {
	float max = std::numeric_limits<float>::min();
	int result = -1;
	for (auto i = 0; i < 3; i++)
	{
		float diff = tmax[i] - tmin[i];
		if (max <= diff) {
			max = diff;
			result = i;
		}
	}
	return result;
}

bool Flyscene::isInBox(Box& box, Tucano::Face& face) {
	Eigen::Vector3f tmax = box.tmax;
	Eigen::Vector3f tmin = box.tmin;

	Eigen::Vector3f v0 = (mesh.getShapeModelMatrix() * mesh.getVertex(face.vertex_ids[0])).head<3>();
	Eigen::Vector3f v1 = (mesh.getShapeModelMatrix() * mesh.getVertex(face.vertex_ids[1])).head<3>();
	Eigen::Vector3f v2 = (mesh.getShapeModelMatrix() * mesh.getVertex(face.vertex_ids[2])).head<3>();

	vector<Eigen::Vector3f> vertices;
	vertices.push_back(v0);
	vertices.push_back(v1);
	vertices.push_back(v2);

	for (auto i = 0; i < vertices.size(); i++)
	{
		Eigen::Vector3f vertex = vertices.at(i);
		if (tmax.x() >= vertex.x() && tmin.x() <= vertex.x() && tmax.y() >= vertex.y()
			&& tmin.y() <= vertex.y() && tmax.z() >= vertex.z() && tmin.z() <= vertex.z())
		{
			return true;
		}
	}

	return false;
}

void Flyscene::divideBox_KD(int max_numberFaces) {

	std::queue<Box> list_box;
	list_box.push(root_box);
	int number_leafs = 0;
	while (list_box.size() > 0)
	{
		Box box = list_box.front();

		if (box.triangles.size() <= max_numberFaces && box.triangles.size() > 0) {
#ifdef show_KD
			this->boxes.push_back(box);
#endif // show_KD
			list_box.pop();
			number_leafs++;
		}
		else
		{
			int axis = axisToDivide(box.tmax, box.tmin);

			Eigen::Vector3f average_point = averagePoint(box);

			Eigen::Vector3f midMax;
			Eigen::Vector3f midMin;

			if (axis == 0) {
				float x = average_point.x();

				midMax = Eigen::Vector3f(x, box.tmax.y(), box.tmax.z());
				midMin = Eigen::Vector3f(x, box.tmin.y(), box.tmin.z());

			}
			else if (axis == 1) {
				float y = average_point.y();

				midMax = Eigen::Vector3f(box.tmax.x(), y, box.tmax.z());
				midMin = Eigen::Vector3f(box.tmin.x(), y, box.tmin.z());
			}
			else {
				float z = average_point.z();

				midMax = Eigen::Vector3f(box.tmax.x(), box.tmax.y(), z);
				midMin = Eigen::Vector3f(box.tmin.x(), box.tmin.y(), z);

			}

			Box box1 = Box(box.tmin, midMax);
			Box box2 = Box(midMin, box.tmax);

			// adding the children of the current box.
			box.children.push_back(box1);
			box.children.push_back(box1);

			for (auto i = 0; i < box.triangles.size(); i++)
			{
				Tucano::Face face = box.triangles.at(i);
				bool isInBox1 = isInBox(box1, face);
				bool isInBox2 = isInBox(box2, face);

				// adding the triangle to the list of trinagles if it lies in box1 or box2
				if (isInBox1) {
					box1.triangles.push_back(face);
				}
				if (isInBox2) {
					box2.triangles.push_back(face);
				}
			}

			list_box.push(box1);
			list_box.push(box2);

			list_box.pop();
		}
	}
	std::cout << "Number of leafs: " << number_leafs << std::endl;
}
