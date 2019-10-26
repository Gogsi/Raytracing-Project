#include "flyscene.hpp"
#include <GLFW/glfw3.h>
#include <queue>
#include <thread>
#include <ctime>

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
                                    "resources/models/bunny.obj");


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

  // Flat structure:
  this->boxes = divideBox(root_box, 1000);
  
  // if u want to visualize the bounding boxeswith flat structure
   #define show_flat
   showBoxes();

  // if u want to visualize the bounding boxes with the KD trees 
  //#define show_KD
  //showBoxes();

}

void Flyscene::paintGL(void) {
  

  // update the camera view matrix with the last mouse interactions
  flycamera.updateViewMatrix();
  Eigen::Vector4f viewport = flycamera.getViewport();

  // clear the screen and set background color
  glClearColor(0.9, 0.9, 0.9, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // position the scene light at the last ray-tracing light source
  scene_light.resetViewMatrix();
  scene_light.viewMatrix()->translate(-lights.back());



  // render the scene using OpenGL and one light source
  phong.render(mesh, flycamera, scene_light);

  // render the bounding boxes:
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
	
  // render the ray and camera representation for ray debug
  ray.render(flycamera, scene_light);
  camerarep.render(flycamera, scene_light);

  // render ray tracing light sources as yellow spheres
  for (int i = 0; i < lights.size(); ++i) {
    lightrep.resetModelMatrix();
    lightrep.modelMatrix()->translate(lights[i]);
    lightrep.render(flycamera, scene_light);
  }

  // render coordinate system at lower right corner
  flycamera.renderAtCorner();
}

void Flyscene::simulate(GLFWwindow *window) {
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

void Flyscene::createDebugRay(const Eigen::Vector2f &mouse_pos) {
  ray.resetModelMatrix();
  // from pixel position to world coordinates
  Eigen::Vector3f screen_pos = flycamera.screenToWorld(mouse_pos);

  // direction from camera center to click position
  Eigen::Vector3f dir = (screen_pos - flycamera.getCenter()).normalized();
  
  // position and orient the cylinder representing the ray
  ray.setOriginOrientation(flycamera.getCenter(), dir);

  // place the camera representation (frustum) on current camera location, 
  camerarep.resetModelMatrix();
  camerarep.setModelMatrix(flycamera.getViewMatrix().inverse());

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
  std::cout << "ray tracing done! " << std::endl;
}

void Flyscene::updating_pixels(vector<vector<Eigen::Vector3f>>& pixel_data, Eigen::Vector3f& origin, Eigen::Vector2i& image_size, int number_threads, int thread_id) {
	// for every pixel shoot a ray from the origin through the pixel coords
	for (int j = thread_id; j < image_size[1]; j += number_threads) {
		for (int i = 0; i < image_size[0]; i++) {
			// create a ray from the camera passing through the pixel (i,j)
			Eigen::Vector3f screen_coords = flycamera.screenToWorld(Eigen::Vector2f(i, j));
			// launch raytracing for the given ray and write result to pixel data
			pixel_data[i][j] = traceRay(0, Ray(origin, screen_coords - origin));
		}
	}
}

void Flyscene::updating_pixels_KD(vector<vector<Eigen::Vector3f>>& pixel_data, Eigen::Vector3f& origin, Eigen::Vector2i& image_size, int number_threads, int thread_id) {
	// for every pixel shoot a ray from the origin through the pixel coords
	for (int j = thread_id; j < image_size[1]; j += number_threads) {
		for (int i = 0; i < image_size[0]; i++) {
			// create a ray from the camera passing through the pixel (i,j)
			Eigen::Vector3f screen_coords = flycamera.screenToWorld(Eigen::Vector2f(i, j));
			// launch raytracing for the given ray and write result to pixel data
			pixel_data[i][j] = traceRay_KD(root_box ,0 ,Ray(origin, screen_coords - origin));
		}
	}
}


void Flyscene::showBoxes() {

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
	
	//std::cout << boxes.size() << std::endl;
	Tucano::Face closest_triangle;
	HitInfo smallestHit;
	float smallestT = INFINITY;
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
	//std::cout << "Box hit: " << i << std::endl;
	if (smallestT != INFINITY) {
		return Shader(bounce, closest_triangle, smallestHit, ray);
	}

	if(bounce == 0) return Eigen::Vector3f(1.0, 1.0, 1.0);
	return Eigen::Vector3f(0.0, 0.0, 0.0);
}

Eigen::Vector3f Flyscene::Shader(int bounce, Tucano::Face face, HitInfo hit, Ray ray) {

	//Tucano::Face face = mesh.getFace(hit.faceId);
	auto mat = materials[face.material_id];

	Eigen::Vector3f normalN = hit.normal;
	Eigen::Vector3f totalColor = Eigen::Vector3f(0.0, 0.0, 0.0);

	for (auto lightPosition : lights)
	{
		// LIGHT
		Eigen::Vector3f lightIntensity = Eigen::Vector3f(1.0, 1.0, 1.0);
		Eigen::Vector3f lightDirection = (lightPosition - hit.point).normalized();
		Eigen::Vector3f reflectedLight = reflect(-lightDirection, normalN);

		// EYE
		Eigen::Vector3f eyeDirection = (ray.getOrigin() - hit.point).normalized();
		float dotted = eyeDirection.dot(reflectedLight.normalized());

		Eigen::Vector3f reflectedColor = Eigen::Vector3f(0, 0, 0);


		if (bounce < MAX_BOUNCES) {
			// calc reflectedRay
			Ray reflectedRay = ray.reflectRay(hit.normal, hit.point);
			reflectedColor = traceRay(bounce + 1, reflectedRay);
		}

		if (!canSeeLight(lightPosition, hit.point)) {
			totalColor += mat.getAmbient().cwiseProduct(lightIntensity) + mat.getSpecular().cwiseProduct(reflectedColor);
			continue;
		}



		Eigen::Vector3f ambient = mat.getAmbient().cwiseProduct(lightIntensity);

		Eigen::Vector3f diffuse = std::max(normalN.dot(lightDirection), 0.0f) * mat.getDiffuse().cwiseProduct(lightIntensity);

		Eigen::Vector3f specular = multiply(lightIntensity, mat.getSpecular()) * std::pow(std::max(dotted, 0.0f), 33);

		Eigen::Vector3f color = ambient + diffuse + specular;

		totalColor += color + mat.getSpecular().cwiseProduct(reflectedColor); // Not sure what the reflection factor is. So any bugs could be caused by this
	}
	
	return totalColor;
}

Eigen::Vector3f Flyscene::reflect(Eigen::Vector3f light, Eigen::Vector3f normal) {
	return  light - 2 * normal.dot(light) * normal;
}

Eigen::Vector3f Flyscene::multiply(Eigen::Vector3f a, Eigen::Vector3f b) {
	float x = a.x() * b.x();
	float y = a.y() * b.y();
	float z = a.z() * b.z();
	return Eigen::Vector3f(x,y,x);
}

bool Flyscene::canSeeLight(Eigen::Vector3f lightPos, Eigen::Vector3f position)
{
	Eigen::Vector3f toLight = lightPos - position;
	HitInfo hit = intersectTriangle(triangles, position, toLight);
	float toLightLength = toLight.norm();
	bool behind = false;

	if (hit.t != INFINITY)
	{

		if (hit.t > toLightLength)
		{
			behind = true;
		}

	}
	return (hit.t == INFINITY || behind == true);
}

vector<Eigen::Vector3f> Flyscene::getNPointsOnCircle(Eigen::Vector3f center, float radius, Eigen::Vector3f normal, int n)
{
	Eigen::Vector3f notNormal = Eigen::Vector3f(1.0, 0.0, 0.0);
	if (notNormal == normal || notNormal == -normal) {
		notNormal = Eigen::Vector3f(1.0, 1.0, 0.0).normalized();
	}

	Eigen::Vector3f radiusVector = normal.cross(notNormal).normalized() * radius;
	float theta = M_PI * 2 / n;
	Eigen::AngleAxisf rotation = Eigen::AngleAxisf(theta, normal);

	vector<Eigen::Vector3f> res;
	res.push_back(radiusVector);

	for (int i = 0; i < n - 1; i++) {
		radiusVector = rotation * radiusVector;
		res.push_back(radiusVector);
	}

	return res;
}

HitInfo Flyscene::intersectPlane(Eigen::Vector3f& origin, Eigen::Vector3f& dir) {

	int max = mesh.getNumberOfFaces();

	float smallestT = INFINITY;
	int smallestFace = -1;
	// for all faces....
	for (size_t i = 0; i < max; i++)
	{
		Tucano::Face curFace = mesh.getFace(i);

		Eigen::Vector3f v0 = (mesh.getModelMatrix() * mesh.getVertex(curFace.vertex_ids.at(0))).head<3>(); // You can pick any vertex

		Eigen::Vector3f normalizeNormal = curFace.normal.normalized();

		auto D = normalizeNormal.dot(v0); // any point on the plane can be used to calculate D
		float denom = normalizeNormal.dot(dir);

		if (denom != 0.0) {
			float t = (D - origin.dot(normalizeNormal)) / denom;

			Eigen::Vector3f hit = origin + t * dir;

			if (t >= 1E-6) {
				if (t < smallestT) {
					smallestT = t;
					smallestFace = i;
				}
			}
		}
	}

	return HitInfo{ smallestT, smallestFace};
}

HitInfo Flyscene::intersectTriangle(vector<Tucano::Face>& faces, Eigen::Vector3f origin, Eigen::Vector3f dir) {

	int max = faces.size();

	float smallestT = INFINITY;
	int smallestFace = -1;
	Eigen::Vector3f smallestHitPoint = Eigen::Vector3f(0, 0, 0);
	Eigen::Vector3f smallestNormal = Eigen::Vector3f(0, 0, 0);

	// for all faces....
	for (size_t i = 0; i < max; i++)
	{
		Tucano::Face curFace = faces.at(i);

		Eigen::Vector3f v0 = (mesh.getShapeModelMatrix() * mesh.getVertex(curFace.vertex_ids[0])).head<3>();
		Eigen::Vector3f v1 = (mesh.getShapeModelMatrix() * mesh.getVertex(curFace.vertex_ids[1])).head<3>();
		Eigen::Vector3f v2 = (mesh.getShapeModelMatrix() * mesh.getVertex(curFace.vertex_ids[2])).head<3>();

		Eigen::Vector3f normalizeNormal = (mesh.getShapeModelMatrix().inverse().matrix().transpose() * Eigen::Vector4f(curFace.normal.x(), curFace.normal.y(), curFace.normal.z(), 0)).head<3>().normalized();

		float D = normalizeNormal.dot(v0);

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
		}
	}

	return HitInfo { smallestT, smallestFace, smallestNormal, smallestHitPoint };
}

HitInfo Flyscene::intersectBox(Box& box, Eigen::Vector3f origin, Eigen::Vector3f dest) {

	Eigen::Vector3f dir = dest;
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

	return HitInfo{ tmin, -1 };
}

bool Flyscene::isInTriangle(Eigen::Vector3f& hit, Eigen::Vector3f& v0, Eigen::Vector3f& v1, Eigen::Vector3f& v2)
{
	Eigen::Vector3f u = v1 - v0;
	Eigen::Vector3f v = v2 - v0;
	Eigen::Vector3f w = hit - v0;

	float uDotU = u.dot(u);
	float uDotV = u.dot(v);
	float vDotV = v.dot(v);

	float wDotU = w.dot(u);
	float wDotV = w.dot(v);

	float denom = uDotU * vDotV - uDotV * uDotV;
	float s1 = vDotV * wDotU - uDotV * wDotV;
	float t1 = uDotU * wDotV - uDotV * wDotU;

	float s = s1 / denom;
	float t = t1 / denom;

	if (s < 0 || t < 0 || s + t > 1) return false;
	return true;
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
	while (list_box.size() > 0 )
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

Eigen::Vector3f Flyscene::traceRay_KD(Box& big_box, int bounce, Ray ray) {


	
	Tucano::Face closest_triangle;
	HitInfo smallestHit;
	float smallestT = INFINITY;
	std::queue<Box> queue_boxes;
	queue_boxes.push(big_box);

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
}
