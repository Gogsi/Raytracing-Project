#include "flyscene.hpp"
#include <GLFW/glfw3.h>

#define MAX_BOUNCES 2

void Flyscene::initialize(int width, int height) {
  // initiliaze the Phong Shading effect for the Opengl Previewer
  phong.initialize();

  // set the camera's projection matrix
  flycamera.setPerspectiveMatrix(60.0, width / (float)height, 0.1f, 100.0f);
  flycamera.setViewport(Eigen::Vector2f((float)width, (float)height));

  // load the OBJ file and materials
  Tucano::MeshImporter::loadObjFile(mesh, materials,
                                    "resources/models/cube.obj");


  // normalize the model (scale to unit cube and center at origin)
  mesh.normalizeModelMatrix();

  // pass all the materials to the Phong Shader
  for (int i = 0; i < materials.size(); ++i)
    phong.addMaterial(materials[i]);



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
  int size = image_size[1];
  // for every pixel shoot a ray from the origin through the pixel coords
  for (int j = 0; j < image_size[1]; ++j) {
    for (int i = 0; i < image_size[0]; ++i) {
      // create a ray from the camera passing through the pixel (i,j)
      screen_coords = flycamera.screenToWorld(Eigen::Vector2f(i, j));
      // launch raytracing for the given ray and write result to pixel data
      pixel_data[i][j] = traceRay(0 , Ray(origin, screen_coords - origin));
    }
	std::cout << j*100/size << std::endl;
  }

  // write the ray tracing result to a PPM image
  Tucano::ImageImporter::writePPMImage("result.ppm", pixel_data);
  std::cout << "ray tracing done! " << std::endl;
}


Eigen::Vector3f Flyscene::traceRay(int bounce, Ray ray) {

	// "dest" is the location of the current pixel in world space. Subtracting camera origin from it gives the ray direction.
	Box box = Box(mesh);
	HitInfo result = intersectTriangle(ray.getOrigin(), ray.getDirection());

	if (result.t != INFINITY) {
		
		// recursive call
		return Shader(bounce, result, ray);
		//return Eigen::Vector3f(1.0, 0, 0);
	}

	return Eigen::Vector3f(0, 1.0, 0);
}

Eigen::Vector3f Flyscene::Shader(int bounce, HitInfo hit, Ray ray) {

	Tucano::Face face = mesh.getFace(hit.faceId);
	auto mat = materials[face.material_id];

	Eigen::Vector3f normalN = hit.normal;

	// LIGHT
	Eigen::Vector3f lightIntensity = Eigen::Vector3f(1.0, 1.0, 1.0);
	Eigen::Vector3f lightPosition = lightrep.getCentroid(); //for now
	Eigen::Vector3f lightDirection = (lightPosition - hit.point).normalized();
	Eigen::Vector3f reflectedLight = reflect(-lightDirection, normalN);

	// EYE
	Eigen::Vector3f eyeDirection = (ray.getOrigin() - hit.point).normalized();
	float dotted = eyeDirection.dot(reflectedLight.normalized());

	Eigen::Vector3f reflectedColor = Eigen::Vector3f(0, 0, 0);

	// RECURSION
	if (bounce < MAX_BOUNCES) {
		// calc reflectedRay
		Ray reflectedRay = ray.reflectRay(hit.normal, hit.point);
		reflectedColor = traceRay(bounce + 1, reflectedRay);
	}

	// AMBIENT DIFFUSE SPECULAR

	Eigen::Vector3f diffuse = std::max(normalN.dot(lightDirection), 0.0f) * mat.getDiffuse().cwiseProduct( lightIntensity);

	Eigen::Vector3f ambient = mat.getAmbient().cwiseProduct(lightIntensity);

	Eigen::Vector3f specular = multiply(lightIntensity, mat.getSpecular()) * std::pow(std::max(dotted, 0.0f), mat.getShininess());

	Eigen::Vector3f color = ambient + diffuse + specular;

	return color + mat.getSpecular().cwiseProduct(reflectedColor); // Not sure what the reflection factor is. So any bugs could be caused by this
	
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

HitInfo Flyscene::intersectPlane(Eigen::Vector3f& origin,
	Eigen::Vector3f& dir) {

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

			if (t >= 0) {
				if (t < smallestT) {
					smallestT = t;
					smallestFace = i;
				}
			}
		}
	}

	return HitInfo{ smallestT, smallestFace};
}

HitInfo Flyscene::intersectTriangle(Eigen::Vector3f origin,
	Eigen::Vector3f dir) {

	int max = mesh.getNumberOfFaces();

	float smallestT = INFINITY;
	int smallestFace = -1;
	Eigen::Vector3f smallestHitPoint = Eigen::Vector3f(0, 0, 0);
	Eigen::Vector3f smallestNormal = Eigen::Vector3f(0, 0, 0);

	// for all faces....
	for (size_t i = 0; i < max; i++)
	{
		Tucano::Face curFace = mesh.getFace(i);

		Eigen::Vector3f v0 = (mesh.getShapeModelMatrix() * mesh.getVertex(curFace.vertex_ids[0])).head<3>();
		Eigen::Vector3f v1 = (mesh.getShapeModelMatrix() * mesh.getVertex(curFace.vertex_ids[1])).head<3>();
		Eigen::Vector3f v2 = (mesh.getShapeModelMatrix() * mesh.getVertex(curFace.vertex_ids[2])).head<3>();

		Eigen::Vector3f normalizeNormal = curFace.normal.normalized();

		float D = normalizeNormal.dot(v0);

		float denom = normalizeNormal.dot(dir);

		if (denom != 0.0) {
			float t = (D - origin.dot(normalizeNormal)) / denom;

			Eigen::Vector3f hitPoint = origin + t * dir;

			if (t >= 0 && isInTriangle(hitPoint, v0, v1, v2)) {
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

HitInfo Flyscene::intersectBox(Box& box, Eigen::Vector3f& origin, Eigen::Vector3f& dest) {

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
		return HitInfo{ INFINITY, -1};
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
