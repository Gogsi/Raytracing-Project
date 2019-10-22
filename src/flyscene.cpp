#include "flyscene.hpp"
#include <GLFW/glfw3.h>
#include <queue>

using namespace std;

void Flyscene::initialize(int width, int height) {
  // initiliaze the Phong Shading effect for the Opengl Previewer
  phong.initialize();

  // set the camera's projection matrix
  flycamera.setPerspectiveMatrix(60.0, width / (float)height, 0.1f, 100.0f);
  flycamera.setViewport(Eigen::Vector2f((float)width, (float)height));

  // load the OBJ file and materials
  Tucano::MeshImporter::loadObjFile(mesh, materials,
                                    "resources/models/dodgeColorTest.obj");


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

  // create the array of boxes
  Box box = Box(mesh);
  this->boxes = divideBox(box, 1000);

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

  std::cout << boxes.size() << std::endl;

  // for every pixel shoot a ray from the origin through the pixel coords
  for (int j = 0; j < image_size[1]; ++j) {
    for (int i = 0; i < image_size[0]; ++i) {
      // create a ray from the camera passing through the pixel (i,j)
      screen_coords = flycamera.screenToWorld(Eigen::Vector2f(i, j));
      // launch raytracing for the given ray and write result to pixel data
      pixel_data[i][j] = traceRay(origin, screen_coords);
    }
  }

  // write the ray tracing result to a PPM image
  Tucano::ImageImporter::writePPMImage("result.ppm", pixel_data);
  std::cout << "ray tracing done! " << std::endl;
}


Eigen::Vector3f Flyscene::traceRay(Eigen::Vector3f &origin,
                                   Eigen::Vector3f &dest) {
	Eigen::Vector3f newOrigin = origin;

	// "dest" is the location of the current pixel in world space. Subtracting camera origin from it gives the ray direction.
	Eigen::Vector3f newDir = dest - origin; 
	

	//std::cout << boxes.size() << std::endl;

	for (auto i = 0; i < boxes.size() ; i++)
	{
		HitInfo result = intersectBox(boxes.at(i), origin, newDir);

		if (result.t != INFINITY)
		{
			std::cout << "Box ID:" + i << std::endl;
			std::cout << "Number of faces" + boxes.at(i).triangles.size() << std::endl;
			return Eigen::Vector3f(i/30, 1.0, i/30);
		}
	}

/*	if (result.t != INFINITY) {
		Tucano::Face face = mesh.getFace(result.faceId);
		auto mat = phong.getMaterial(face.material_id);
		return mat.getDiffuse();
	}*/

	return Eigen::Vector3f(1.0, 0, 0);
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

HitInfo Flyscene::intersectTriangle(Eigen::Vector3f& origin,
	Eigen::Vector3f& dir) {

	int max = mesh.getNumberOfFaces();

	float smallestT = INFINITY;
	int smallestFace = -1;

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

			Eigen::Vector3f hit = origin + t * dir;

			if (t >= 0 && isInTriangle(hit, v0, v1, v2)) {
				if (t < smallestT) {
					smallestT = t;
					smallestFace = i;
				}
			}
		}
	}

	return HitInfo { smallestT, smallestFace };
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
	int n = 10;
	while (list_box.size() > 0 )
	{
		/*std::cout << "size of the queue before taking the first element" << std::endl;
		std::cout << list_box.size() << std::endl;*/
		Box box = list_box.front();
		std::cout << " " << std::endl;
		std::cout << "number of triangles : "<< box.triangles.size() << std::endl;
		
		if (box.triangles.size() <= max_numberFaces && box.triangles.size() > 0) {
			result.push_back(box);
			list_box.pop();
			std::cout << "box added to the final result"  << std::endl;
			std::cout << result.size() << std::endl;
			
		}
		else 
		{
			int axis = axisToDivide(box.tmax, box.tmin);
			std::cout << "division on axis : " << axis << std::endl;

			Eigen::Vector3f average_point = averagePoint(box);
			std::cout << "average_point : " << average_point << std::endl;

			if (axis == 0) {
				float x = average_point.x();

				Eigen::Vector3f midMax = Eigen::Vector3f(x, box.tmax.y(), box.tmax.z());
				Eigen::Vector3f midMin = Eigen::Vector3f(x, box.tmin.y(), box.tmin.z());

				Box box1 = Box(box.tmin, midMax);
				Box box2 = Box(midMin, box.tmax);

				for (auto i = 0; i < box.triangles.size(); i++)
				{
					Tucano::Face face = box.triangles.at(i);
					bool isInBox1 = isInBox(box1, face);
					bool isInBox2 = isInBox(box2, face);
					if (isInBox1) {
						box1.triangles.push_back(face);
					}
					if (isInBox2) {
						box2.triangles.push_back(face);
					}
				}

				std::cout << "number of triangles of box1: " << box1.triangles.size() << std::endl;
				std::cout << "number of triangles of box2: " << box2.triangles.size() << std::endl;

				list_box.push(box1);
				list_box.push(box2);
				/*result.push_back(box1);
				result.push_back(box2);*/
			}
			else if (axis == 1) {
				float y = average_point.y();

				Eigen::Vector3f midMax = Eigen::Vector3f(box.tmax.x(), y, box.tmax.z());
				Eigen::Vector3f midMin = Eigen::Vector3f(box.tmin.x(), y, box.tmin.z());

				Box box1 = Box(box.tmin, midMax);
				Box box2 = Box(midMin, box.tmax);

				for (auto i = 0; i < box.triangles.size(); i++)
				{
					Tucano::Face face = box.triangles.at(i);
					bool isInBox1 = isInBox(box1, face);
					bool isInBox2 = isInBox(box2, face);
					if (isInBox1) {
						box1.triangles.push_back(face);
					}
					if (isInBox2) {
						box2.triangles.push_back(face);
					}
				}

				std::cout << "number of triangles of box1: " << box1.triangles.size() << std::endl;
				std::cout << "number of triangles of box2: " << box2.triangles.size() << std::endl;

				list_box.push(box1);
				list_box.push(box2);
				/*result.push_back(box1);
				result.push_back(box2);*/
			}
			else {
				float z = average_point.z();

				Eigen::Vector3f midMax = Eigen::Vector3f(box.tmax.x(), box.tmax.y(), z);
				Eigen::Vector3f midMin = Eigen::Vector3f(box.tmin.x(), box.tmin.y(), z);

				Box box1 = Box(box.tmin, midMax);
				Box box2 = Box(midMin, box.tmax);

				for (auto i = 0; i < box.triangles.size(); i++)
				{
					Tucano::Face face = box.triangles.at(i);
					bool isInBox1 = isInBox(box1, face);
					bool isInBox2 = isInBox(box2, face);
					if (isInBox1) {
						box1.triangles.push_back(face);
					}
					if (isInBox2) {
						box2.triangles.push_back(face);
					}
				}

				std::cout << "number of triangles of box1: " << box1.triangles.size() << std::endl;
				std::cout << "number of triangles of box2: " << box2.triangles.size() << std::endl;

				list_box.push(box1);
				list_box.push(box2);
				/*result.push_back(box1);
				result.push_back(box2);*/
			}
			std::cout << "size of the queue after adding wo boxes" << std::endl;
			std::cout << list_box.size() << std::endl;
			list_box.pop();
			std::cout << "size of the queue after poping" << std::endl;
			std::cout << list_box.size() << std::endl;
		}
	}
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
