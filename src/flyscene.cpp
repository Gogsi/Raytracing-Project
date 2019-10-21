#include "flyscene.hpp"
#include <GLFW/glfw3.h>

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
  float dx = (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS ? 1.0 : 0.0) -
             (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS ? 1.0 : 0.0);
  float dy = (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS
                  ? 1.0
                  : 0.0) -
             (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS
                  ? 1.0
                  : 0.0);
  float dz = (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS ? 1.0 : 0.0) -
             (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS ? 1.0 : 0.0);
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
  float size = image_size[1];
  // for every pixel shoot a ray from the origin through the pixel coords
  for (int j = 0; j < image_size[1]; ++j) {
    for (int i = 0; i < image_size[0]; ++i) {
      // create a ray from the camera passing through the pixel (i,j)
      screen_coords = flycamera.screenToWorld(Eigen::Vector2f(i, j));
      // launch raytracing for the given ray and write result to pixel data
      pixel_data[i][j] = traceRay(origin, screen_coords);
    }
	std::cout << j*100/size << std::endl;
  }

  // write the ray tracing result to a PPM image
  Tucano::ImageImporter::writePPMImage("result.ppm", pixel_data);
  std::cout << "ray tracing done! " << std::endl;
}

class Box {
public:
	Eigen::Vector3f tmax, tmin;
	Box(Eigen::Vector3f tmin, Eigen::Vector3f tmax) {
		this->tmax = tmax;
		this->tmin = tmin;
	}
	Box(Tucano::Mesh mesh, int start, int end) {

		float tx_max = std::numeric_limits<float>::min();
		float ty_max = std::numeric_limits<float>::min();
		float tz_max = std::numeric_limits<float>::min();
		float tx_min = std::numeric_limits<float>::max();
		float ty_min = std::numeric_limits<float>::max();
		float tz_min = std::numeric_limits<float>::max();

		for (int i = start; i < end; i++) {
			Eigen::Vector4f v = mesh.getVertex(i);

			float x = v.x();
			float y = v.y();
			float z = v.z();

			if (x > tx_max) {
				tx_max = x;
			}
			if (x < tx_min) {
				tx_min = x;
			}
			if (y > ty_max) {
				ty_max = y;
			}
			if (y < ty_min) {
				ty_min = y;
			}
			if (z > tz_max) {
				tz_max = z;
			}
			if (z < tz_min) {
				tz_min = z;
			}
		};



		tmax = Eigen::Vector3f(tx_max, ty_max, tz_max);
		tmin = Eigen::Vector3f(tx_min, ty_min, tz_min);

		/*std::cout << tmax << std::endl;
		std::cout << tmin << std::endl;*/

	};
};



bool intersectionBox(Box& box, Eigen::Vector3f& origin, Eigen::Vector3f& dest) {

	Eigen::Vector3f dir = dest - origin;
	Eigen::Vector3f invDir = Eigen::Vector3f(1/dir.x(), 1/dir.y(), 1/dir.z());

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
		return false;
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
		return false;
	}
	if (tzmin > tmin) {
		tmin = tzmin;
	}
	if (tzmax < tmax) {
		tmax = tzmax;
	}

	return true;





	/*
	float tx_min = (box.tmin.x() - origin.x()) / dir.x();
	float tx_max = (box.tmax.x() - origin.x()) / dir.x();
	float ty_min = (box.tmin.y() - origin.y()) / dir.y();
	float ty_max = (box.tmax.y() - origin.y()) / dir.y();
	float tz_min = (box.tmin.z() - origin.z()) / dir.z();
	float tz_max = (box.tmax.z() - origin.z()) / dir.z();

	float tx_in = std::min(tx_min, tx_max);
	float tx_out = std::max(tx_min, tx_max);
	float ty_in = std::min(ty_min, ty_max);
	float ty_out = std::max(ty_min, ty_max);
	float tz_in = std::min(tz_min, tz_max);
	float tz_out = std::max(tz_min, tz_max);

	float t_in = std::max(tx_in, ty_in);
	float t_win = std::max(t_in, tz_in);

	float t_out = std::min(tx_out, tx_out);
	float t_wout = std::min(t_out, tz_out);


	if (t_win > t_wout || t_wout < 0) {
	
		return false;

	}

	std::cout << t_win << std::endl;
	// std::cout << t_wout << std::endl;

	return true;
	*/
};

Eigen::Vector3f Flyscene::traceRay(Eigen::Vector3f& origin,
	Eigen::Vector3f& dest) {
	Box box = Box(mesh, 0, mesh.getNumberOfVertices() / 2);

	bool ray_box = intersectionBox(box, origin, dest);

	if (ray_box) {
		//std::cout << "hit" << std::endl;
		return Eigen::Vector3f(1.0, 0.0, 0.0);
	}

	return Eigen::Vector3f(0.0, 1.0, 0.0);

	//// just some fake random color per pixel until you implement your ray tracing
	//// remember to return your RGB values as floats in the range [0, 1]!!!
	//return Eigen::Vector3f(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX,
	//	rand() / (float)RAND_MAX);
}