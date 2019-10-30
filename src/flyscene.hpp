#ifndef __FLYSCENE__
#define __FLYSCENE__

// Must be included before glfw.
#include <GL/glew.h>

#include <GLFW/glfw3.h>

#include <tucano/effects/phongmaterialshader.hpp>
#include <tucano/mesh.hpp>
#include <tucano/shapes/camerarep.hpp>
#include <tucano/shapes/cylinder.hpp>
#include <tucano/shapes/sphere.hpp>
#include <tucano/shapes/box.hpp>
#include <tucano/utils/flycamera.hpp>
#include <tucano/utils/imageIO.hpp>
#include <tucano/utils/mtlIO.hpp>
#include <tucano/utils/objimporter.hpp>
#include <thread> 

#include "hitInfo.h"
#include "Ray.h"
#include "box.h"
#include "Intersect.h"
#include "Scene.h"

class Flyscene {

public:
	Flyscene(void) { this->root_box = Box(); }

  /**
   * @brief Initializes the shader effect
   * @param width Window width in pixels
   * @param height Window height in pixels
   */
  void initialize(int width, int height);

  /**
   * Repaints screen buffer.
   **/
  virtual void paintGL();

  /**
   * Perform a single simulation step.
   **/
  virtual void simulate(GLFWwindow *window);

  /**
   * Returns the pointer to the flycamera instance
   * @return pointer to flycamera
   **/
  Tucano::Flycamera *getCamera(void) { return &flycamera; }

  /**
   * @brief Add a new light source
   */
  void addLight(void) { lights.push_back(flycamera.getCenter()); }
  void addSphericalLight(void) { sphericalLights.push_back(make_pair(flycamera.getCenter(), 0.01)); }


  /**
   * @brief Create a debug ray at the current camera location and passing
   * through pixel that mouse is over
   * @param mouse_pos Mouse cursor position in pixels
   */
  void createDebugRay(const Eigen::Vector2f &mouse_pos);

  /**
   * @brief raytrace your scene from current camera position   
   */
  void raytraceScene(int width = 0, int height = 0);

  /**
   * @brief trace a single ray from the camera passing through dest
   * @param origin Ray origin
   * @param dest Other point on the ray, usually screen coordinates
   * @return a RGB color
   */
  Eigen::Vector3f traceRay(int bounce, Ray ray, bool insideObject);

  HitInfo intersectPlane(Eigen::Vector3f& origin, Eigen::Vector3f& dir);

  HitInfo intersectTriangle(vector<Tucano::Face>& faces, Eigen::Vector3f origin, Eigen::Vector3f dir);

  HitInfo intersectBox(Box& box, Eigen::Vector3f origin, Eigen::Vector3f dest);

  bool isInTriangle(Eigen::Vector3f& hit, Eigen::Vector3f& v0, Eigen::Vector3f& v1, Eigen::Vector3f& v2);

  vector<Box> divideBox(Box& box, int max_numberFaces);

  int axisToDivide(Eigen::Vector3f& tmax, Eigen::Vector3f& tmin);

  bool isInBox(Box& box, Tucano::Face& face);

  Eigen::Vector3f averagePoint(Box& box);

  Eigen::Vector3f Shader(int bounce, Tucano::Face face, HitInfo hit, Ray ray, bool insideObject);

  bool canSeeLight(Eigen::Vector3f lightPos, Eigen::Vector3f position);

  void updating_pixels(vector<vector<Eigen::Vector3f>>& pixel_data, Eigen::Vector3f& origin, Eigen::Vector2i& image_size, int number_threads, int thread_id);

  void divideBox_KD(int max_numberFaces);

  Eigen::Vector3f calculateColor(int bounce, Eigen::Vector3f lightPosition, HitInfo hit, Eigen::Vector3f normalN, Ray ray, Tucano::Material::Mtl mat, bool insideObject);

  void initBoundingBoxes();

  void ReflectDebugRay(Eigen::Vector3f origin, Eigen::Vector3f dir, int bounce, int jumps);

  void renderBoundingBoxes();

  std::pair<HitInfo,Tucano::Face> getIntersections(Ray ray);
private:
  // A simple phong shader for rendering meshes
  Tucano::Effects::PhongMaterial phong;

  int jumps;

  // A fly through camera
  Tucano::Flycamera flycamera;

  // the size of the image generated by ray tracing
  Eigen::Vector2i raytracing_image_size;

  // A camera representation for animating path (false means that we do not
  // render front face)
  Tucano::Shapes::CameraRep camerarep = Tucano::Shapes::CameraRep(false);

  // a frustum to represent the camera in the scene
  Tucano::Shapes::Sphere lightrep;

  // light sources for ray tracing
  vector<Eigen::Vector3f> lights;
  vector<pair<Eigen::Vector3f, float>> sphericalLights;

  // Scene light represented as a camera
  Tucano::Camera scene_light;

  /// A very thin cylinder to draw a debug ray
  Tucano::Shapes::Cylinder ray = Tucano::Shapes::Cylinder(0.1, 1.0, 16, 64);
  Tucano::Shapes::Cylinder lightRay = Tucano::Shapes::Cylinder(0.1, 1.0, 16, 64);

  // List of bounding boxes to show:
  vector<Tucano::Shapes::Box> bounding_boxes;

  // Scene meshes
  Tucano::Mesh mesh;

  /// MTL materials
  vector<Tucano::Material::Mtl> materials;

  // list of bounding boxes
  vector<Box> boxes;

  vector<Tucano::Face> triangles;

  vector<Tucano::Shapes::Cylinder> rays;

  vector<Tucano::Shapes::Box> ray_hitbox;

  Tucano::Face sphereFace;

  Scene currentScene = Scene();

  bool renderBox;

  bool renderAllBox;


public:
  // Root box
  Box root_box;

};

#endif // FLYSCENE
