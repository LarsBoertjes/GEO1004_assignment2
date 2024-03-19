#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

//-- https://github.com/nlohmann/json
//-- used to read and write (City)JSON
#include "json.hpp" //-- it is in the /include/ folder

// CGAL libraries
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/linear_least_squares_fitting_3.h>

using json = nlohmann::json;


int   get_no_roof_surfaces(json &j);
int   get_no_ground_surfaces(json &j);
void  visit_roofsurfaces(json &j);
void  list_all_vertices(json &j);

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point3;
typedef Kernel::Plane_3 Plane3;
typedef Kernel::Vector_3 Vector3;


int main(int argc, const char * argv[]) {
  //-- will read the file passed as argument or twobuildings.city.json if nothing is passed
  const char* filename = (argc > 1) ? argv[1] : "../data/twobuildings.city.json";
  std::cout << "Processing: " << filename << std::endl;
  std::ifstream input(filename);
  json j;
  input >> j; //-- store the content of the file in a nlohmann::json object
  input.close();

  for (auto& co : j["CityObjects"].items()) {
      std::vector<std::vector<int>> groundSurfaceBoundaries;
      std::vector<std::vector<int>> roofSurfaceBoundaries;

      for (auto& geom : co.value()["geometry"]) {
          if (geom.contains("boundaries")) {
              for (auto& shell : geom["boundaries"]) {
                  for (auto& surface : shell) {
                      std::vector<Point3> vertex_coord;
                      for (auto& boundary : surface) {
                          for (auto& v : boundary) {
                              std::vector<int> vi = j["vertices"][v.get<int>()];
                              double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
                              double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
                              double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                              vertex_coord.emplace_back(x, y, z);
                          }
                      }

                      Plane3 plane;
                      CGAL::linear_least_squares_fitting_3(vertex_coord.begin(), vertex_coord.end(), plane, CGAL::Dimension_tag<0>());

                      Vector3 normal = plane.orthogonal_vector();
                      // Normalizing the vector for convenience?
                      //normal = normal / std::sqrt(normal.squared_length());

                      std::cout << std::setprecision(2) << std::fixed << "Normal Vector: (" << normal.x() << ", " << normal.y() << ", " << normal.z() << ")" << std::endl;
                  }
              }
          }
      }
      std::cout << "------------------------------------------------------" << std::endl;
  }

  for (auto& co : j["CityObjects"].items()) {
      if (co.value().contains("geometry") && !co.value()["geometry"].empty()) {
          std::vector<std::vector<int>> groundSurfaceVertices;
          for (auto& geom : co.value()["geometry"]) {
              if (geom["type"] == "Solid" || geom["type"] == "MultiSurface") {
                  if (geom.contains("semantics") && geom["semantics"].contains("values")) {
                      for (size_t surfaceIndex = 0; surfaceIndex < geom["semantics"]["values"].size(); ++surfaceIndex) {
                          for (size_t ringIndex = 0; ringIndex < geom["semantics"]["values"][surfaceIndex].size(); ++ringIndex) {
                              int sem_index = geom["semantics"]["values"][surfaceIndex][ringIndex];
                              if (geom["semantics"]["surfaces"][sem_index]["type"] == "GroundSurface") {
                                  std::vector<int> surfaceVertices;
                                  for (auto& index : geom["boundaries"][surfaceIndex][ringIndex][0]) {
                                      surfaceVertices.push_back(index);
                                  }
                                  groundSurfaceVertices.push_back(surfaceVertices);
                              }
                          }
                      }
                  }
              }
          }
          // Create a new MultiSurface geometry if there are ground surfaces
          if (!groundSurfaceVertices.empty()) {
              json newGeometry = {
                      {"type", "MultiSurface"},
                      {"lod", "0.2"},
                      {"boundaries", groundSurfaceVertices}
              };
              // Append the new geometry to the CityObject's geometry list
              co.value()["geometry"].push_back(newGeometry);
          }
      }
  }



  //-- write to disk the modified city model (out.city.json)
  std::ofstream o("out.city.json");
  o << j.dump(5) << std::endl;
  o.close();


  return 0;
}


// Visit every 'RoofSurface' in the CityJSON model and output its geometry (the arrays of indices)
// Useful to learn to visit the geometry boundaries and at the same time check their semantics.
void visit_roofsurfaces(json &j) {
  for (auto& co : j["CityObjects"].items()) {
    for (auto& g : co.value()["geometry"]) {
      if (g["type"] == "Solid") {
        for (int i = 0; i < g["boundaries"].size(); i++) {
          for (int j = 0; j < g["boundaries"][i].size(); j++) {
            int sem_index = g["semantics"]["values"][i][j];
            if (g["semantics"]["surfaces"][sem_index]["type"].get<std::string>().compare("RoofSurface") == 0) {
              std::cout << "RoofSurface: " << g["boundaries"][i][j] << std::endl;
            }
          }
        }
      }
    }
  }
}


// Returns the number of 'RooSurface' in the CityJSON model
int get_no_roof_surfaces(json &j) {
  int total = 0;
  for (auto& co : j["CityObjects"].items()) {
    for (auto& g : co.value()["geometry"]) {
      if (g["type"] == "Solid") {
        for (auto& shell : g["semantics"]["values"]) {
          for (auto& s : shell) {
            if (g["semantics"]["surfaces"][s.get<int>()]["type"].get<std::string>().compare("RoofSurface") == 0) {
              total += 1;
            }
          }
        }
      }
    }
  }
  return total;
}


// CityJSON files have their vertices compressed: https://www.cityjson.org/specs/1.1.1/#transform-object
// this function visits all the surfaces and print the (x,y,z) coordinates of each vertex encountered
void list_all_vertices(json& j) {
  for (auto& co : j["CityObjects"].items()) {
    std::cout << "= CityObject: " << co.key() << std::endl;
    for (auto& g : co.value()["geometry"]) {
      if (g["type"] == "Solid") {
        for (auto& shell : g["boundaries"]) {
          for (auto& surface : shell) {
            for (auto& ring : surface) {
              std::cout << "---" << std::endl;
              for (auto& v : ring) { 
                std::vector<int> vi = j["vertices"][v.get<int>()];
                double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
                double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
                double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
                std::cout << std::setprecision(2) << std::fixed << v << " (" << x << ", " << y << ", " << z << ")" << std::endl;                
              }
            }
          }
        }
      }
    }
  }
}

int get_no_ground_surfaces(json &j) {
    int total = 0;
    for (auto& co : j["CityObjects"].items()) {
        for (auto& g : co.value()["geometry"]) {
            if (g["type"] == "Solid") {
                for (auto& shell : g["semantics"]["values"]) {
                    for (auto& s : shell) {
                        if (g["semantics"]["surfaces"][s.get<int>()]["type"].get<std::string>().compare("GroundSurface") == 0) {
                            total += 1;
                        }
                    }
                }
            }
        }
    }
    return total;
}

