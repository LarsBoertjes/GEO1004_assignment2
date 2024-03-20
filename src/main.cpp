#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <limits>

//-- https://github.com/nlohmann/json
//-- used to read and write (City)JSON
#include "json.hpp" //-- it is in the /include/ folder

// CGAL libraries
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_tag Tag;
struct FaceInfo {
    bool interior, processed;
    FaceInfo() {
        processed = false;
        interior = false;
    }
};

typedef CGAL::Triangulation_vertex_base_2<Kernel> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, Kernel, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TriangulationDataStructure, Tag> Triangulation;
typedef Kernel::Point_2 Point2;
typedef Kernel::Point_3 Point3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Plane_3 Plane3;
typedef Kernel::Vector_3 Vector3;

using json = nlohmann::json;

struct Vertex {
    double x, y, z;
};

struct Face {
    std::list<int> boundary;
    Kernel::Plane_3 best_plane;
    Triangulation triangulation;
};

int   get_no_roof_surfaces(json &j);
int   get_no_ground_surfaces(json &j);
void  visit_roofsurfaces(json &j);
void  visit_groundsurfaces(json &j);
void  list_all_vertices(json &j);
double findSurfaceArea(json j, json &surface, Plane3 bestplane);


int main(int argc, const char * argv[]) {
    //-- will read the file passed as argument or twobuildings.city.json if nothing is passed
    const char* filename = (argc > 1) ? argv[1] : "../data/specialcase_2.city.json";
    std::cout << "Processing: " << filename << std::endl;
    std::ifstream input(filename);
    json j;
    input >> j; //-- store the content of the file in a nlohmann::json object
    input.close();

    // iterate over all the cityobjects
    for (auto& co : j["CityObjects"].items()) {
        std::vector<std::vector<std::vector<int>>> groundSurfaceBoundaries;
        std::vector<std::vector<std::vector<int>>> roofSurfaceBoundaries;
        std::vector<std::pair<std::vector<std::vector<int>>, double>> potentialSurfaces;
        std::vector<double> surfaceAreasRoofSurfaces;
        std::vector<double> maxZRoofSurfaces;
        std::vector<double> minZRoofSurfaces;

        // iterating over list of surfaces
        for (auto &geom: co.value()["geometry"]) {
            if (geom.contains("boundaries")) {

                double overallMinZ = std::numeric_limits<double>::infinity();
                double overallMaxZ = -std::numeric_limits<double>::infinity();
                double area = 0;

                // iterating over surfaces collection
                for (auto &shell: geom["boundaries"]) {

                    // iterate over surfaces
                    // find out if it is potentially roof or ground
                    for (int i = 0; i < shell.size(); ++i) {
                        std::vector<Point3> vertex_coord;
                        double totalHeight = 0;

                        for (auto vertex : shell[i][0]) {
                            std::vector<int> vi = j["vertices"][vertex.get<int>()];
                            double x = (vi[0] * j["transform"]["scale"][0].get<double>()) +
                                       j["transform"]["translate"][0].get<double>();
                            double y = (vi[1] * j["transform"]["scale"][1].get<double>()) +
                                       j["transform"]["translate"][1].get<double>();
                            double z = (vi[2] * j["transform"]["scale"][2].get<double>()) +
                                       j["transform"]["translate"][2].get<double>();
                            vertex_coord.emplace_back(x, y, z);
                            totalHeight += z;

                            if (z < overallMinZ) overallMinZ = z;
                            if (z > overallMaxZ) overallMaxZ = z;
                        }

                        Plane3 plane;
                        CGAL::linear_least_squares_fitting_3(vertex_coord.begin(), vertex_coord.end(), plane,
                                                             CGAL::Dimension_tag<0>());

                        Vector3 normal = plane.orthogonal_vector();
                        normal = normal / std::sqrt(normal.squared_length());
                        double averageHeight = totalHeight / vertex_coord.size();

                        // if normal vector points upwards or downwards it is a potential roof/ground surface
                        if (normal.z() > 0.25 || normal.z() < -0.25) {
                            std::vector<std::vector<int>> surfaces;
                            for (auto surface : shell[i]) {
                                std::vector<int> exterior_interiors;
                                for (auto intext : surface) {
                                    exterior_interiors.push_back(intext);
                                }
                                surfaces.push_back(exterior_interiors);
                            }
                            potentialSurfaces.push_back(std::make_pair(surfaces, averageHeight));

                            // find area to assign as weight for eave/ridge calculation
                            area = findSurfaceArea(j, shell[i], plane);
                            std::cout << "Area: " << area << std::endl;
                            std::cout << "averageHeight: " << averageHeight << std::endl;
                        }
                    }
                }

                // potential surface with lowest elevation is stored as ground, the others as roof
                if (!potentialSurfaces.empty()) {
                    std::sort(potentialSurfaces.begin(), potentialSurfaces.end(), [](const auto &a, const auto &b) {
                        return a.second < b.second;
                    });


                    double minZ = potentialSurfaces[0].second;
                    double maxZ = potentialSurfaces[potentialSurfaces.size() - 1].second;
                    double threshold = minZ + (0.1 * (maxZ - minZ));

                    for (const auto &surface: potentialSurfaces) {
                        if (surface.second <= threshold) {
                            groundSurfaceBoundaries.push_back(surface.first);
                        } else if (surface.second > threshold) {
                            roofSurfaceBoundaries.push_back(surface.first);
                            surfaceAreasRoofSurfaces.push_back(area);
                            maxZRoofSurfaces.push_back(overallMaxZ);
                            minZRoofSurfaces.push_back(overallMinZ);
                        }
                    }
                }
            }
        }


        json MultiSurfaceLOD02 = {
                {"type", "MultiSurface"},
                {"lod", "0.2"},
                {"boundaries", groundSurfaceBoundaries}
        };

        co.value()["geometry"].push_back(MultiSurfaceLOD02);

        std::cout << surfaceAreasRoofSurfaces.size() << std::endl;

        for (int i = 0; i < surfaceAreasRoofSurfaces.size(); ++i) {
            std::cout << surfaceAreasRoofSurfaces[i] << " " << maxZRoofSurfaces[i] << " " << minZRoofSurfaces[i] << std::endl;
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


void visit_groundsurfaces(json &j) {
    for (auto& co : j["CityObjects"].items()) {
        for (auto& g : co.value()["geometry"]) {
            if (g["type"] == "Solid") {
                for (int i = 0; i < g["boundaries"].size(); i++) {
                    for (int j = 0; j < g["boundaries"][i].size(); j++) {
                        int sem_index = g["semantics"]["values"][i][j];
                        if (g["semantics"]["surfaces"][sem_index]["type"].get<std::string>().compare("GroundSurface") == 0) {
                            std::cout << "GroundSurface: " << g["boundaries"][i][j] << std::endl;
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

double findSurfaceArea(json j, json &surface, Plane3 bestplane) {
    std::cout << surface << std::endl;
    // Initialize triangulation
    Triangulation triangulation;
    double surfaceArea = 0;

    std::vector<double> areas;

    for (auto &face : surface) {
        double areaFace = 0;
        std::vector<Point2> facePoints2D;

        for (auto vertex : face) {
            std::vector<int> vi = j["vertices"][vertex.get<int>()];
            double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
            double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
            double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
            Point3 point3d(x, y, z);
            Point2 point2d = bestplane.to_2d(point3d);
            facePoints2D.push_back(point2d);
        }

        for (const Point2& pt : facePoints2D) {
            triangulation.insert(pt);
        }

        facePoints2D.push_back(facePoints2D[0]);

        for (int i = 0; i < facePoints2D.size() - 1; ++i) {
            triangulation.insert_constraint(facePoints2D[i], facePoints2D[i + 1]);
        }

        // Label triangulation with interior exterior
        std::list<Triangulation::Face_handle> to_check;
        triangulation.infinite_face()->info().processed = true;
        CGAL_assertion(triangulation.infinite_face()->info().processed == true);
        CGAL_assertion(triangulation.infinite_face()->info().interior == false);
        to_check.push_back(triangulation.infinite_face());
        while (!to_check.empty()) {
            CGAL_assertion(to_check.front()->info().processed == true);
            for (int neighbour = 0; neighbour < 3; ++neighbour) {
                if (to_check.front()->neighbor(neighbour)->info().processed) {

                } else {
                    to_check.front()->neighbor(neighbour)->info().processed = true;
                    CGAL_assertion(to_check.front()->neighbor(neighbour)->info().processed == true);
                    if (triangulation.is_constrained(Triangulation::Edge(to_check.front(), neighbour))) {
                        to_check.front()->neighbor(neighbour)->info().interior = !to_check.front()->info().interior;
                        to_check.push_back(to_check.front()->neighbor(neighbour));
                    } else {
                        to_check.front()->neighbor(neighbour)->info().interior = to_check.front()->info().interior;
                        to_check.push_back(to_check.front()->neighbor(neighbour));
                    }
                }
            } to_check.pop_front();
        }

        for (auto it = triangulation.finite_faces_begin(); it != triangulation.finite_faces_end(); ++it) {
            if (it->info().interior) {
                Point2 v0 = it->vertex(0)->point();
                Point2 v1 = it->vertex(1)->point();
                Point2 v2 = it->vertex(2)->point();

                double area = CGAL::area(v0, v1, v2);
                areaFace += area;
            }

        }
        areas.push_back(areaFace);
    }

    if (areas.size() == 1) {
        return areas[0];
    } else if (areas.size() > 1) {
        double sumOfRest = std::accumulate(areas.begin() + 1, areas.end(), 0.0);
        return areas[0] - sumOfRest;
    }
}