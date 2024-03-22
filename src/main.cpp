#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <numeric>

//-- https://github.com/nlohmann/json
//-- used to read and write (City)JSON
#include "json.hpp" //-- it is in the /include/ folder

// CGAL libraries
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
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
typedef CGAL::Polygon_2<Kernel> Polygon_2;

using json = nlohmann::json;

double findSurfaceArea(json j, json surface, Plane3 bestplane);
bool arePointsCollinear(const std::vector<Point2>& points);
double calculateRidge(std::vector<double> areas, std::vector<double> maxZvalues);
double calculateEave(std::vector<double> areas, std::vector<double> minZvalues);
std::vector<std::vector<std::vector<int>>> extrudeGroundSurface(json& j,
                                                                std::vector<std::vector<std::vector<int>>> groundSurface,
                                                                double &newZ);
std::vector<std::vector<int>> constructWallBoundaries(
        const std::vector<std::vector<std::vector<int>>>& groundSurfaceBoundaries,
        const std::vector<std::vector<std::vector<int>>>& extrudedGroundSurfaceBoundaries);
std::vector<std::vector<std::vector<std::vector<int>>>> createSolid(
        const std::vector<std::vector<std::vector<int>>>& groundSurfaceBoundaries,
        const std::vector<std::vector<std::vector<int>>>& extrudedGroundSurfaceBoundaries,
        const std::vector<std::vector<int>>& wallBoundaries);
std::vector<std::vector<std::vector<int>>> flipOrientation(std::vector<std::vector<std::vector<int>>>& surfaces);
std::tuple<std::vector<std::pair<std::vector<std::vector<int>>, double>>,
        std::vector<double>, std::vector<double>, std::vector<double>>
findPotentialSurfaces(json& j, json& surfacesArray);


int main(int argc, const char * argv[]) {
    //-- will read the file passed as argument or twobuildings.city.json if nothing is passed
    const char* filename = (argc > 1) ? argv[1] : "../data/tudcampus.city.json";
    std::cout << "Processing: " << filename << std::endl;
    std::ifstream input(filename);
    json j;
    input >> j; //-- store the content of the file in a nlohmann::json object
    input.close();

    // iterate over all the city objects
    for (auto& co : j["CityObjects"].items()) {
        // initialize arrays for storing of surfaces
        std::vector<std::vector<std::vector<int>>> groundSurfaceBoundaries;
        std::vector<std::vector<std::vector<int>>> roofSurfaceBoundaries;
        std::vector<std::pair<std::vector<std::vector<int>>, double>> potentialSurfaces;
        std::vector<double> surfaceAreasRoofSurfaces;
        std::vector<double> maxZRoofSurfaces;
        std::vector<double> minZRoofSurfaces;

        // iterating over boundaries
        for (auto &geom: co.value()["geometry"]) {
            if (geom.contains("boundaries")) {

                if (geom["type"] == "Solid") {
                    // go one layer deeper for solids
                    // iterating over surfaces collection
                    for (auto &shell: geom["boundaries"]) {

                        auto [potentialRoofGround, surfaceAreas, maxZs, minZs] = findPotentialSurfaces(j, shell);

                        // for each collection of surfaces extract potentials, areas, min and max values.
                        potentialSurfaces.insert(potentialSurfaces.end(), potentialRoofGround.begin(), potentialRoofGround.end());
                        surfaceAreasRoofSurfaces.insert(surfaceAreasRoofSurfaces.end(), surfaceAreas.begin(), surfaceAreas.end());
                        maxZRoofSurfaces.insert(maxZRoofSurfaces.end(), maxZs.begin(), maxZs.end());
                        minZRoofSurfaces.insert(minZRoofSurfaces.end(), minZs.begin(), minZs.end());

                    }

                } else if (geom["type"] == "MultiSurface") {
                    auto [potentialRoofGround, surfaceAreas, maxZs, minZs] = findPotentialSurfaces(j, geom["boundaries"]);

                    // for each collection of surfaces extract potentials, areas, min and max values.
                    potentialSurfaces.insert(potentialSurfaces.end(), potentialRoofGround.begin(), potentialRoofGround.end());
                    surfaceAreasRoofSurfaces.insert(surfaceAreasRoofSurfaces.end(), surfaceAreas.begin(), surfaceAreas.end());
                    maxZRoofSurfaces.insert(maxZRoofSurfaces.end(), maxZs.begin(), maxZs.end());
                    minZRoofSurfaces.insert(minZRoofSurfaces.end(), minZs.begin(), minZs.end());
                }


                // potential surfaces within the lowest 10 percent are stored as ground
                if (!potentialSurfaces.empty()) {
                    std::sort(potentialSurfaces.begin(), potentialSurfaces.end(), [](const auto &a, const auto &b) {
                        return a.second < b.second;
                    });

                    // calculate 10 percent threshold
                    double minZ = potentialSurfaces[0].second;
                    double maxZ = potentialSurfaces[potentialSurfaces.size() - 1].second;
                    double threshold = minZ + (0.1 * (maxZ - minZ));

                    for (const auto &surface: potentialSurfaces) {
                        if (surface.second <= threshold) {
                            groundSurfaceBoundaries.push_back(surface.first);
                        } else if (surface.second > threshold) {
                            roofSurfaceBoundaries.push_back(surface.first);
                        }
                    }
                }
            }
        }

        // remove groundSurfaces from area, minZ, maxZ vectors
        int n = groundSurfaceBoundaries.size();

        std::vector<int> indices(minZRoofSurfaces.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::nth_element(indices.begin(), indices.begin() + n, indices.end(),
                         [&minZRoofSurfaces](int i1, int i2) { return minZRoofSurfaces[i1] < minZRoofSurfaces[i2]; });
        std::sort(indices.begin(), indices.begin() + n, std::greater<>()); // Sort the first n indices in descending order for safe removal

        // remove elements from vectors using the identified indices
        for (int i = 0; i < n; ++i) {
            int idxToRemove = indices[i];
            surfaceAreasRoofSurfaces.erase(surfaceAreasRoofSurfaces.begin() + idxToRemove);
            maxZRoofSurfaces.erase(maxZRoofSurfaces.begin() + idxToRemove);
            minZRoofSurfaces.erase(minZRoofSurfaces.begin() + idxToRemove);
        }

        // calculate ridge, eave and lod1.2 z-value
        double ridge = calculateRidge(surfaceAreasRoofSurfaces, maxZRoofSurfaces);
        double eave = calculateEave(surfaceAreasRoofSurfaces, minZRoofSurfaces);
        double lod12z = ((ridge - eave) * 0.7) + eave;

        // Construct LoD1.2
        std::vector<std::vector<std::vector<int>>> extrudedGroundSurfaceBoundaries = extrudeGroundSurface(j, groundSurfaceBoundaries, lod12z);

        // Extract wall boundaries
        std::vector<std::vector<int>> wallBoundaries = constructWallBoundaries(groundSurfaceBoundaries, extrudedGroundSurfaceBoundaries);

        // reverse roof for better view results in CityJSON Ninja
        std::vector<std::vector<std::vector<int>>> reversedRoof = flipOrientation(extrudedGroundSurfaceBoundaries);

        // create Solid
        auto solid = createSolid(groundSurfaceBoundaries, reversedRoof, wallBoundaries);

        // write lod 1.2 to file
        json SolidLOD12 = {
                {"type", "Solid"},
                {"lod", "1.2"},
                {"boundaries", solid}
        };

        co.value()["geometry"].push_back(SolidLOD12);

        // reverse ground as well for better view in CityJSON Ninja
        std::vector<std::vector<std::vector<int>>> reversedGround = flipOrientation(groundSurfaceBoundaries);

        // write lod 0.2 to file
        json MultiSurfaceLOD02 = {
                {"type", "MultiSurface"},
                {"lod", "0.2"},
                {"boundaries", reversedGround}
        };

        co.value()["geometry"].push_back(MultiSurfaceLOD02);
    }

    //-- write to disk the modified city model (out.city.json)
    std::ofstream o("out.city.json");
    o << j.dump(1) << std::endl;
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


double findSurfaceArea(json j, json surface, Plane3 bestplane) {
    // input is surface with possible extra loops for interior rings.

    // initialize triangulation
    Triangulation triangulation;

    // store the areas of all loops
    double faceArea = 0;

    for (auto &loop : surface) {
        std::vector<Point2> loopPoints2D;

        for (const auto& vertex : loop) {
            std::vector<int> vi = j["vertices"][vertex.get<int>()];
            double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
            double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
            double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
            Point3 point3d(x, y, z);
            Point2 point2d = bestplane.to_2d(point3d);
            loopPoints2D.push_back(point2d);
        }

        // return zero area if points are collinear
        if (arePointsCollinear(loopPoints2D)) {
            return 0.0;
        }

        // insert points in triangulation object
        for (const Point2& pt : loopPoints2D) {
            triangulation.insert(pt);
        }

        if (!loopPoints2D.empty()) {
            loopPoints2D.push_back(loopPoints2D[0]);
        }

        // insert constraints to triangulation
        for (int i = 0; i < loopPoints2D.size() - 1; ++i) {
            triangulation.insert_constraint(loopPoints2D[i], loopPoints2D[i + 1]);
        }

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
        // if it belongs to the interior
        if (it->info().interior) {
            Point2 v0 = it->vertex(0)->point();
            Point2 v1 = it->vertex(1)->point();
            Point2 v2 = it->vertex(2)->point();

            double areaTriangle = CGAL::area(v0, v1, v2);
            faceArea += areaTriangle;
        }
    }
    return faceArea;
}

bool arePointsCollinear(const std::vector<Point2>& points) {
    if (points.size() < 3) {
        return true;
    }

    // Check every consecutive triplet for collinearity
    for (size_t i = 2; i < points.size(); ++i) {
        if (!CGAL::collinear(points[i - 2], points[i - 1], points[i])) {
            return false;
        }
    }

    return true;
}

double calculateRidge(std::vector<double> areas, std::vector<double> maxZvalues) {
    double totalArea = 0;
    double ridge = 0;

    for (double area : areas) {
        totalArea += area;
    }

    for (int i = 0; i < maxZvalues.size(); ++i) {
        ridge += ((areas[i]/totalArea) * maxZvalues[i]);

    }

    return ridge;
}

double calculateEave(std::vector<double> areas, std::vector<double> minZvalues) {
    double totalArea = 0;
    double eave = 0;

    for (double area : areas) {
        totalArea += area;
    }

    for (int i = 0; i < minZvalues.size(); ++i) {
        eave += ((areas[i]/totalArea) * minZvalues[i]);
    }

    return eave;
}


std::vector<std::vector<std::vector<int>>> extrudeGroundSurface(json& j,
                                                                std::vector<std::vector<std::vector<int>>> groundSurface,
                                                                double &newZ) {
    std::vector<std::vector<std::vector<int>>> extrudedGroundSurfaceBoundaries;

    for (const auto& surface : groundSurface) {
        std::vector<std::vector<int>> newSurface;
        for (const auto& loop : surface) {
            std::vector<int> newLoop;
            for (const auto vertexIndex : loop) {
                auto& vi = j["vertices"][vertexIndex];
                double x = (vi[0].get<double>() * j["transform"]["scale"][0].get<double>()) +
                           j["transform"]["translate"][0].get<double>();
                double y = (vi[1].get<double>() * j["transform"]["scale"][1].get<double>()) +
                           j["transform"]["translate"][1].get<double>();
                double z = newZ;

                // inverse scale and translate
                double newX = (x - j["transform"]["translate"][0].get<double>()) / j["transform"]["scale"][0].get<double>();
                double newY = (y - j["transform"]["translate"][1].get<double>()) / j["transform"]["scale"][1].get<double>();
                double newZ = (z - j["transform"]["translate"][2].get<double>()) / j["transform"]["scale"][2].get<double>();

                int newVertexIndex = j["vertices"].size();
                j["vertices"].push_back({newX, newY, newZ});

                newLoop.push_back(newVertexIndex);
            }
            newSurface.push_back(newLoop);
        }
        extrudedGroundSurfaceBoundaries.push_back(newSurface);
    }
    return extrudedGroundSurfaceBoundaries;
}

std::vector<std::vector<int>> constructWallBoundaries(
        const std::vector<std::vector<std::vector<int>>>& groundSurfaceBoundaries,
        const std::vector<std::vector<std::vector<int>>>& extrudedGroundSurfaceBoundaries) {

    std::vector<std::vector<int>> wallBoundaries;

    for (size_t surfaceIndex = 0; surfaceIndex < groundSurfaceBoundaries.size(); ++surfaceIndex) {
        for (size_t loopIndex = 0; loopIndex < groundSurfaceBoundaries[surfaceIndex].size(); ++loopIndex) {
            const auto& groundLoop = groundSurfaceBoundaries[surfaceIndex][loopIndex];
            const auto& extrudedLoop = extrudedGroundSurfaceBoundaries[surfaceIndex][loopIndex];

            for (size_t vertexIndex = 0; vertexIndex < groundLoop.size(); ++vertexIndex) {
                size_t nextVertexIndex = (vertexIndex + 1) % groundLoop.size();

                std::vector<int> wall = {
                        groundLoop[vertexIndex],                    // Current ground vertex
                        groundLoop[nextVertexIndex],                // Next ground vertex
                        extrudedLoop[nextVertexIndex],              // Next extruded ground vertex
                        extrudedLoop[vertexIndex]                   // Current extruded ground vertex
                };
                std::reverse(wall.begin(), wall.end());
                wallBoundaries.push_back(wall);
            }
        }
    }

    return wallBoundaries;
}


std::vector<std::vector<std::vector<std::vector<int>>>> createSolid(
        const std::vector<std::vector<std::vector<int>>>& groundSurfaceBoundaries,
        const std::vector<std::vector<std::vector<int>>>& extrudedGroundSurfaceBoundaries,
        const std::vector<std::vector<int>>& wallBoundaries) {

    std::vector<std::vector<std::vector<std::vector<int>>>> solid;

    std::vector<std::vector<std::vector<int>>> shell;

    for (const auto& loop : groundSurfaceBoundaries) {
        shell.push_back(loop);
    }

    // Wall surfaces
    for (const auto& wall : wallBoundaries) {
        shell.push_back({wall}); // Wrap wall in another vector to match the nesting
    }

    // Extruded surface
    for (const auto& loop : extrudedGroundSurfaceBoundaries) {
        shell.push_back(loop);
    }

    // Add the assembled shell to the solid
    solid.push_back(shell);

    return solid;
}

std::vector<std::vector<std::vector<int>>> flipOrientation(std::vector<std::vector<std::vector<int>>>& surfaces) {
    std::vector<std::vector<std::vector<int>>> flippedSurfaces;

    for (auto &surface : surfaces) {
        std::vector<std::vector<int>> reversedSurface;
        for (auto loop : surface) {
            std::reverse(loop.begin(), loop.end());
            reversedSurface.push_back(loop);
        }
        flippedSurfaces.push_back(reversedSurface);
    }

    return flippedSurfaces;
}


std::tuple<std::vector<std::pair<std::vector<std::vector<int>>, double>>,
        std::vector<double>, std::vector<double>, std::vector<double>>
findPotentialSurfaces(json& j, json& surfacesArray) {
    std::vector<std::pair<std::vector<std::vector<int>>, double>> potentialRoofGround;
    std::vector<double> surfaceAreasRoofSurfaces;
    std::vector<double> maxZRoofSurfaces;
    std::vector<double> minZRoofSurfaces;

    for (int i = 0; i < surfacesArray.size(); ++i) {
        double overallMinZ = std::numeric_limits<double>::infinity();
        double overallMaxZ = -std::numeric_limits<double>::infinity();
        double area;
        double totalHeight = 0;
        std::vector<Point3> vertex_coord;

        // retrieve coordinates from indices to find normal later
        // only access first surface because if it has interior loops it will be the same plane
        for (const auto& vertices : surfacesArray[i][0]) {
            std::vector<int> vi = j["vertices"][vertices.get<int>()];
            double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
            double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
            double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
            vertex_coord.emplace_back(x, y, z);
            totalHeight += z;

            // update min and max z-values for all surfaces
            if (z < overallMinZ) overallMinZ = z;
            if (z > overallMaxZ) overallMaxZ = z;

        }

        // find best-fitting plane for normal vector
        Plane3 plane;
        CGAL::linear_least_squares_fitting_3(vertex_coord.begin(), vertex_coord.end(), plane, CGAL::Dimension_tag<0>());

        // check if polygon is simple by converting to 2d
        Polygon_2 polygon;
        for (Point3 &vertex : vertex_coord) {
            Point2 vertex_2d = plane.to_2d(vertex);
            polygon.push_back(vertex_2d);
        }

        // if polygon is not simple, skip it
        if (!polygon.is_simple()) {
            continue;
        }

        // extract and normalize normal vector
        Vector3 normal = plane.orthogonal_vector();
        normal = normal / std::sqrt(normal.squared_length());

        // calculating average height of plane
        double averageHeight = totalHeight / vertex_coord.size();

        // if normal vector points upwards or downwards store it as potential roof or ground surface
        if (normal.z() > 0.25 || normal.z() < -0.25) {
            std::vector<std::vector<int>> surfaces;
            for (const auto& surface : surfacesArray[i]) {
                std::vector<int> exterior_interiors;
                for (const auto& loop : surface) {
                    exterior_interiors.push_back(loop);
                }
                surfaces.push_back(exterior_interiors);
            }
            potentialRoofGround.emplace_back(surfaces, averageHeight);

            // compute area to assign as weight for eave and ridge calculation
            area = findSurfaceArea(j, surfacesArray[i], plane);

            // store areas, min and max z-values
            surfaceAreasRoofSurfaces.push_back(area);
            maxZRoofSurfaces.push_back(overallMaxZ);
            minZRoofSurfaces.push_back(overallMinZ);

        }

    }

    return std::make_tuple(potentialRoofGround, surfaceAreasRoofSurfaces, maxZRoofSurfaces, minZRoofSurfaces);
}