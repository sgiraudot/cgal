#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Geographic_information_system/Triangulated_irregular_network.h>
#include <CGAL/boost/graph/graph_traits_TIN.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <CGAL/Surface_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Point_set_3<Point_3> Point_set;
typedef typename Point_set::Point_map Point_map;
typedef CGAL::Surface_mesh<Point_3> Mesh;

namespace GIS = CGAL::Geographic_information_system;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef GIS::Triangulated_irregular_network<Kernel, Point_set, Point_map> TIN;

int main (int argc, char** argv)
{
  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " points.ply" << std::endl;
    return EXIT_FAILURE;
  }

  // Read points
  std::ifstream ifile (argv[1], std::ios_base::binary);
  Point_set points;
  ifile >> points;
  std::cerr << points.size() << " point(s) read" << std::endl;

  // Create TIN
  TIN tin (points, points.point_map());

  // Save as Mesh
  Mesh mesh;
  CGAL::copy_face_graph (tin, mesh);
  std::ofstream ofile ("out.off", std::ios_base::binary);
  CGAL::set_binary_mode (ofile);
  ofile << mesh;

  return EXIT_SUCCESS;
}
  
