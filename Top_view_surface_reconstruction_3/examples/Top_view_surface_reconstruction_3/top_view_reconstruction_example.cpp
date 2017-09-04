#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/top_view_surface_reconstruction.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <vector>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Surface_mesh_on_cdt<Kernel> Mesh_on_cdt;

int main (int argc, char** argv)
{
  if (argc != 2)
    return EXIT_FAILURE;
  
  std::ifstream stream(argv[1]);
  Point_set points;
  stream >> points;

  Mesh_on_cdt mesh;

  CGAL::top_view_surface_reconstruction<Kernel>
    (points.begin(), points.end(), points.point_map(), mesh, 0.30);


  return EXIT_SUCCESS;
}
