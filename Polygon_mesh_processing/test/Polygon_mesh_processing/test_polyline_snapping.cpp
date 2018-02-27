#include <cstdlib>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/polyline_snapping.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::FT FT;
typedef CGAL::cpp11::tuple<std::size_t, std::size_t, std::size_t, std::size_t, FT, FT> Snapping_point;

using CGAL::cpp11::get;

int main (int argc, char** argv)
{
  std::vector<std::vector<Point_3> > polylines;
  
  if (argc == 1)
  {
    std::cerr << "Usage: " << argv[0] << " input1.polylines.txt input2.polylines.txt (...)" << std::endl;
    return EXIT_SUCCESS;
  }

  for (int i = 1; i < argc; ++ i)
  {
    std::ifstream stream(argv[std::size_t(i)]);
    if (!stream)
    {
      std::cerr << "Error: can't read " << argv[std::size_t(i)] << std::endl;
    }

    std::string line;
    std::istringstream iss;
    while(getline(stream,line))
    {
      iss.clear();
      iss.str(line);
      int size = 0;
      if (iss >> size && (size > 0))
      {
        polylines.push_back (std::vector<Point_3>());
        polylines.back().reserve (size);
        for (int n = 0; n < size; ++ n)
        {
          double x, y, z;
          if (iss >> CGAL::iformat(x) >> CGAL::iformat(y) >> CGAL::iformat(z))
            polylines.back().push_back (Point_3 (x, y, z));
        }
      }
    }
  }

  std::cerr << polylines.size() << " polyline(s) read" << std::endl;

  std::vector<Snapping_point> snapping_points;
  CGAL::Polygon_mesh_processing::polyline_snapping (polylines, 0.001,
                                                    std::back_inserter (snapping_points),
                                                    Kernel());
  std::cerr << snapping_points.size() << " snapping point(s) found" << std::endl;

  std::ofstream output ("snapping_points.xyz");
  output.precision(18);
  
  for (std::size_t i = 0; i < snapping_points.size(); ++ i)
  {
    const Snapping_point& s = snapping_points[i];
    output << CGAL::barycenter (polylines[get<0>(s)][get<1>(s)],     (1. - get<4>(s)),
                                polylines[get<0>(s)][get<1>(s) + 1],       get<4>(s)) << std::endl
           << CGAL::barycenter (polylines[get<2>(s)][get<3>(s)],     (1. - get<5>(s)),
                                polylines[get<2>(s)][get<3>(s) + 1],       get<5>(s)) << std::endl;
  }
  return EXIT_SUCCESS;
}
