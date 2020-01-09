#ifndef CGAL_GIS_TIN_H
#define CGAL_GIS_TIN_H

#include <fstream>
#include <queue>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

namespace CGAL
{

namespace Geographic_information_system
{

template <typename Kernel,
          typename PointRange,
          typename PointMap,
          typename VertexInfo = Default,
          typename FaceInfo = Default>
class Triangulated_irregular_network
  : public Delaunay_triangulation_2<Kernel,
                                    Triangulation_data_structure_2
                                    <Triangulation_vertex_base_with_info_2
                                     <std::pair<typename Kernel::FT, VertexInfo>, Kernel>,
                                     Triangulation_face_base_with_info_2<FaceInfo, Kernel> > >
{
public:
  
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Point_3 Point_3;
  
  typedef typename PointRange::const_iterator const_iterator;
  typedef typename boost::property_traits<PointMap>::key_type Point_item;
  typedef typename boost::property_traits<PointMap>::reference reference;

  typedef Triangulation_vertex_base_with_info_2<std::pair<FT, VertexInfo>, Kernel> Vbi;
  typedef Triangulation_face_base_with_info_2<FaceInfo, Kernel> Fbi;
  typedef Triangulation_data_structure_2<Vbi, Fbi> Tds;
  typedef Delaunay_triangulation_2<Kernel,Tds> Base;

  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Face_handle Face_handle;

private:

  struct Item_to_point_and_info_pair
  {
    PointMap point_map;

    Item_to_point_and_info_pair (PointMap point_map)
      : point_map (point_map) { }

    std::pair<Point_2, std::pair<FT, VertexInfo> > operator() (const Point_item& item) const
    {
      const Point_3& p = get (point_map, item);
      return std::make_pair (Point_2 (p.x(), p.y()), std::make_pair (p.z(), VertexInfo()));
    }
  };

  const PointRange& m_points;
  PointMap m_point_map;

public:

  Triangulated_irregular_network (const PointRange& points, PointMap point_map)
    : m_points(points), m_point_map(point_map)
    , Base (boost::make_transform_iterator (points.begin(), Item_to_point_and_info_pair(point_map)),
            boost::make_transform_iterator (points.end(), Item_to_point_and_info_pair(point_map)))
  {

  }

  const VertexInfo& operator[] (Vertex_handle vh) const { return vh->info().second; }
  VertexInfo& operator[] (Vertex_handle vh) { return vh->info().second; }
  const FaceInfo& operator[] (Face_handle fh) const { return fh->info(); }
  FaceInfo& operator[] (Face_handle fh) { return fh->info(); }

  Point_3 point_3 (Vertex_handle vh) const
  {
    const Point_2& p = vh->point();
    return Point_3 (p.x(), p.y(), vh->info().first);
  }

};

} // namespace Geographic_information_system

} // namespace CGAL

#endif // CGAL_GIS_TIN_H
