#ifndef CGAL_LEVEL_OF_DETAIL_TRIANGULATION_MEGAFACET_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_TRIANGULATION_MEGAFACET_CREATOR_H

// STL includes.
#include <utility>

// CGAL includes.
#include <CGAL/Level_of_detail/internal/Vegetation/Tree.h>
#include <CGAL/Level_of_detail/internal/utils.h>

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>

namespace CGAL {

namespace Level_of_detail {

template<class GeomTraits, class InputRange, class PointMap, class Triangulation>
class Triangulation_megafacet_creator {
			
public:
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Triangle_2 Triangle_2;
  typedef typename InputRange::const_iterator const_iterator;
  typedef Tree<GeomTraits, InputRange, PointMap> Tree_item;
  typedef Building<GeomTraits, typename Triangulation::Face_handle> Building_item;

  typedef typename Triangulation::Face_handle Face_handle;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Edge Edge;

  using Triangulation_vertex_handle  = typename Triangulation::Vertex_handle;
  
private:

  const std::vector<Tree_item>& m_trees;
  const std::vector<Building_item>& m_buildings;
  const Triangulation& m_triangulation;

  std::vector<std::vector<Point_3> > m_trees_outline;
  std::vector<std::vector<Point_3> > m_buildings_outline;

  std::set<Face_handle> m_done;

public:

  Triangulation_megafacet_creator(const std::vector<Tree_item>& trees,
                                  const std::vector<Building_item>& buildings,
                                  const Triangulation &triangulation)
    : m_trees (trees), m_buildings (buildings),
      m_triangulation(triangulation)
  {

    for (typename Triangulation::Finite_faces_iterator it = triangulation.finite_faces_begin();
         it != triangulation.finite_faces_end(); ++ it)
    {
      if (it->info().visibility_label() == Visibility_label::INSIDE)
      {
        int idx = it->info().group_number();
        CGAL_assertion (idx != -1);

        if (!m_done.insert (it).second)
          continue;

        construct_megafacet (m_buildings_outline, it);
      }
      else if (it->info().visibility_label() == Visibility_label::VEGETATION)
      {
        int idx = it->info().group_number();
        CGAL_assertion (idx != -1);

        if (!m_done.insert (it).second)
          continue;

        construct_megafacet (m_trees_outline, it);
      }
    }

    // test to remove
    std::ofstream file ("outlines.polylines.txt");
    file.precision(18);

    for (std::size_t i = 0; i < m_trees_outline.size(); ++ i)
    {
      file << m_trees_outline[i].size() + 1;
      for (std::size_t j = 0; j <= m_trees_outline[i].size(); ++ j)
        file << " " << m_trees_outline[i][j % m_trees_outline[i].size()];
      file << std::endl;

      if (m_trees_outline[i].size() > 4)
        output_to_off (m_trees_outline[i]);
    }
    for (std::size_t i = 0; i < m_buildings_outline.size(); ++ i)
    {
      file << m_buildings_outline[i].size() + 1;
      for (std::size_t j = 0; j <= m_buildings_outline[i].size(); ++ j)
        file << " " << m_buildings_outline[i][j % m_buildings_outline[i].size()];
      file << std::endl;
      
      if (m_buildings_outline[i].size() > 4)
        output_to_off (m_buildings_outline[i]);
    }
  }

  std::size_t number_of_buildings() const { return m_buildings_outline.size(); }
  std::size_t number_of_trees() const { return m_trees_outline.size(); }

  const std::vector<Point_3>& building(std::size_t i) const { return m_buildings_outline[i]; }
  const std::vector<Point_3>& tree(std::size_t i) const { return m_trees_outline[i]; }

private:

  void output_to_off (const std::vector<Point_3>& outline) const
  {
    static int idx = 0;

    std::string fname = "facet";
    if (idx < 100)
      fname = fname + "0";
    if (idx < 10)
      fname = fname + "0";
    fname = fname + std::to_string(idx) + ".off";

    std::ofstream f (fname);
    f.precision(18);
    f << "OFF" << std::endl << outline.size() << " 1 0" << std::endl;

    for (std::size_t i = 0; i < outline.size(); ++ i)
      f << outline[i] << std::endl;

    f << outline.size();
    for (std::size_t i = 0; i < outline.size(); ++ i)
      f << " " << i;
    f << std::endl;
    
    ++ idx;
  }

  bool construct_megafacet (std::vector<std::vector<Point_3> >& outline,
                            Face_handle face)
  {

    Visibility_label label = face->info().visibility_label();
    int group = face->info().group_number();

    int ipivot = -1;
    
    for (std::size_t i = 0; i < 3; ++ i)
    {
      if (!(face->neighbor(i)->info().visibility_label() == label
            && face->neighbor(i)->info().group_number() == group))
      {
        ipivot = i;
        break;
      }
    }
    
    if (ipivot == -1)
      return false;

    Vertex_handle vsource = face->vertex(m_triangulation.ccw(ipivot));
    Vertex_handle vtarget = face->vertex(m_triangulation.cw(ipivot));
    Vertex_handle vsource_start = vsource;
    Vertex_handle vtarget_start = vtarget;

    outline.push_back (std::vector<Point_3>());
    do
    {
      outline.back().push_back (internal::point_3<Point_3>(face, m_triangulation.ccw(ipivot)));
      m_done.insert(face);
    
      Face_handle fcandidate = face;
      int icandidate = face->index(vsource);

      while (fcandidate->neighbor(icandidate)->info().visibility_label() == label
             && fcandidate->neighbor(icandidate)->info().group_number() == group)
      {
        Vertex_handle vpivot = fcandidate->vertex(ipivot);
        Face_handle fcandidate_next = fcandidate->neighbor(icandidate);
        int icandidate_next = fcandidate_next->index(vpivot);
        int ipivot_next = fcandidate_next->index(fcandidate);

        if (m_triangulation.is_infinite(fcandidate_next))
        {
          outline.pop_back();
          return false;
        }

        fcandidate = fcandidate_next;
        icandidate = icandidate_next;
        ipivot = ipivot_next;
      }

      face = fcandidate;
      ipivot = icandidate;
      vsource = face->vertex(m_triangulation.ccw(ipivot));
      CGAL_assertion(vsource == vtarget);

      vtarget = face->vertex(m_triangulation.cw(ipivot));
    }
    while (!(vtarget == vtarget_start && vsource == vsource_start));

    return true;
  }

};

} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_TREE_FACE_TAGGER_H
