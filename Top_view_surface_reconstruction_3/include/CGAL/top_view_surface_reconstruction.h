#ifndef CGAL_TOP_VIEW_SURFACE_RECONSTRUCTION_H
#define CGAL_TOP_VIEW_SURFACE_RECONSTRUCTION_H

#include <CGAL/Top_view_surface_reconstruction_3/Surface_mesh_on_cdt.h>
#include <CGAL/Top_view_surface_reconstruction_3/Border_graph.h>

#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <queue>

namespace CGAL
{

namespace Top_view_surface_reconstruction_3
{
  
namespace internal
{


  template <typename Point>
  struct Sort_by_z
  {
    bool operator()(const Point& a, const Point& b)
    {
      return a.z() < b.z();
    }
  };

  template <typename GeomTraits, typename PointInputIterator, typename PointMap>
  void insert_filtered_points (PointInputIterator begin,
                               PointInputIterator end,
                               PointMap point_map,
                               Surface_mesh_on_cdt<GeomTraits>& output_mesh,
                               double spacing,
                               double quantile)
  {
    typedef typename GeomTraits::Point_3 Point_3;
    typedef std::map<std::pair<int, int>, std::vector<Point_3> > Grid;
    Grid grid;

    for (PointInputIterator it = begin; it != end; ++ it)
    {
      std::pair<int, int> coord = std::make_pair ((int)(get(point_map, *it).x() / spacing),
                                                  (int)(get(point_map, *it).y() / spacing));

      typename Grid::iterator found = grid.find (coord);
      if (found == grid.end())
        grid.insert (std::make_pair (coord, std::vector<Point_3> (1, get(point_map, *it))));
      else
        found->second.push_back (get (point_map, *it));
    }
    
    for (typename Grid::iterator it = grid.begin(); it != grid.end(); ++ it)
    {
      std::vector<Point_3>& cell_pts = it->second;
      std::sort (cell_pts.begin(), cell_pts.end(), Sort_by_z<Point_3>());
      output_mesh.insert (cell_pts[std::min (cell_pts.size() - 1, std::size_t(cell_pts.size() * quantile))]);
    }
  }

  template <typename GeomTraits>
  bool face_with_height_under_tolerance (Surface_mesh_on_cdt<GeomTraits>& mesh,
                                         typename Surface_mesh_on_cdt<GeomTraits>::Face_handle f,
                                         double epsilon)
  {
    for (std::size_t i = 0; i < 3; ++ i)
    {
      const typename GeomTraits::Point_3& a = mesh.point(f->vertex(i)->info()[0]);
      const typename GeomTraits::Point_3& b = mesh.point(f->vertex((i+1)%3)->info()[0]);
      if (std::fabs (a.z() - b.z()) > epsilon)
        return false;
    }
    return true;
  }

  template <typename GeomTraits>
  bool cluster_is_large_enough (std::vector<typename Surface_mesh_on_cdt<GeomTraits>::Face_handle>& faces,
                                double epsilon)
  {
    std::vector<typename GeomTraits::Triangle_2> tri;
    for (std::size_t i = 0; i < faces.size(); ++ i)
      tri.push_back (typename GeomTraits::Triangle_2 (faces[i]->vertex(0)->point(),
                                                      faces[i]->vertex(1)->point(),
                                                      faces[i]->vertex(2)->point()));
    typename GeomTraits::Line_2 line;
    typename GeomTraits::Point_2 centroid;
    CGAL::linear_least_squares_fitting_2 (tri.begin(), tri.end(),
                                          line, centroid, CGAL::Dimension_tag<2>());

    typename GeomTraits::Vector_2 ref = line.to_vector().perpendicular(CGAL::POSITIVE);
    ref = ref / std::sqrt (ref*ref);

    double min = std::numeric_limits<double>::max();
    double max = -std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < faces.size(); ++ i)
      for (std::size_t j = 0; j < 3; ++ j)
        {
          typename GeomTraits::Vector_2 v (centroid, faces[i]->vertex(j)->point());
          double coord = v * ref;
          min = (std::min)(coord, min);
          max = (std::max)(coord, max);
        }
    return ((max - min) > epsilon);
  }

  template <typename GeomTraits>
  void filter_faces (Surface_mesh_on_cdt<GeomTraits>& mesh,
                     double spacing, double meshing_factor)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    double epsilon = spacing * meshing_factor;

    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin();
         it != mesh.finite_faces_end (); ++ it)
      if (!(face_with_height_under_tolerance<GeomTraits> (mesh, it, epsilon)))
        mesh.skip (it);
    
    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end (); ++ it)
    {
      typename SMCDT::Face_circulator circ = mesh.incident_faces (it);
      typename SMCDT::Face_circulator start = circ;

      std::size_t nb_layers = 0;
      bool in = false;
      bool first_in = false, last_in = false;
      do
      {
        if (mesh.is_infinite(circ))
        {
          ++ circ;
          continue;
        }

        if (mesh.is_skiped(circ))
        {
          if (!in)
            ++ nb_layers;
          in = true;
        }
        else
        {
          in = false;
        }
        last_in = in;
        if (circ == start)
          first_in = in;
        
        ++ circ;
      }
      while (circ != start);

      if (first_in && last_in)
        nb_layers --;
      if (nb_layers > 1)
      {
        circ = start;
        do
        {
          mesh.skip(circ);
          ++ circ;
        }
        while (circ != start);
      }
    }
    
    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin();
         it != mesh.finite_faces_end (); ++ it)
    {
      if (mesh.has_mesh_face (it) || mesh.is_skiped (it))
        continue;

      std::queue<typename SMCDT::Face_handle> todo;
      std::vector<typename SMCDT::Face_handle> done;
      todo.push (it);

      while (!(todo.empty()))
      {
        typename SMCDT::Face_handle current = todo.front();
        todo.pop();

        if (mesh.has_mesh_face (current) || mesh.is_skiped (current))
          continue;
        
        current->info() = typename SMCDT::Face_index(0);

        done.push_back (current);
        for (std::size_t i = 0; i < 3; ++ i)
        {
          if (mesh.is_infinite (current->neighbor (i))
              || mesh.has_mesh_face (current->neighbor(i))
              || mesh.is_skiped (current->neighbor(i)))
            continue;

          todo.push (current->neighbor(i));
          
        }
      }

      if (cluster_is_large_enough<GeomTraits> (done, meshing_factor * epsilon))
        for (std::size_t i = 0; i < done.size(); ++ i)
          mesh.add_face (done[i]);
      else
        for (std::size_t i = 0; i < done.size(); ++ i)
          mesh.skip (done[i]);
    }

  }
    
  template <typename GeomTraits>
  void cleanup_buffer_zone (Surface_mesh_on_cdt<GeomTraits>& mesh)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    // Constrain edges
    for (typename SMCDT::Finite_edges_iterator it = mesh.finite_edges_begin ();
         it != mesh.finite_edges_end(); ++ it)
    {
      typename SMCDT::Face_handle f0 = it->first;
      typename SMCDT::Face_handle f1 = it->first->neighbor (it->second);
        
      if (!(mesh.is_infinite (f0)) && !(mesh.is_infinite (f1))
          && ((mesh.is_skiped(f0) && !(mesh.is_skiped(f1))) ||
              (!(mesh.is_skiped(f0)) && mesh.is_skiped(f1))))
        mesh.insert_constraint (it->first->vertex ((it->second + 1)%3),
                                it->first->vertex ((it->second + 2)%3));
    }

    // Remove isolated vertices in buffer zone
    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end();)
    {
      typename SMCDT::Finite_vertices_iterator current = it ++;
      
      if (mesh.are_there_incident_constraints(current))
        continue;
        
      typename SMCDT::Face_circulator circ = mesh.incident_faces (current);
      typename SMCDT::Face_circulator start = circ;
      bool rem = true;
      do
      {
        if (!mesh.is_infinite (circ) && mesh.has_mesh_face(circ))
        {
          rem = false;
          break;
        }
        ++ circ;
      }
      while (circ != start);

      if (rem)
        mesh.remove (current);
    }

    // Mark new faces as skiped
    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin ();
         it != mesh.finite_faces_end(); ++ it)
      if (!mesh.has_mesh_face(it))
        mesh.skip(it);
  }

  template <typename GeomTraits>
  typename GeomTraits::Point_2 midpoint (typename Surface_mesh_on_cdt<GeomTraits>::Face_handle f)
  {
    return typename GeomTraits::Point_2
      ((f->vertex(0)->point().x() + f->vertex(1)->point().x() + f->vertex(2)->point().x()) / 3.,
       (f->vertex(0)->point().y() + f->vertex(1)->point().y() + f->vertex(2)->point().y()) / 3.);
  }
  
  template <typename GeomTraits>
  void build_border_graph (Surface_mesh_on_cdt<GeomTraits>& mesh,
                           Border_graph<GeomTraits>& graph)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
    typedef Border_graph<GeomTraits> Graph;
    
    typedef std::map<typename SMCDT::Face_handle, typename Graph::vertex_descriptor> Map_v2v;
    
    Map_v2v map_v2v;

    for (typename SMCDT::Finite_edges_iterator it = mesh.finite_edges_begin();
         it != mesh.finite_edges_end(); ++ it)
    {
      if (mesh.is_infinite (it->first) ||
          mesh.is_infinite (it->first->neighbor(it->second)))
        continue;

      if (mesh.is_skiped(it->first) && mesh.is_skiped(it->first->neighbor(it->second)))
      {
        typename Map_v2v::iterator map0, map1;
        bool inserted = false;
        boost::tie (map0, inserted) = map_v2v.insert (std::make_pair (it->first, typename Graph::vertex_descriptor()));
        if (inserted)
        {
          typename GeomTraits::Point_2 p0 = midpoint<GeomTraits> (it->first);
          map0->second = boost::add_vertex(p0, graph);
        }

        boost::tie (map1, inserted) = map_v2v.insert (std::make_pair (it->first->neighbor(it->second),
                                                                      typename Graph::vertex_descriptor()));
        if (inserted)
        {
          typename GeomTraits::Point_2 p1 = midpoint<GeomTraits> (it->first->neighbor(it->second));
          map1->second = boost::add_vertex(p1, graph);
        }

        boost::add_edge (map0->second, map1->second, graph);
      }
    }
  }

  template <typename GeomTraits>
  void simplify_border_graph (Surface_mesh_on_cdt<GeomTraits>& mesh,
                              Border_graph<GeomTraits>& graph,
                              double epsilon)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
    typedef Border_graph<GeomTraits> Graph;

    for (std::size_t i = 0; i < graph.size(); ++ i)
    {
      std::vector<typename GeomTraits::Triangle_2> facets;
      typename SMCDT::Face_handle hint = typename SMCDT::Face_handle();
      for (std::size_t j = 0; j < graph[i].size(); ++ j)
      {
        hint = mesh.locate (graph.point(graph[i][j]));
        facets.push_back (mesh.triangle (hint));
      }

      std::vector<std::size_t> tokeep;
      tokeep.push_back (0);
      tokeep.push_back (graph[i].size() - 1);
      std::queue<std::pair<std::size_t, std::size_t> > todo;
      todo.push (std::make_pair (0, graph[i].size() - 1));

      while (!todo.empty())
      {
        std::size_t first = todo.front().first;
        std::size_t last = todo.front().second;
        todo.pop();
        
        std::size_t ind_max = std::size_t(-1);
        double max = 0.;

        // Handle closed polylines
        if (graph[i][first] == graph[i][last])
          for (std::size_t j = first+1; j < last; ++ j)
          {
            double dist = CGAL::squared_distance (graph.point(graph[i][j]),
                                                  graph.point(graph[i][first]));
            if (dist > max)
            {
              max = dist;
              ind_max = j;
            }
          }
        else
        {
          typename GeomTraits::Segment_2 s (graph.point(graph[i][first]),
                                            graph.point(graph[i][last]));
    
          for (std::size_t j = first+1; j < last; ++ j)
          {
            if (CGAL::do_intersect (s, facets[j]))
              continue;
        
            double dist = CGAL::squared_distance (s, facets[j]);
            if (dist > max)
            {
              max = dist;
              ind_max = j;
            }
          }
        }
    
        if (max > epsilon * epsilon)
        {
          tokeep.push_back (ind_max);

          if (first < ind_max)
            todo.push (std::make_pair (first, ind_max));
          if (ind_max < last)
            todo.push (std::make_pair (ind_max, last));
        }
      }

      std::sort (tokeep.begin(), tokeep.end());
    }
  }
  
} // namespace internal

} // namespace Top_view_surface_reconstruction_3
  
template <typename GeomTraits, typename PointInputIterator, typename PointMap>
void top_view_surface_reconstruction (PointInputIterator begin,
                                      PointInputIterator end,
                                      PointMap point_map,
                                      Surface_mesh_on_cdt<GeomTraits>& output_mesh,
                                      double spacing,
                                      double meshing_factor = 3.,
                                      double quantile = 0.5)
{
//  typedef typename GeomTraits::Point_3 Point_3;
//  typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
  typedef Border_graph<GeomTraits> Graph;
  
  Top_view_surface_reconstruction_3::internal::insert_filtered_points<GeomTraits>
    (begin, end, point_map, output_mesh, spacing, quantile);

  Top_view_surface_reconstruction_3::internal::filter_faces<GeomTraits>
    (output_mesh, spacing, meshing_factor);

  output_mesh.DEBUG_dump_off_0();
  
  Top_view_surface_reconstruction_3::internal::cleanup_buffer_zone<GeomTraits>
    (output_mesh);

  output_mesh.DEBUG_dump_off_2();

  Graph graph;

  Top_view_surface_reconstruction_3::internal::build_border_graph<GeomTraits>
    (output_mesh, graph);

  graph.filter_small_terminal_borders(spacing * meshing_factor);
  graph.split_into_polylines();

  Top_view_surface_reconstruction_3::internal::simplify_border_graph<GeomTraits>
    (output_mesh, graph, spacing * meshing_factor);
  
  graph.DEBUG_dump_poly();

}

} // namespace CGAL

#endif // CGAL_TOP_VIEW_SURFACE_RECONSTRUCTION_H
