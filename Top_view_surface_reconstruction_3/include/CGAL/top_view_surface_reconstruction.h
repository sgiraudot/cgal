#ifndef CGAL_TOP_VIEW_SURFACE_RECONSTRUCTION_H
#define CGAL_TOP_VIEW_SURFACE_RECONSTRUCTION_H

//#define TOP_VIEW_FIX_DUPLICATE_VERTICES

#define TOP_VIEW_DEBUG
#define TOP_VIEW_LOG
//#define TOP_VIEW_CHECK_STRUCTURE

#ifdef TOP_VIEW_DEBUG
#define TOP_VIEW_SILENT false
#else
#define TOP_VIEW_SILENT true
#endif

#define TOP_VIEW_CERR \
  if(TOP_VIEW_SILENT) {} else std::cerr

#include <CGAL/Top_view_surface_reconstruction_3/Surface_mesh_on_cdt.h>
#include <CGAL/Top_view_surface_reconstruction_3/Border_graph.h>
#include <CGAL/Top_view_surface_reconstruction_3/Line_detector.h>

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
  bool face_with_height_over_tolerance (Surface_mesh_on_cdt<GeomTraits>& mesh,
                                        typename Surface_mesh_on_cdt<GeomTraits>::Face_handle f,
                                        double epsilon)
  {
    for (std::size_t i = 0; i < 3; ++ i)
    {
      const typename GeomTraits::Point_3& a = mesh.point(f->vertex(i));
      const typename GeomTraits::Point_3& b = mesh.point(f->vertex((i+1)%3));
      if (std::fabs (a.z() - b.z()) > epsilon)
        return true;

    }
    return false;
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
  void filter_closed_polygons_with_no_3D_point_inside (Surface_mesh_on_cdt<GeomTraits>& mesh,
                                                       bool test_full_face)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    std::vector<std::vector<typename SMCDT::Face_handle> > polygons;

    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin();
         it != mesh.finite_faces_end (); ++ it)
    {
      if (!mesh.is_default (it))
        continue;

      if (mesh.has_at_least_one_mesh_vertex (it))
        continue;

      polygons.push_back (std::vector<typename SMCDT::Face_handle>());

      bool okay = false;
      
      std::queue<typename SMCDT::Face_handle> todo;
      todo.push(it);
      while (!todo.empty())
      {
        typename SMCDT::Face_handle current = todo.front();
        todo.pop();

        if (!mesh.is_default(current))
          continue;

        if ((test_full_face && mesh.has_three_mesh_vertices (current))
            ||(!test_full_face && mesh.has_at_least_one_mesh_vertex (current)))
          okay = true;

        current->info().index = typename SMCDT::Face_index(0); // Make non-default
        polygons.back().push_back (current);

        for (std::size_t i = 0; i < 3; ++ i)
        {
          if (mesh.is_constrained (std::make_pair (current, i)))
            continue;
          todo.push(current->neighbor(i));
        }
      }

      if (okay)
        polygons.pop_back();
    }

    std::ofstream file ("removed.xyz");
    file.precision(18);
    
    for (std::size_t i = 0; i < polygons.size(); ++ i)
    {
      typename SMCDT::Face_handle f = polygons[i][0];
      int current = -1;
      int pivot = -1;

      for (std::size_t k = 0; k < polygons[i].size(); ++ k)
      {
        f = polygons[i][k];
        for (std::size_t j = 0; j < 3; ++ j)
          if (mesh.is_constrained (std::make_pair (f, j)))
          {
            current = (j + 1)%3;
            pivot = (j + 2)%3;
            break;
          }
        if (current != -1)
          break;
      }
      CGAL_assertion (current != -1);
      
      std::vector<std::pair<typename SMCDT::Vertex_handle, std::size_t> > border_vertices;

      while (true)
      {
        int e = 3 - (current + pivot);
        if (mesh.is_constrained(std::make_pair(f, e))) // Current edge is constrained
        {
          if (!border_vertices.empty() &&
              f->vertex(current) == border_vertices.front().first)
            break;

          std::vector<typename SMCDT::Edge> edges;
          mesh.incident_constraints (f->vertex(current), std::back_inserter(edges));
          border_vertices.push_back (std::make_pair (f->vertex(current), edges.size()));

          if (mesh.is_constrained(std::make_pair(f,pivot))) // Next edge is contrained
          {
            pivot = current;
            current = e;
          }
          else // Next edge is not constrained
          {
            typename SMCDT::Face_handle prev = f;
            f = f->neighbor(pivot);
            pivot = f->index(border_vertices.back().first);
            current = f->index(prev);
          }
        }
        else // Current edge is not constrained
        {
          typename SMCDT::Face_handle prev = f;
          f = f->neighbor(e);
          pivot = f->index(border_vertices.back().first);
          current = f->index(prev);
        }
      }

      std::size_t first = std::size_t(-1);
      std::size_t latest = 0;
      std::size_t begin = 0;
      std::size_t end = 0;
      std::size_t longest = 0;
      std::size_t j = 0;

      while (true)
      {
        if (border_vertices[j].second > 2)
        {
          if (first == std::size_t(-1))
            first = j;
          else 
          {
            if (latest < j && j - latest > longest)
            {
              begin = latest;
              end = j;
              longest = j - latest;
            }
            if (latest > j && (j + border_vertices.size()) - latest > longest) // Loop case
            {
              begin = latest;
              end = j;
              longest = (j + border_vertices.size()) - latest;
            }
          }
          latest = j;
        }

        ++ j;
        if (j == border_vertices.size())
        {
          if (first == std::size_t(-1))
            break; // Closed loop
          else
            j = 0;
        }

        if (j == first)
          break;
      }

      if (first == std::size_t(-1)) // Remove all
      {
        begin = 0;
        end = 0;
      }

      j = begin;

      typename SMCDT::Vertex_handle prev = border_vertices[j].first;

      do
      {
        ++ j;
        if (j == border_vertices.size())
          j = 0;
        typename SMCDT::Vertex_handle curr = border_vertices[j].first;

        mesh.remove_constraint (prev, curr);
        
      }
      while (j != end);
      
      j = begin;
      do
      {
        ++ j;
        if (j == border_vertices.size())
          j = 0;

        if (begin == end ||
            (j != begin && j != end))
        {
          if (mesh.are_there_incident_constraints (border_vertices[j].first))
            mesh.remove_incident_constraints (border_vertices[j].first);
          file << border_vertices[j].first->point() << " 0" << std::endl;
          mesh.remove (border_vertices[j].first);
        }
        
      }
      while (j != end);
    }

    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin();
         it != mesh.finite_faces_end (); ++ it)
      mesh.make_default (it);
  }
  
  template <typename GeomTraits>
  void ignore_faces_close_to_infinite_vertex (Surface_mesh_on_cdt<GeomTraits>& mesh,
                                              double epsilon)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
    
    typename SMCDT::Vertex_handle vinf = mesh.infinite_vertex();
    typename SMCDT::Face_circulator circ = mesh.incident_faces(vinf);
    typename SMCDT::Face_circulator start = circ;
    do
    {
      mesh.make_ignored (circ);
        
      std::size_t indinf = circ->index(vinf);
      typename GeomTraits::Segment_2 seg (circ->vertex((indinf + 1)%3)->point(),
                                          circ->vertex((indinf + 2)%3)->point());

      std::queue<typename SMCDT::Edge> todo;
      todo.push(std::make_pair (circ, indinf));
      std::set<typename SMCDT::Face_handle> done;
      while (!(todo.empty()))
      {
        typename SMCDT::Edge edge = todo.front();
        todo.pop();

        typename SMCDT::Face_handle f = edge.first->neighbor(edge.second);
        if (mesh.is_infinite(f))
          continue;
        if (!(done.insert(f).second))
          continue;

        bool skip = false;
        for (std::size_t k = 0; k < 3; ++ k)
          if (CGAL::squared_distance (f->vertex(k)->point(), seg) < epsilon * epsilon
              || CGAL::squared_distance (f->vertex((k+1)%3)->point(),
                                         f->vertex((k+2)%3)->point()) > 100. * epsilon * epsilon)
          {
            skip = true;
            break;
          }

        if (skip)
        {
          typename SMCDT::Vertex_handle next = f->vertex(f->index(edge.first));
          std::size_t ind = f->index(next);
          mesh.make_ignored (f);
          todo.push (std::make_pair (f, (ind+1)%3));
          todo.push (std::make_pair (f, (ind+2)%3));
        }
      }
            
      ++ circ;
    }
    while (circ != start);

    std::ofstream file ("ignored.xyz");
    file.precision(18);
    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end(); ++ it)
    {
      if (!mesh.has_unique_mesh_vertex (it))
        continue;
      
      bool ignored = true;
      typename SMCDT::Face_circulator start = mesh.incident_faces(it);
      typename SMCDT::Face_circulator circ = start;
      do
      {
        if (!mesh.is_ignored(circ))
        {
          ignored = false;
          break;
        }
        ++ circ;
      }
      while (circ != start);

      if (ignored)
      {
        file << it->point() << " 0" << std::endl;
        mesh.remove_mesh_vertex (it);
      }
    }
    
    // for (typename SMCDT::Finite_edges_iterator it = mesh.finite_edges_begin();
    //      it != mesh.finite_edges_end(); ++ it)
    //   if (mesh.is_ignored(it->first) != mesh.is_ignored(it->first->neighbor(it->second)))
    //     mesh.insert_constraint (it->first->vertex ((it->second + 1) % 3),
    //                             it->first->vertex ((it->second + 2) % 3));
  }

  
  template <typename GeomTraits>
  void filter_faces (Surface_mesh_on_cdt<GeomTraits>& mesh,
                     double epsilon)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    ignore_faces_close_to_infinite_vertex<GeomTraits> (mesh, epsilon);

    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin();
         it != mesh.finite_faces_end (); ++ it)
      if (!mesh.is_ignored(it) && face_with_height_over_tolerance<GeomTraits> (mesh, it, epsilon))
        mesh.make_wall_buffer (it);
    
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
        if (mesh.is_ignored(circ))
        {
          ++ circ;
          continue;
        }

        if (mesh.is_buffer(circ))
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
          if (!mesh.is_ignored(circ))
              mesh.make_wall_buffer(circ);
          ++ circ;
        }
        while (circ != start);
      }
    }
    
    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin();
         it != mesh.finite_faces_end (); ++ it)
    {
      if (!mesh.is_default(it))
        continue;

      std::queue<typename SMCDT::Face_handle> todo;
      std::vector<typename SMCDT::Face_handle> done;
      todo.push (it);

      while (!(todo.empty()))
      {
        typename SMCDT::Face_handle current = todo.front();
        todo.pop();

        if (!mesh.is_default(current))
          continue;
        
        current->info().index = typename SMCDT::Face_index(0);

        done.push_back (current);
        for (std::size_t i = 0; i < 3; ++ i)
        {
          if (mesh.is_infinite (current->neighbor (i))
              || !mesh.is_default (current->neighbor(i)))
            continue;

          todo.push (current->neighbor(i));
          
        }
      }

      if (cluster_is_large_enough<GeomTraits> (done, 3. * epsilon))
        for (std::size_t i = 0; i < done.size(); ++ i)
          mesh.make_pending (done[i]);
      else
        for (std::size_t i = 0; i < done.size(); ++ i)
          mesh.make_wall_buffer (done[i]);
    }

  }

  template <typename GeomTraits>
  void filter_vertices (Surface_mesh_on_cdt<GeomTraits>& mesh,
                        double epsilon)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    std::ofstream file ("removed.xyz");
    file.precision(18);
    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end (); ++ it)
    {
      if (!mesh.has_mesh_vertex (it))
        continue;

      std::queue<typename SMCDT::Vertex_handle> todo;
      todo.push(it);

      std::set<typename SMCDT::Vertex_handle> removed;
      while (!todo.empty())
      {
        typename SMCDT::Vertex_handle current = todo.front();
        todo.pop();

        if (removed.find(current) != removed.end())
          continue;
        
        const typename GeomTraits::Point_3& a = mesh.point(current);

        typename SMCDT::Vertex_circulator circ = mesh.incident_vertices (current);
        typename SMCDT::Vertex_circulator start = circ;

        std::vector<typename SMCDT::Vertex_handle> neighbors;
        
        bool incident_to_border = false;
        bool too_long_edge = false;
        do
        {
          if (!mesh.has_mesh_vertex (circ))
            incident_to_border = true;
          else
          {
            neighbors.push_back (circ);
            
            const typename GeomTraits::Point_3& b = mesh.point(circ);

            if (CGAL::abs (a.z() - b.z()) > epsilon)
              too_long_edge = true;
          }
          ++ circ;
        }
        while (circ != start);

        if (incident_to_border && too_long_edge)
        {
          typename SMCDT::Vertex_handle vh = it;
          if (vh == current)
            ++ it;
          file << a << std::endl;
          removed.insert (current);
          mesh.remove (current);
          for (std::size_t i = 0; i < neighbors.size(); ++ i)
            todo.push(neighbors[i]);
        }
      }
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
        
      if (!(mesh.is_ignored (f0)) && !(mesh.is_ignored (f1))
          && (mesh.is_buffer(f0) != mesh.is_buffer(f1)))
        mesh.insert_constraint (it->first->vertex ((it->second + 1)%3),
                                it->first->vertex ((it->second + 2)%3));
    }


    // Remove isolated vertices in rest buffer
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
        if (mesh.is_pending(circ) || mesh.is_ignored(circ)
            || mesh.is_ridge_buffer(circ))
        {
          rem = false;
          break;
        }

        ++ circ;
      }
      while (circ != start);

      if (!rem)
        continue;

      typename GeomTraits::Point_2 point = current->point();
      mesh.remove (current);

      std::vector<typename SMCDT::Face_handle> conflict;
      mesh.get_conflicts (point, std::back_inserter (conflict));
      for (std::size_t i = 0; i < conflict.size(); ++ i)
        mesh.make_wall_buffer(conflict[i]);
    }

    // Remove isolated vertices in ridge buffer
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
        if (mesh.is_pending(circ) || mesh.is_ignored(circ))
//            || mesh.is_wall_buffer(circ))
        {
          rem = false;
          break;
        }

        ++ circ;
      }
      while (circ != start);

      if (!rem)
        continue;

      typename GeomTraits::Point_2 point = current->point();
      mesh.remove (current);

      std::vector<typename SMCDT::Face_handle> conflict;
      mesh.get_conflicts (point, std::back_inserter (conflict));
      for (std::size_t i = 0; i < conflict.size(); ++ i)
        mesh.make_ridge_buffer(conflict[i]);
    }

  }

  template <typename GeomTraits>
  void build_border_graph (Surface_mesh_on_cdt<GeomTraits>& mesh,
                           Border_graph<GeomTraits>& graph,
                           double epsilon)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
    typedef Border_graph<GeomTraits> Graph;
    
    typedef std::map<typename SMCDT::Face_handle, typename Graph::vertex_descriptor> Map_f2v;
    
    Map_f2v map_f2v;

    for (typename SMCDT::Finite_edges_iterator it = mesh.finite_edges_begin();
         it != mesh.finite_edges_end(); ++ it)
    {
      if (mesh.is_ignored (it->first) ||
          mesh.is_ignored (it->first->neighbor(it->second)))
        continue;

      if (mesh.is_buffer(it->first) && mesh.is_buffer(it->first->neighbor(it->second)))
      {
        typename Map_f2v::iterator map0, map1;
        bool inserted = false;
        boost::tie (map0, inserted) = map_f2v.insert (std::make_pair (it->first, typename Graph::vertex_descriptor()));
        if (inserted)
        {
          typename GeomTraits::Point_2 p0 = mesh.midpoint (it->first);
          map0->second = boost::add_vertex(p0, graph);
        }

        boost::tie (map1, inserted) = map_f2v.insert (std::make_pair (it->first->neighbor(it->second),
                                                                      typename Graph::vertex_descriptor()));
        if (inserted)
        {
          typename GeomTraits::Point_2 p1 = mesh.midpoint (it->first->neighbor(it->second));
          map1->second = boost::add_vertex(p1, graph);
        }

        boost::add_edge (map0->second, map1->second, graph);
      }
    }

    std::vector<typename GeomTraits::Point_2> removed_points;

    graph.filter_small_terminal_borders(epsilon, std::back_inserter (removed_points));

    typename SMCDT::Face_handle fhint;
    
    for (std::size_t i = 0; i < removed_points.size(); ++ i)
    {
      fhint = mesh.locate(removed_points[i]);
      mesh.make_ignored_buffer(fhint);
    }
    
    graph.split_into_polylines();
  }

  template <typename GeomTraits>
  void simplify_border_graph (Surface_mesh_on_cdt<GeomTraits>& mesh,
                              Border_graph<GeomTraits>& graph,
                              double epsilon)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
    typedef Border_graph<GeomTraits> Graph;

    std::vector<std::vector<typename GeomTraits::Triangle_2> > all_facets (graph.size());

    std::vector<std::vector<typename GeomTraits::Line_2> > all_supports (graph.size());

    TOP_VIEW_CERR << "  Simplification" << std::endl;
    for (std::size_t i = 0; i < graph.size(); ++ i)
    {
      std::vector<typename GeomTraits::Triangle_2>& facets
        = all_facets[i];
      
      std::vector<typename GeomTraits::Line_2>& support
        = all_supports[i];
      
      typename SMCDT::Face_handle hint = typename SMCDT::Face_handle();
      for (std::size_t j = 0; j < graph[i].size(); ++ j)
      {
        hint = mesh.locate (graph.point(graph[i][j]), hint);
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
      for (std::size_t j = 0; j < graph[i].size() - 1; ++ j)
        remove_edge (graph[i][j], graph[i][j+1], graph);

      std::vector<typename Graph::vertex_descriptor> new_polyline;
      for (std::size_t j = 0; j < tokeep.size(); ++ j)
      {
        new_polyline.push_back (graph[i][tokeep[j]]);
        if (j != tokeep.size() - 1)
        {
          add_edge (graph[i][tokeep[j]], graph[i][tokeep[j+1]], graph);

          // Compute support
          typename GeomTraits::Line_2 line;
          typename GeomTraits::Point_2 centroid;
          CGAL::linear_least_squares_fitting_2 (facets.begin() + tokeep[j],
                                                facets.begin() + tokeep[j+1],
                                                line, centroid,
                                                CGAL::Dimension_tag<2>());
          support.push_back (line);
        }
      }

      for (std::size_t j = 0; j < graph[i].size(); ++ j)
        if (degree (graph[i][j], graph) == 0)
          remove_vertex (graph[i][j], graph);

      graph[i].swap (new_polyline);
    }

  }

  template <typename GeomTraits>
  void create_mesh_with_borders (Surface_mesh_on_cdt<GeomTraits>& input,
                                 Line_detector<GeomTraits>& lines,
                                 Surface_mesh_on_cdt<GeomTraits>& output,
                                 double epsilon,
                                 double faces_epsilon,
                                 bool planes_computed)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
    typedef Line_detector<GeomTraits> Detector;

    typename SMCDT::CDT map;

    // Insert borders in CDT and map
    for (std::size_t i = 0; i < lines.size(); ++ i)
    {
      typename Detector::iterator it = lines.begin(i);
      typename SMCDT::Vertex_handle previous_cdt = output.insert (*it);
      typename SMCDT::Vertex_handle previous_map = map.insert (*it);

      ++ it;
      for (; it != lines.end(i); ++ it)
      {
        const typename GeomTraits::Point_2& next = *it;
        typename GeomTraits::Vector_2 vec (previous_cdt->point(), next);
        std::size_t nb_pts = std::max(std::size_t(1), std::size_t(std::sqrt(vec*vec) / epsilon));

        typename SMCDT::Vertex_handle current_cdt = previous_cdt;

        for (std::size_t k = 1; k <= nb_pts; ++ k)
        {
          typename GeomTraits::Point_2 point
            = previous_cdt->point() + (k / double(nb_pts)) * vec;

          if (k == nb_pts)
            point = next;
          
          typename SMCDT::Vertex_handle next_cdt = output.insert (point);
          output.insert_constraint (current_cdt, next_cdt);
          current_cdt = next_cdt;
        }
        previous_cdt = current_cdt;

        typename SMCDT::Vertex_handle current_map = map.insert (next);
        map.insert_constraint (previous_map, current_map);
        previous_map = current_map;
      }
    }

    // Insert vertices not too close to borders
    typename SMCDT::Face_handle located = typename SMCDT::Face_handle();

    for (typename SMCDT::Finite_vertices_iterator it = input.finite_vertices_begin();
         it != input.finite_vertices_end(); ++ it)
    {
      located = map.locate (it->point(), located);

      bool far_enough = true;

      std::queue<typename SMCDT::Edge> todo;
      for (std::size_t i = 0; i < 3; ++ i)
        todo.push(std::make_pair (located, i));

      while (!todo.empty())
      {
        typename SMCDT::Edge e = todo.front();
        todo.pop();
        
        typename GeomTraits::Segment_2 s (e.first->vertex ((e.second+1)%3)->point(),
                                          e.first->vertex ((e.second+2)%3)->point());
        if (!map.is_constrained (e))
        {
          if (CGAL::squared_distance (it->point(), s) < epsilon * epsilon)
          {
            typename SMCDT::Face_handle next_face = e.first->neighbor(e.second);
            int next_idx = next_face->index(e.first);
            todo.push(std::make_pair (next_face, (next_idx + 1) % 3));
            todo.push(std::make_pair (next_face, (next_idx + 2) % 3));
          }
        }
        else
        {
          if (CGAL::squared_distance (it->point(), s) < epsilon * epsilon)
          {
            far_enough = false;
            break;
          }
        }
      }

      if (far_enough)
        output.insert (input.point(it));
    }

    filter_vertices<GeomTraits> (output, faces_epsilon);
    
    filter_closed_polygons_with_no_3D_point_inside<GeomTraits> (output, planes_computed);

    ignore_faces_close_to_infinite_vertex<GeomTraits> (output, epsilon);
    
    // Create existing faces
    for (typename SMCDT::Finite_faces_iterator it = output.finite_faces_begin();
         it != output.finite_faces_end (); ++ it)
      if (!output.is_ignored(it))
      {
        bool to_insert = true;
        for (std::size_t i = 0; i < 3; ++ i)
          if (!output.has_mesh_vertex (it->vertex(i)))
          {
            to_insert = false;
            break;
          }
        if (to_insert)
          output.add_face(it);
      }
    
  }

  template <typename GeomTraits>
  void create_3d_regions (Surface_mesh_on_cdt<GeomTraits>& input,
                          Line_detector<GeomTraits>& lines,
                          Surface_mesh_on_cdt<GeomTraits>& output,
                          double epsilon)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
    typedef Line_detector<GeomTraits> Detector;

    typename SMCDT::CDT map;

    // Insert borders in CDT
    for (std::size_t i = 0; i < lines.size(); ++ i)
    {
      if (lines.is_ridge(i))
        continue;
      
      typename Detector::iterator it = lines.begin(i);
      typename SMCDT::Vertex_handle previous_cdt = output.insert (*it);

      ++ it;
      for (; it != lines.end(i); ++ it)
      {
        typename SMCDT::Vertex_handle current_cdt = output.insert (*it);
        output.insert_constraint (previous_cdt, current_cdt);
        previous_cdt = current_cdt;
      }
    }
    
    output.DEBUG_dump_poly();

    // Create 3D vertices    
    for (typename SMCDT::Finite_vertices_iterator it = output.finite_vertices_begin();
         it != output.finite_vertices_end(); ++ it)
    {
      typename SMCDT::Edge_circulator circ = output.incident_edges(it), start = circ;
      do
      {
        if (output.is_constrained(*circ))
        {
          typename SMCDT::Vertex_handle vertex = it;
          typename SMCDT::Vertex_handle neighbor = circ->first->vertex((circ->second + 1)%3);
          if (neighbor == vertex)
            neighbor = circ->first->vertex((circ->second + 2)%3);

          output.insert (it,
                       typename GeomTraits::Point_3 (it->point().x(), it->point().y(),
                                                     std::numeric_limits<double>::quiet_NaN()),
                       typename GeomTraits::Direction_2
                       (typename GeomTraits::Vector_2 (it->point(), neighbor->point())));
        }
        
        ++ circ;
      }
      while (circ != start);
    }

    // Insert ridges in CDT
    for (std::size_t i = 0; i < lines.size(); ++ i)
    {
      if (!lines.is_ridge(i))
        continue;
      
      typename Detector::iterator it = lines.begin(i);
      typename SMCDT::Vertex_handle previous_cdt = output.insert (*it);

      ++ it;
      for (; it != lines.end(i); ++ it)
      {
        typename SMCDT::Vertex_handle current_cdt = output.insert (*it);
        output.insert_constraint (previous_cdt, current_cdt);
        previous_cdt = current_cdt;
      }
    }
    
    std::copy (input.planes().begin(), input.planes().end(),
               std::back_inserter (output.planes()));


    std::size_t nb_superfacets = output.segment_into_superfacets();
    TOP_VIEW_CERR << "  " << nb_superfacets << " superfacet(s) found" << std::endl;

    std::vector<std::vector<std::size_t> > plane_scores (nb_superfacets);
    

    typename SMCDT::Face_handle fhint;
    for (typename SMCDT::Finite_faces_iterator it = input.finite_faces_begin();
         it != input.finite_faces_end(); ++ it)
    {
      if (!input.is_pending(it) ||
          !it->info().has_plane())
        continue;

      fhint = output.locate (input.midpoint(it), fhint);

      if (output.is_default(fhint))
        continue;

      plane_scores[std::size_t(fhint->info().index)].push_back (it->info().plane_index);
    }

    std::vector<std::size_t> chosen_plane(nb_superfacets);
    for (std::size_t i = 0; i < plane_scores.size(); ++ i)
    {
      std::sort (plane_scores[i].begin(), plane_scores[i].end());

      std::size_t current = std::size_t(-1);
      std::size_t nb = 0;
      std::size_t nb_max = 0;
      std::size_t chosen = 0;
      for (std::size_t j = 0; j < plane_scores[i].size(); ++ j)
        if (plane_scores[i][j] != current)
        {
          if (nb > nb_max)
          {
            nb_max = nb;
            chosen = current;
          }
          nb = 0;
          current = plane_scores[i][j];
        }
        else
          ++ nb;
      if (nb > nb_max)
        chosen = current;
      CGAL_assertion (chosen != std::size_t(-1));
      chosen_plane[i] = chosen;
    }
    output.apply_planes_and_reset_superfacets (chosen_plane);
    
    // Estimate missing heights from planes
    std::ofstream file1 ("heights_1.xyz");
    file1.precision(18);
    std::ofstream file2 ("heights_2.xyz");
    file2.precision(18);
    for (typename SMCDT::Finite_vertices_iterator it = output.finite_vertices_begin();
         it != output.finite_vertices_end(); ++ it)
    {
      if (it->info().empty())
        output.insert (it);
                       
      std::vector<std::set<std::size_t> > section_incident_planes (it->info().size());
      
      typename SMCDT::Face_circulator circ = output.incident_faces(it);
      typename SMCDT::Face_circulator start = circ;
      do
      {
        if (!output.is_infinite(circ))
        {
          typename GeomTraits::Point_2 ref = output.midpoint(circ);
        
          std::size_t section = output.find_section_of_point_from_vertex_view (it, ref);

          if (!circ->info().has_plane())
            std::cerr << "Warning plane!" << std::endl;
          else
            section_incident_planes[section].insert (circ->info().plane_index);
        }
        ++ circ;
      }
      while (circ != start);

      for (std::size_t i = 0; i < section_incident_planes.size(); ++ i)
      {
        if (section_incident_planes[i].empty())
          continue;
        if (section_incident_planes[i].size() == 1)
        {
          double z = output.estimate_height_with_plane(it->point(), *section_incident_planes[i].begin());
          output.point(it, i) = typename GeomTraits::Point_3 (it->point().x(), it->point().y(), z);
          file1 << output.point(it, i) << std::endl;
        }
        else
        {
          double z = output.estimate_height_with_planes(it->point(), section_incident_planes[i], epsilon);
          output.point(it, i) = typename GeomTraits::Point_3 (it->point().x(), it->point().y(), z);
          file2 << output.point(it, i) << std::endl;
        }
      }
    }
    
    for (typename SMCDT::Finite_faces_iterator it = output.finite_faces_begin();
         it != output.finite_faces_end (); ++ it)
      if (!output.has_mesh_face(it))
      {
        typename SMCDT::Vertex_index v[3];
        bool okay = true;
        for (std::size_t i = 0; i < 3; ++ i)
        {
          if (output.has_unique_mesh_vertex(it->vertex(i)))
            v[i] = output.mesh_vertex (it->vertex(i));
          else if (output.has_mesh_vertex (it->vertex(i)))
          {
            typename GeomTraits::Point_2 ref = output.midpoint(it);
            std::size_t idx = output.find_section_of_point_from_vertex_view (it->vertex(i), ref);
            v[i] = output.mesh_vertex(it->vertex(i), idx);
          }
          else
          {
            okay = false;
            break;
          }

          if (!output.has_defined_height (v[i]))
          {
            okay = false;
            break;
          }
        }
        
        if (okay)
        {
          output.add_face(it, v[0], v[1], v[2]);
        }
        else
        {
          output.make_ignored (it);
        }
      }
  }
  
  template <typename GeomTraits>
  void region_growing_fill_voids (Surface_mesh_on_cdt<GeomTraits>& mesh)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    std::vector<typename GeomTraits::Plane_3>& planes = mesh.planes();

    std::set<std::pair<double, typename SMCDT::Edge> > todo;

    for (typename SMCDT::Finite_edges_iterator it = mesh.finite_edges_begin();
         it != mesh.finite_edges_end(); ++ it)
    {
      typename SMCDT::Face_handle fref = it->first;
      typename SMCDT::Face_handle fcurrent = it->first->neighbor (it->second);

      if (!mesh.is_pending(fref) || !mesh.is_pending(fcurrent))
        continue;

      if (fcurrent->info().has_plane() == fref->info().has_plane())
        continue;
      
      if (fcurrent->info().has_plane())
        std::swap (fref, fcurrent);

      typename GeomTraits::Vector_3 v0 = mesh.normal_vector(fcurrent);
      typename GeomTraits::Vector_3 v1 = planes[fref->info().plane_index].orthogonal_vector();
      todo.insert (std::make_pair (CGAL::abs(v0 * v1), *it));
    }

    while (!todo.empty())
    {
      typename SMCDT::Edge current = todo.begin()->second;
      todo.erase (todo.begin());

      typename SMCDT::Face_handle fref = current.first;
      typename SMCDT::Face_handle fcurrent = current.first->neighbor(current.second);
      if (fcurrent->info().has_plane())
      {
        std::swap (fref, fcurrent);
        if (fcurrent->info().has_plane())
          continue;
      }

      fcurrent->info().plane_index = fref->info().plane_index;

      for (std::size_t i = 0; i < 3; ++ i)
      {
        typename SMCDT::Face_handle neigh = fcurrent->neighbor(i);

        if (!mesh.is_pending(neigh))
          continue;

        if (!neigh->info().has_plane())
        {
          typename GeomTraits::Vector_3 v0 = mesh.normal_vector(neigh);
          typename GeomTraits::Vector_3 v1 = planes[fref->info().plane_index].orthogonal_vector();
          todo.insert (std::make_pair (CGAL::abs(v0 * v1), std::make_pair (fcurrent, i)));
        }
      }
    }
  }
    
  
  template <typename GeomTraits>
  void region_growing_fill_voids_on_meshless_faces (Surface_mesh_on_cdt<GeomTraits>& mesh)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin();
         it != mesh.finite_faces_end(); ++ it)
    {
      if (mesh.is_ignored(it) || mesh.has_mesh_face(it))
        continue;

      std::set<typename SMCDT::Face_handle> done;
      
      std::queue<typename SMCDT::Face_handle> todo;
      todo.push(it);
      while (!todo.empty())
      {
        typename SMCDT::Face_handle current = todo.front();
        todo.pop();

        if (mesh.has_mesh_face(current))
        {
          it->info().plane_index = current->info().plane_index;
          break;
        }

        for (std::size_t i = 0; i < 3; ++ i)
          if (!mesh.is_constrained (std::make_pair (current, i))
              && done.insert (current->neighbor(i)).second)
            todo.push (current->neighbor(i));
      }

      if (!it->info().has_plane())
      {
        std::ofstream f("face.xyz");
        f.precision(18);
        for (std::size_t i = 0; i < 3; ++ i)
          f << it->vertex(i)->point() << " 0" << std::endl;
        abort();
      }
    }    
  }
    
  template <typename GeomTraits>
  void region_growing_one_pass (Surface_mesh_on_cdt<GeomTraits>& mesh,
                                double epsilon,
                                std::size_t nb_min,
                                double angle_max)
  {  
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
    
    std::vector<typename GeomTraits::Plane_3>& planes = mesh.planes();

    std::vector<typename SMCDT::Face_handle> faces;
    
    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin();
         it != mesh.finite_faces_end(); ++ it)
      if (mesh.is_pending(it) && !it->info().has_plane())
        faces.push_back (it);

    std::sort (faces.begin(), faces.end(), typename SMCDT::Sort_faces_by_planarity(&mesh));

    std::size_t nb_failed_trials = 0;

    for (std::size_t fi = 0; fi < faces.size(); ++ fi)
    {
      if (faces[fi]->info().has_plane())
        continue;

      faces[fi]->info().plane_index = planes.size();
      
      //characteristics of the seed
      typename GeomTraits::Vector_3 normal_seed = mesh.normal_vector (faces[fi]);
      if (nb_min == 1)
        normal_seed = typename GeomTraits::Vector_3(0., 0., 1.);
      
      typename GeomTraits::Point_3 pt_seed = mesh.midpoint_3(faces[fi]);
      typename GeomTraits::Plane_3 optimal_plane (pt_seed, normal_seed);

      //initialization containers
      std::vector<typename SMCDT::Face_handle> index_container (1, faces[fi]);
      std::vector<typename SMCDT::Face_handle> index_container_former_ring (1, faces[fi]);
      std::list<typename SMCDT::Face_handle> index_container_current_ring;

      std::set<typename SMCDT::Face_handle> current_overlaps;
      
      std::vector<typename GeomTraits::Triangle_3> support (1, mesh.triangle_3(faces[fi]));
      //propagation
      bool propagation = true;
      do{

        propagation = false;

        for (std::size_t k = 0; k < index_container_former_ring.size(); k++)
        {
          typename SMCDT::Face_handle current = index_container_former_ring[k];
          for (std::size_t i = 0; i < 3; ++ i)
          {
            typename SMCDT::Face_handle neighbor = current->neighbor(i);
            if (!mesh.is_pending (neighbor)
                || neighbor->info().has_plane())
              continue;

            typename GeomTraits::Vector_3 normal = mesh.normal_vector (neighbor);
            if (std::fabs(normal * optimal_plane.orthogonal_vector()) > angle_max)
            {
              typename GeomTraits::Triangle_3 candidate = mesh.triangle_3 (neighbor);
              bool okay = true;
              
              for (std::size_t j = 0; j < 3; ++ j)
                if (CGAL::squared_distance (optimal_plane, candidate[j]) > epsilon * epsilon)
                {
                  okay = false;
                  break;
                }
              if (!okay)
                continue;

              neighbor->info().plane_index = planes.size();
              propagation = true;
              index_container_current_ring.push_back(neighbor);
              support.push_back (mesh.triangle_3 (neighbor));
            }
          }
        }
			
        //update containers
        index_container_former_ring.clear();
        BOOST_FOREACH (typename SMCDT::Face_handle fh, index_container_current_ring)
        {
          index_container_former_ring.push_back(fh);
          index_container.push_back(fh);
        }
        index_container_current_ring.clear();

        if (nb_min != 1)
          CGAL::linear_least_squares_fitting_3 (support.begin(), support.end(),
                                                optimal_plane,
                                                CGAL::Dimension_tag<2>());

      } while (propagation);


      //Test the number of inliers -> reject if inferior to Nmin
      if (index_container.size() < nb_min)
      {
        faces[fi]->info().erase_plane();
        
        for (std::size_t k = 0; k < index_container.size(); k++)
          index_container[k]->info().erase_plane();
        ++ nb_failed_trials;
      }
      else
        planes.push_back (optimal_plane);
    }

    TOP_VIEW_CERR << planes.size() << " planes detected ("
                  << nb_failed_trials << " failed trials)" << std::endl;
  }

  template <typename GeomTraits>
  void region_growing_create_ridge_buffer (Surface_mesh_on_cdt<GeomTraits>& mesh,
                                           std::size_t nb_min)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
    
    std::vector<typename SMCDT::Face_handle> faces_with_plane;
    
    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end(); ++ it)
    {
      if (!mesh.has_unique_mesh_vertex(it))
        continue;

      bool is_ridge = false;

      std::size_t plane_idx = std::size_t(-1);

      faces_with_plane.clear();
      
      typename SMCDT::Face_circulator circ = mesh.incident_faces(it);
      typename SMCDT::Face_circulator start = circ;
      do
      {
        if (circ->info().has_plane())
        {
          std::size_t idx = circ->info().plane_index;
          if (plane_idx == std::size_t(-1))
            plane_idx = idx;
          else if (plane_idx != idx)
            is_ridge = true;
          faces_with_plane.push_back (circ);
        }
        ++ circ;
      }
      while (circ != start);

      if (is_ridge)
        for (std::size_t i = 0; i < faces_with_plane.size(); ++ i)
          mesh.make_ridge_buffer(faces_with_plane[i]);
    }

    std::set<typename SMCDT::Face_handle> done;
    
    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin();
         it != mesh.finite_faces_end(); ++ it)
    {
      if (mesh.is_buffer(it) ||
          !it->info().has_plane())
        continue;

      std::vector<typename SMCDT::Face_handle> region;
      std::queue<typename SMCDT::Face_handle> todo;
      todo.push(it);
      while (!todo.empty())
      {
        typename SMCDT::Face_handle current = todo.front();
        todo.pop();

        if (!done.insert(current).second)
          continue;

        region.push_back (current);

        for (std::size_t i = 0; i < 3; ++ i)
        {
          if (!mesh.is_buffer(current->neighbor(i)) &&
              current->info().plane_index == current->neighbor(i)->info().plane_index)
            todo.push(current->neighbor(i));
        }
      }

      if (region.size() < nb_min)
        for (std::size_t i = 0; i < region.size(); ++ i)
          mesh.make_ridge_buffer(region[i]);
    }

  }
  
  template <typename GeomTraits>
  void region_growing (Surface_mesh_on_cdt<GeomTraits>& mesh,
                       double epsilon,
                       std::size_t nb_min,
                       double angle_max,
                       bool fill_voids_on_meshless_faces = true)
  {

    region_growing_one_pass<GeomTraits> (mesh, epsilon, nb_min, angle_max);
    region_growing_fill_voids<GeomTraits> (mesh);
    region_growing_one_pass<GeomTraits> (mesh, std::numeric_limits<double>::max(), 1, 0.);

    mesh.DEBUG_dump_ply("debug0.ply");
    
    if (fill_voids_on_meshless_faces)
      region_growing_fill_voids_on_meshless_faces<GeomTraits> (mesh);
    else
      region_growing_create_ridge_buffer<GeomTraits> (mesh, nb_min / 4);
  }

  template <typename GeomTraits>
  void project_points_on_detected_planes (Surface_mesh_on_cdt<GeomTraits>& mesh,
                                          double epsilon)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end(); ++ it)
    {
      if (!mesh.has_unique_mesh_vertex(it))
        continue;
      typename SMCDT::Vertex_index vi = mesh.mesh_vertex(it);

      std::set<std::size_t> incident_planes;
      typename SMCDT::Line_3 line (mesh.point(vi), typename GeomTraits::Vector_3 (0., 0., 1.));

      typename SMCDT::Face_circulator circ = mesh.incident_faces(it), start = circ;
      do
      {
        if (circ->info().has_plane())
          incident_planes.insert (circ->info().plane_index);
        ++ circ;
      }
      while (circ != start);

      typename GeomTraits::Point_3 point
        = mesh.barycenter_of_projections_on_planes (it, incident_planes);
      
      if (incident_planes.size() > 2)
      {
        typename GeomTraits::Point_3 candidate = mesh.least_squares_plane_intersection (incident_planes);
        if (CGAL::squared_distance (typename GeomTraits::Point_2 (point.x(), point.y()),
                                    typename GeomTraits::Point_2 (candidate.x(), candidate.y()))
            < epsilon * epsilon)
          point = candidate;
      }
      else if (incident_planes.size() == 2)
      {
        typename std::set<std::size_t>::iterator sit = incident_planes.begin();
        ++ sit;
        typename GeomTraits::Line_3 line;
        if (mesh.intersection_line_of_2_planes (*(incident_planes.begin()), *sit, line))
        {
          typename GeomTraits::Point_3 candidate = line.projection (point);
          if (CGAL::squared_distance (typename GeomTraits::Point_2 (point.x(), point.y()),
                                      typename GeomTraits::Point_2 (candidate.x(), candidate.y()))
              < epsilon * epsilon)
            point = candidate;
        }
      }
      mesh.point(vi) = point;
    }

  }

  template <typename GeomTraits>
  void project_points_on_detected_planes (Surface_mesh_on_cdt<GeomTraits>& mesh)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    std::vector<typename GeomTraits::Plane_3>& planes = mesh.planes();
      
    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end(); ++ it)
    {
      if (!mesh.has_unique_mesh_vertex(it))
        continue;
      typename SMCDT::Vertex_index vi = mesh.mesh_vertex(it);

      std::set<std::size_t> incident_planes;
      typename SMCDT::Line_3 line (mesh.point(vi), typename GeomTraits::Vector_3 (0., 0., 1.));

      typename SMCDT::Face_circulator circ = mesh.incident_faces(it), start = circ;
      do
      {
        if (circ->info().has_plane())
          incident_planes.insert (circ->info().plane_index);
        ++ circ;
      }
      while (circ != start);

      typename GeomTraits::Point_3 new_point;
      std::size_t nb = 0;
      BOOST_FOREACH (std::size_t idx, incident_planes)
      {
        const typename GeomTraits::Plane_3& plane = planes[idx];
        typename CGAL::cpp11::result_of<typename GeomTraits::Intersect_3
                                        (typename GeomTraits::Line_3,
                                         typename GeomTraits::Plane_3)>::type
          result = CGAL::intersection(line, plane);
        typename GeomTraits::Point_3* inter;
        if (result && (inter = boost::get<typename GeomTraits::Point_3>(&*result)))
        {
          new_point = CGAL::barycenter (*inter, 1, new_point, nb);
          ++ nb;
        }
      }
      if (nb != 0)
        mesh.point(vi) = new_point;
    }

  }
                         
  template <typename GeomTraits>
  void generate_missing_3d_points_from_planes (Surface_mesh_on_cdt<GeomTraits>& mesh)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end(); ++ it)
    {
      if (mesh.has_mesh_vertex(it))
        continue;

      typename SMCDT::Edge_circulator circ = mesh.incident_edges(it), start = circ;
      do
      {
        if (mesh.is_constrained(*circ))
        {
          typename SMCDT::Vertex_handle vertex = it;
          typename SMCDT::Vertex_handle neighbor = circ->first->vertex((circ->second + 1)%3);
          if (neighbor == vertex)
            neighbor = circ->first->vertex((circ->second + 2)%3);

          mesh.insert (it,
                       typename GeomTraits::Point_3 (it->point().x(), it->point().y(),
                                                     std::numeric_limits<double>::quiet_NaN()),
                       typename GeomTraits::Direction_2
                       (typename GeomTraits::Vector_2 (it->point(), neighbor->point())));
        }
        
        ++ circ;
      }
      while (circ != start);

      mesh.estimate_missing_heights_from_planes(it);
    }
  }
  
  template <typename GeomTraits>
  void generate_missing_3d_points (Surface_mesh_on_cdt<GeomTraits>& mesh,
                                   double epsilon)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    std::ofstream f1 ("f1.xyz");
    f1.precision(18);
    std::ofstream f2 ("f2.xyz");
    f2.precision(18);
    std::ofstream f3 ("f3.xyz");
    f3.precision(18);
    
    TOP_VIEW_CERR << "  First pass: propagating heights of neighbor 3D points" << std::endl;
    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end(); ++ it)
    {
      if (mesh.has_mesh_vertex(it))
        continue;

      typename SMCDT::Edge_circulator circ = mesh.incident_edges(it), start = circ;
      do
      {
        if (mesh.is_constrained(*circ))
        {
          typename SMCDT::Vertex_handle vertex = it;
          typename SMCDT::Vertex_handle neighbor = circ->first->vertex((circ->second + 1)%3);
          if (neighbor == vertex)
            neighbor = circ->first->vertex((circ->second + 2)%3);

          mesh.insert (it,
                       typename GeomTraits::Point_3 (it->point().x(), it->point().y(),
                                                     std::numeric_limits<double>::quiet_NaN()),
                       typename GeomTraits::Direction_2
                       (typename GeomTraits::Vector_2 (it->point(), neighbor->point())));
        }
        
        ++ circ;
      }
      while (circ != start);

      mesh.estimate_missing_heights(it, epsilon);
      for (std::size_t i = 0; i < mesh.number_of_mesh_vertices (it); ++ i)
      {
        const typename GeomTraits::Point_3& p = mesh.point (it, i);
        if (p.z() == p.z())
          f1 << p << std::endl;
      }
    }

    TOP_VIEW_CERR << "  Second pass: if no 3D point close, propagate heights of neighbor borders" << std::endl;
    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end(); ++ it)
    {
      bool to_do = false;
      std::vector<bool> todo (mesh.number_of_mesh_vertices(it), false);
      for (std::size_t i = 0; i < mesh.number_of_mesh_vertices(it); ++ i)
        if (!mesh.has_defined_height(mesh.mesh_vertex(it, i)))
        {
          to_do = true;
          todo[i] = true;
//          break;
        }
      if (to_do)
        mesh.estimate_missing_heights(it, epsilon, true);

      for (std::size_t i = 0; i < mesh.number_of_mesh_vertices (it); ++ i)
      {
        if (todo[i])
        {
          const typename GeomTraits::Point_3& p = mesh.point (it, i);
          if (p.z() == p.z())
            f2 << p << std::endl;
        }
      }
    }

    TOP_VIEW_CERR << "  Third pass: if no point close, propagate from closest" << std::endl;
    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end(); ++ it)
    {
      bool to_do = false;
      std::vector<bool> todo (mesh.number_of_mesh_vertices(it), false);

      for (std::size_t i = 0; i < mesh.number_of_mesh_vertices(it); ++ i)
        if (!mesh.has_defined_height(mesh.mesh_vertex(it, i)))
        {
          to_do = true;
          todo[i] = true;
//          break;
        }
      
      if (to_do)
        mesh.estimate_missing_heights_no_limit(it);
      for (std::size_t i = 0; i < mesh.number_of_mesh_vertices (it); ++ i)
      {
        if (todo[i])
        {
          typename GeomTraits::Point_3& p = mesh.point (it, i);
          if (p.z() == p.z())
            f3 << p << std::endl;
          else
          {
            mesh.remove (mesh.mesh_vertex(it, i));
            todo.erase(todo.begin() + i);
            -- i;
          }
        }
      }
    }

    mesh.check_structure_integrity();
  }
  
  template <typename GeomTraits>
  void generate_missing_3d_faces (Surface_mesh_on_cdt<GeomTraits>& mesh)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
    
    TOP_VIEW_CERR << "  Generating faces" << std::endl;
    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin();
         it != mesh.finite_faces_end (); ++ it)
      if (!mesh.has_mesh_face(it) && !mesh.is_ignored(it))
      {
        typename SMCDT::Vertex_index v[3];
        bool okay = true;
        for (std::size_t i = 0; i < 3; ++ i)
        {
          if (mesh.has_unique_mesh_vertex(it->vertex(i)))
            v[i] = mesh.mesh_vertex (it->vertex(i));
          else if (mesh.has_mesh_vertex (it->vertex(i)))
          {
            typename GeomTraits::Point_2 ref = mesh.midpoint(it);
            std::size_t idx = mesh.find_section_of_point_from_vertex_view (it->vertex(i), ref);
            v[i] = mesh.mesh_vertex(it->vertex(i), idx);
          }
          else
          {
            okay = false;
            break;
          }

          if (!mesh.has_defined_height (v[i]))
          {
            okay = false;
            break;
          }
        }
        
        if (okay)
        {
          mesh.add_face(it, v[0], v[1], v[2]);
        }
        else
        {
          mesh.make_ignored (it);
        }
      }
  }

  template <typename GeomTraits>
  void snap_intersecting_borders (Surface_mesh_on_cdt<GeomTraits>& mesh)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    std::ofstream file ("inter.polylines.txt");
    file.precision(18);

    std::queue<typename SMCDT::Edge> todo;
    for (typename SMCDT::Finite_edges_iterator it = mesh.finite_edges_begin();
         it != mesh.finite_edges_end(); ++ it)
      if (mesh.is_constrained(*it))
        todo.push (*it);

    while (!todo.empty())
    {
      int idx = todo.front().second;
      typename SMCDT::Face_handle
        fe = todo.front().first,
        fo = fe->neighbor (idx);

      todo.pop();
      if (!mesh.has_mesh_face (fe) || !mesh.has_mesh_face (fo))
        continue;
      
      typename GeomTraits::Point_2
        re = mesh.midpoint(fe),
        ro = mesh.midpoint(fo);

      typename SMCDT::Vertex_handle
        va = fe->vertex ((idx + 1)%3),
        vb = fe->vertex ((idx + 2)%3);

      std::size_t
        sea = mesh.find_section_of_point_from_vertex_view (va, re),
        seb = mesh.find_section_of_point_from_vertex_view (vb, re),
        soa = mesh.find_section_of_point_from_vertex_view (va, ro),
        sob = mesh.find_section_of_point_from_vertex_view (vb, ro);

      typename GeomTraits::Point_3
        &pea = mesh.point(va, sea),
        &peb = mesh.point(vb, seb),
        &poa = mesh.point(va, soa),
        &pob = mesh.point(vb, sob);

      if (pea == poa || peb == pob)
        continue;

      if ((pea.z() - poa.z()) * (peb.z() - pob.z()) > 0.) // Do borders intersect
        continue; // If no, no problem, continue

      file << "2 " << pea << " " << peb << std::endl
           << "2 " << poa << " " << pob << std::endl;

      typename SMCDT::Vertex_handle merged = va;
      typename SMCDT::Vertex_handle other = vb;
      if (CGAL::abs(pea.z() - poa.z()) < CGAL::abs(peb.z() - pob.z())) // snap to A
        mesh.merge_vertices (mesh.mesh_vertex (va, sea),
                             mesh.mesh_vertex (va, soa));
      else // snap to B
      {
        mesh.merge_vertices (mesh.mesh_vertex (vb, seb),
                             mesh.mesh_vertex (vb, sob));
        merged = vb;
        other = va;
      }

      // Recheck neighbor edges
      std::vector<typename SMCDT::Edge> incident;
      mesh.incident_constraints (merged, std::back_inserter (incident));
      if (incident.size() == 1)
        continue;

      for (std::size_t i = 0; i < incident.size(); ++ i)
        if (incident[i].first->vertex ((incident[i].second + 1)%3) != other &&
            incident[i].first->vertex ((incident[i].second + 2)%3) != other)
          todo.push (incident[i]);
    }
  }
  
  template <typename GeomTraits>
  void generate_vertical_walls (Surface_mesh_on_cdt<GeomTraits>& mesh,
                                double epsilon)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    std::ofstream file ("vertical.xyz");
    file.precision(18);

    TOP_VIEW_CERR << "Points" << std::endl;

    // Generated support vertices
    for (typename SMCDT::Finite_vertices_iterator it = mesh.finite_vertices_begin();
         it != mesh.finite_vertices_end(); ++ it)
    {
      if (mesh.number_of_mesh_vertices (it) == 1)
      {
        mesh.set_next_vertex (mesh.mesh_vertex(it), mesh.mesh_vertex(it));
        continue;
      }
      if (mesh.number_of_mesh_vertices (it) == 0
          || !mesh.is_valid (it))
        continue;

      // Get existing heights
      std::vector<std::pair<double, int> > heights;
      heights.reserve (it->info().size());
      for (std::size_t i = 0; i < it->info().size(); ++ i)
        heights.push_back (std::make_pair (mesh.point (it, i).z(), int(i)));

      CGAL_assertion (heights.size() > 1);
      
      std::sort (heights.begin(), heights.end());

      // Refine
      std::size_t size = heights.size() - 1;
      for (std::size_t i = 0; i < size; ++ i)
      {
        double h0 = heights[i].first;
        double h1 = heights[i + 1].first;
        double dh = CGAL::abs (h0 - h1);

        if (dh < epsilon * 1.5)
          continue;

        std::size_t nb_pts = std::size_t (dh / epsilon);
        for (std::size_t j = 1; j < nb_pts; ++ j)
          heights.push_back (std::make_pair (h0 + (h1 - h0) * (j / double(nb_pts)), -1));
      }

      std::sort (heights.begin(), heights.end());

      CGAL_assertion (heights.back().second != -1);
        
      typename SMCDT::Vertex_index latest_vertex = mesh.mesh_vertex (it, std::size_t(heights.back().second));
      mesh.set_next_vertex (latest_vertex, latest_vertex);
      for (std::size_t i = 0; i < heights.size(); ++ i)
      {
        typename SMCDT::Vertex_index vertex;
        if (heights[i].second != -1)
          vertex = mesh.mesh_vertex (it, heights[i].second);
        else
        {
          typename GeomTraits::Point_3 new_point (it->point().x(), it->point().y(), heights[i].first);
          vertex = mesh.insert (it, new_point);
        }
        file << mesh.point (vertex) << std::endl;
        
#ifdef TOP_VIEW_FIX_DUPLICATE_VERTICES
        if (latest_vertex != vertex)
          mesh.set_next_vertex (latest_vertex, vertex);
#else
        mesh.set_next_vertex (latest_vertex, vertex);
#endif
        latest_vertex = vertex;
      }
    }

    mesh.check_structure_integrity();

    TOP_VIEW_CERR << "Faces" << std::endl;
    
    std::ofstream file2 ("faces.xyz");
    file2.precision(18);

    std::ofstream file3 ("notedge.polylines.txt");
    file3.precision(18);

    // Generate walls
    for (typename SMCDT::Finite_edges_iterator it = mesh.finite_edges_begin();
         it != mesh.finite_edges_end(); ++ it)
    {
      if (!mesh.is_constrained(*it))
        continue;

      int idx = it->second;

      typename SMCDT::Face_handle
        fe = it->first,
        fo = fe->neighbor (idx);

      if (!mesh.has_mesh_face (fe) || !mesh.has_mesh_face (fo))
        continue;

      typename GeomTraits::Point_2
        re = mesh.midpoint(fe),
        ro = mesh.midpoint(fo);

      typename SMCDT::Vertex_handle
        va = fe->vertex ((idx + 1)%3),
        vb = fe->vertex ((idx + 2)%3);

      if (!mesh.is_valid (va) || !mesh.is_valid (vb))
        continue;

      std::size_t
        sea = mesh.find_section_of_point_from_vertex_view (va, re),
        seb = mesh.find_section_of_point_from_vertex_view (vb, re),
        soa = mesh.find_section_of_point_from_vertex_view (va, ro),
        sob = mesh.find_section_of_point_from_vertex_view (vb, ro);

      typename GeomTraits::Point_3
        &pea = mesh.point(va, sea),
        &peb = mesh.point(vb, seb),
        &poa = mesh.point(va, soa),
        &pob = mesh.point(vb, sob);

      typename GeomTraits::Plane_3 plane (pea, peb, poa);
      if (pea == poa)
      {
        plane = typename GeomTraits::Plane_3 (peb, pea, pob);
        if (peb == pob)
          continue;
      }

      if (mesh.next_vertex (mesh.mesh_vertex(va)) == typename SMCDT::Vertex_index()
          || mesh.next_vertex (mesh.mesh_vertex(vb)) == typename SMCDT::Vertex_index())
      {
        std::ofstream ff ("pts.xyz");
        ff.precision(18);
        ff << pea << std::endl
           << peb << std::endl
           << poa << std::endl
           << pob << std::endl;
        std::cerr << "EXIT" << std::endl;
        abort();
      }
      
      double hamin = (std::min)(pea.z(), poa.z());
      double hamax = (std::max)(pea.z(), poa.z());
      double hbmin = (std::min)(peb.z(), pob.z());
      double hbmax = (std::max)(peb.z(), pob.z());

      typename SMCDT::CDT cdt;
      typename SMCDT::Vertex_handle vva = typename SMCDT::Vertex_handle();
      typename SMCDT::Vertex_handle vvb = typename SMCDT::Vertex_handle();

      typename SMCDT::Vertex_index start = mesh.mesh_vertex (va);
      typename SMCDT::Vertex_index circ = start;
      do
      {
        typename GeomTraits::Point_3 point_3 = mesh.point (circ);
        if (hamin <= point_3.z() && point_3.z() <= hamax)
        {
          typename GeomTraits::Point_2 point_2 = plane.to_2d (point_3);
          typename SMCDT::Vertex_handle vh = cdt.insert (point_2);
          vh->info().push_back (std::make_pair (typename GeomTraits::Direction_2 (0., 0.), circ));

          if (point_3.z() == pea.z())
            vva = vh;
        }

        circ = mesh.next_vertex (circ);
      }
      while (start != circ);

      start = mesh.mesh_vertex (vb);
      circ = start;
      do
      {
        typename GeomTraits::Point_3 point_3 = mesh.point (circ);
        if (hbmin <= point_3.z() && point_3.z() <= hbmax)
        {
          typename GeomTraits::Point_2 point_2 = plane.to_2d (point_3);
          typename SMCDT::Vertex_handle vh = cdt.insert (point_2);
          vh->info().push_back (std::make_pair (typename GeomTraits::Direction_2 (0., 0.), circ));

          if (point_3.z() == peb.z())
            vvb = vh;
        }
        
        circ = mesh.next_vertex (circ);
      }
      while (start != circ);


      bool inverse_orientation = false;        
      if (vva != typename SMCDT::Vertex_handle() &&
          vvb != typename SMCDT::Vertex_handle() &&
          vva != vvb)
      {
        typename SMCDT::Edge border_edge;

        if (cdt.is_edge (vva, vvb, border_edge.first, border_edge.second))
        {

          if (cdt.is_infinite(border_edge.first))
            border_edge = cdt.mirror_edge(border_edge);

          int idx1 = fe->index (va);
          int idx2 = fe->index (vb);
          int border_idx1 = border_edge.first->index(vva);
          int border_idx2 = border_edge.first->index(vvb);

          if ((cdt.cw(idx1) == idx2 && cdt.cw(border_idx1) == border_idx2) ||
              (cdt.ccw(idx1) == idx2 && cdt.ccw(border_idx1) == border_idx2))
            inverse_orientation = true;
        }
        else
        {
          TOP_VIEW_CERR << "Warning: vertices for orientation test are not on an edge" << std::endl;
          file3 << "2 " << mesh.point(vva->info()[0].second) << " "
                << mesh.point(vvb->info()[0].second) << std::endl;
        }
      }
      else
        TOP_VIEW_CERR << "Warning: vertices for orientation test not found" << std::endl;
      
      for (typename SMCDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
           fit != cdt.finite_faces_end(); ++ fit)
      {
        typename SMCDT::Vertex_index
          v0 = fit->vertex(0)->info()[0].second,
          v1 = fit->vertex(1)->info()[0].second,
          v2 = fit->vertex(2)->info()[0].second;
        if (inverse_orientation)
          std::swap (v0, v1);
        if (!mesh.add_face (v0, v1, v2))
          for (std::size_t k = 0; k < 3; ++ k)
            file2 << mesh.point(fit->vertex(k)->info()[0].second) << std::endl;
      }

    }

    TOP_VIEW_CERR << "DONE" << std::endl;
  }
} // namespace internal


  struct Parameters
  {
    double epsilon;
    double quantile;
    bool estimate_borders_from_planes;
    bool flatten_planar_regions;

    Parameters(double epsilon,
               double quantile = 0.5,
               bool estimate_borders_from_planes = true,
               bool flatten_planar_regions = true)
      : epsilon (epsilon)
      , quantile (quantile)
      , estimate_borders_from_planes (estimate_borders_from_planes)
      , flatten_planar_regions (flatten_planar_regions)
    { }
  };
  
} // namespace Top_view_surface_reconstruction_3

#if 0
template <typename GeomTraits, typename PointInputIterator, typename PointMap>
void top_view_surface_reconstruction_old (PointInputIterator begin,
                                      PointInputIterator end,
                                      PointMap point_map,
                                      Surface_mesh_on_cdt<GeomTraits>& output_mesh,
                                      const Top_view_surface_reconstruction_3::Parameters& parameters)
{
//  typedef typename GeomTraits::Point_3 Point_3;
  typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
  typedef Border_graph<GeomTraits> Graph;
  typedef Line_detector<GeomTraits> Detect;

  CGAL_assertion (parameters.epsilon > 0.);
  
  double threshold_factor = 2.;
  double meshing_factor = std::sqrt(2.);
  
  SMCDT tmp_mesh;

  TOP_VIEW_CERR << "Inserting filtered points" << std::endl;
  Top_view_surface_reconstruction_3::internal::insert_filtered_points<GeomTraits>
    (begin, end, point_map, tmp_mesh, parameters.epsilon, parameters.quantile);

  TOP_VIEW_CERR << "Filtering faces" << std::endl;
  Top_view_surface_reconstruction_3::internal::filter_faces<GeomTraits>
    (tmp_mesh, parameters.epsilon * threshold_factor);

#ifdef TOP_VIEW_LOG
  tmp_mesh.DEBUG_dump_off_0();
#endif

  TOP_VIEW_CERR << "Cleaning buffer zone" << std::endl;
  Top_view_surface_reconstruction_3::internal::cleanup_buffer_zone<GeomTraits>
    (tmp_mesh);

  Graph graph;

  TOP_VIEW_CERR << "Building border graph" << std::endl;
  Top_view_surface_reconstruction_3::internal::build_border_graph<GeomTraits>
    (tmp_mesh, graph);

  graph.filter_small_terminal_borders(10. * parameters.epsilon * threshold_factor);
  graph.split_into_polylines();


#ifdef TOP_VIEW_LOG
  tmp_mesh.DEBUG_dump_poly();
  graph.DEBUG_dump_poly("poly_complex.polylines.txt");
#endif
  
  TOP_VIEW_CERR << "Simplifying border graph" << std::endl;
  Top_view_surface_reconstruction_3::internal::simplify_border_graph<GeomTraits>
    (tmp_mesh, graph, parameters.epsilon * threshold_factor);

#ifdef TOP_VIEW_LOG
  graph.DEBUG_dump_poly("poly_simplified.polylines.txt");
#endif

#ifdef TOP_VIEW_LOG
  tmp_mesh.DEBUG_dump_ply_3();
#endif

  TOP_VIEW_CERR << "Detecting lines in buffer" << std::endl;
  Detect detect (tmp_mesh);

  detect.sort_candidates_with_graph (graph);
  detect.run (parameters.epsilon /* * meshing_factor*/);

#ifdef TOP_VIEW_LOG
  detect.DEBUG_dump_ply("buffer_lines.ply");
  detect.DEBUG_dump_polyline("detected.polylines.txt");
  tmp_mesh.DEBUG_dump_off_2();
#endif

  TOP_VIEW_CERR << "Creating mesh with borders" << std::endl;
  Top_view_surface_reconstruction_3::internal::create_mesh_with_borders<GeomTraits>
    (tmp_mesh, detect, output_mesh, parameters.epsilon * meshing_factor,
     parameters.epsilon * threshold_factor * meshing_factor,
     parameters.estimate_borders_from_planes);

#ifdef TOP_VIEW_LOG
  output_mesh.DEBUG_dump_off_4();
  output_mesh.DEBUG_dump_off_1();
#endif

  if (parameters.estimate_borders_from_planes || parameters.flatten_planar_regions)
  {
    TOP_VIEW_CERR << "Detecting planar regions" << std::endl;
    Top_view_surface_reconstruction_3::internal::region_growing<GeomTraits>
      (output_mesh, parameters.epsilon, 20, 0.9);

#ifdef TOP_VIEW_LOG
    output_mesh.DEBUG_dump_ply_1();
    output_mesh.DEBUG_dump_ply_2();
    output_mesh.DEBUG_dump_planes_inter();
#endif

  }

  
  if (parameters.flatten_planar_regions)
  {
    TOP_VIEW_CERR << "Flattening planar regions" << std::endl;
    Top_view_surface_reconstruction_3::internal::project_points_on_detected_planes<GeomTraits>
      (output_mesh, parameters.epsilon);
  }

  output_mesh.check_structure_integrity();

  {
    std::ofstream file("slivers.off");
    file.precision(18);
    std::vector<typename GeomTraits::Triangle_2> tri;
    for (typename SMCDT::Finite_faces_iterator it = output_mesh.finite_faces_begin();
         it != output_mesh.finite_faces_end(); ++ it)
    {
      if (output_mesh.is_ignored(it))
        continue;

      bool sliver = false;
      for (std::size_t i = 0; i < 3; ++ i)
      {
        typename GeomTraits::Point_2 pt = it->vertex(i)->point();
        typename GeomTraits::Line_2 l (it->vertex((i+1)%3)->point(),
                                       it->vertex((i+2)%3)->point());
        if (CGAL::squared_distance(pt, l) < 0.000001)
          sliver = true;
      }

      if (sliver)
        tri.push_back (typename GeomTraits::Triangle_2
                       (it->vertex(0)->point(),
                        it->vertex(1)->point(),
                        it->vertex(2)->point()));
    }

    file << "OFF" << std::endl
         << tri.size() * 3 << " " << tri.size() << " 0" << std::endl;
    for (std::size_t i = 0; i < tri.size(); ++ i)
      for (std::size_t j = 0; j < 3; ++ j)
        file << tri[i][j] << " 0" << std::endl;
    for (std::size_t i = 0; i < tri.size(); ++ i)
      file << "3 " << 3*i << " " << 3*i + 1 << " " << 3*i + 2 << std::endl;
  }

  TOP_VIEW_CERR << "Generating missing 3D points" << std::endl;
  if (parameters.estimate_borders_from_planes)
    Top_view_surface_reconstruction_3::internal::generate_missing_3d_points_from_planes<GeomTraits>
      (output_mesh);
  else
    Top_view_surface_reconstruction_3::internal::generate_missing_3d_points<GeomTraits>
      (output_mesh, parameters.epsilon * meshing_factor);

  Top_view_surface_reconstruction_3::internal::generate_missing_3d_faces<GeomTraits>
    (output_mesh);
  
  TOP_VIEW_CERR << "Snapping intersecting borders" << std::endl;
  Top_view_surface_reconstruction_3::internal::snap_intersecting_borders<GeomTraits>
    (output_mesh);
  output_mesh.check_structure_integrity();

  TOP_VIEW_CERR << "Generating vertical walls" << std::endl;
  Top_view_surface_reconstruction_3::internal::generate_vertical_walls<GeomTraits>
    (output_mesh, parameters.epsilon * meshing_factor);

#ifdef TOP_VIEW_LOG
  output_mesh.DEBUG_dump_edges();
#endif

  
#ifndef TOP_VIEW_FIX_DUPLICATE_VERTICES
  output_mesh.stitch_borders();
#endif

  output_mesh.check_structure_integrity();
#ifdef TOP_VIEW_LOG
  output_mesh.DEBUG_dump_off_5();
#endif

}
#endif




template <typename GeomTraits, typename PointInputIterator, typename PointMap>
void top_view_surface_reconstruction (PointInputIterator begin,
                                      PointInputIterator end,
                                      PointMap point_map,
                                      Surface_mesh_on_cdt<GeomTraits>& output_mesh,
                                      const Top_view_surface_reconstruction_3::Parameters& parameters)
{
  typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
  typedef Border_graph<GeomTraits> Graph;
  typedef Line_detector<GeomTraits> Detect;

  CGAL_assertion (parameters.epsilon > 0.);
  
  double threshold_factor = 2.;
  double meshing_factor = std::sqrt(2.);
  
  SMCDT tmp_mesh;

  TOP_VIEW_CERR << "Inserting filtered points" << std::endl;
  Top_view_surface_reconstruction_3::internal::insert_filtered_points<GeomTraits>
    (begin, end, point_map, tmp_mesh, parameters.epsilon, parameters.quantile);

  TOP_VIEW_CERR << "Filtering faces" << std::endl;
  Top_view_surface_reconstruction_3::internal::filter_faces<GeomTraits>
    (tmp_mesh, parameters.epsilon * threshold_factor);

  TOP_VIEW_CERR << "Detecting planar regions" << std::endl;
  Top_view_surface_reconstruction_3::internal::region_growing<GeomTraits>
    (tmp_mesh, parameters.epsilon, 20, 0.9, false);

#ifdef TOP_VIEW_LOG
  tmp_mesh.DEBUG_dump_ply("debug1.ply");
#endif
  
  TOP_VIEW_CERR << "Cleaning buffer zone" << std::endl;
  Top_view_surface_reconstruction_3::internal::cleanup_buffer_zone<GeomTraits>
    (tmp_mesh);


  Graph graph;

  TOP_VIEW_CERR << "Building border graph" << std::endl;
  Top_view_surface_reconstruction_3::internal::build_border_graph<GeomTraits>
    (tmp_mesh, graph, 10. * parameters.epsilon * threshold_factor);

#ifdef TOP_VIEW_LOG
  graph.DEBUG_dump_poly("poly_complex.polylines.txt");
#endif
  
  TOP_VIEW_CERR << "Simplifying border graph" << std::endl;
  Top_view_surface_reconstruction_3::internal::simplify_border_graph<GeomTraits>
    (tmp_mesh, graph, parameters.epsilon * threshold_factor);

#ifdef TOP_VIEW_LOG
  tmp_mesh.DEBUG_dump_ply("debug2.ply");
  tmp_mesh.DEBUG_dump_poly();
  graph.DEBUG_dump_poly("poly_simplified.polylines.txt");
#endif


  TOP_VIEW_CERR << "Detecting lines in buffer" << std::endl;
  Detect detect (tmp_mesh);

  detect.sort_candidates_with_graph (graph);
  detect.run (parameters.epsilon /* * meshing_factor*/);

#ifdef TOP_VIEW_LOG
  detect.DEBUG_dump_ply("buffer_lines.ply");
  detect.DEBUG_dump_polyline("detected.polylines.txt");
  tmp_mesh.DEBUG_dump_off_2();
#endif

  TOP_VIEW_CERR << "Creating 3D regions" << std::endl;
  Top_view_surface_reconstruction_3::internal::create_3d_regions<GeomTraits>
    (tmp_mesh, detect, output_mesh, parameters.epsilon * threshold_factor);

#ifdef TOP_VIEW_LOG
  output_mesh.DEBUG_dump_off_0();
  output_mesh.DEBUG_dump_ply("debug3.ply");
#endif

}


  
} // namespace CGAL

#endif // CGAL_TOP_VIEW_SURFACE_RECONSTRUCTION_H
