#ifndef CGAL_TOP_VIEW_SURFACE_RECONSTRUCTION_H
#define CGAL_TOP_VIEW_SURFACE_RECONSTRUCTION_H

//#define TOP_VIEW_FIX_DUPLICATE_VERTICES
#define TOP_VIEW_DEBUG

#ifdef TOP_VIEW_DEBUG
#define TOP_VIEW_CERR std::cerr
#else
#define TOP_VIEW_CERR std::ostream(0)
#endif

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
  void filter_closed_polygons_with_no_3D_point_inside (Surface_mesh_on_cdt<GeomTraits>& mesh)
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

        if (mesh.has_at_least_one_mesh_vertex (current))
          okay = true;

        current->info() = typename SMCDT::Face_index(0); // Make non-default
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
        mesh.make_buffer (it);
    
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
              mesh.make_buffer(circ);
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
        
        current->info() = typename SMCDT::Face_index(0);

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
          mesh.add_face (done[i]);
      else
        for (std::size_t i = 0; i < done.size(); ++ i)
          mesh.make_buffer (done[i]);
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
          && ((mesh.is_buffer(f0) && !(mesh.is_buffer(f1))) ||
              (!(mesh.is_buffer(f0)) && mesh.is_buffer(f1))))
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
        if (mesh.has_mesh_face(circ) || mesh.is_ignored(circ))
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

    // Mark new faces as buffer
    for (typename SMCDT::Finite_faces_iterator it = mesh.finite_faces_begin ();
         it != mesh.finite_faces_end(); ++ it)
      if (mesh.is_default(it))
        mesh.make_buffer(it);
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
      if (mesh.is_ignored (it->first) ||
          mesh.is_ignored (it->first->neighbor(it->second)))
        continue;

      if (mesh.is_buffer(it->first) && mesh.is_buffer(it->first->neighbor(it->second)))
      {
        typename Map_v2v::iterator map0, map1;
        bool inserted = false;
        boost::tie (map0, inserted) = map_v2v.insert (std::make_pair (it->first, typename Graph::vertex_descriptor()));
        if (inserted)
        {
          typename GeomTraits::Point_2 p0 = mesh.midpoint (it->first);
          map0->second = boost::add_vertex(p0, graph);
        }

        boost::tie (map1, inserted) = map_v2v.insert (std::make_pair (it->first->neighbor(it->second),
                                                                      typename Graph::vertex_descriptor()));
        if (inserted)
        {
          typename GeomTraits::Point_2 p1 = mesh.midpoint (it->first->neighbor(it->second));
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
        mesh.make_handled_buffer (hint);
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
                                 Border_graph<GeomTraits>& graph,
                                 Surface_mesh_on_cdt<GeomTraits>& output,
                                 double epsilon)
  {
    typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;

    typename SMCDT::CDT map;

    // Insert borders in CDT and map
    for (std::size_t i = 0; i < graph.size(); ++ i)
    {
      typename SMCDT::Vertex_handle previous_cdt = output.insert (graph.point(graph[i][0]));
      typename SMCDT::Vertex_handle previous_map = map.insert (graph.point(graph[i][0]));
      
      for (std::size_t j = 1; j < graph[i].size(); ++ j)
      {
        const typename GeomTraits::Point_2& next = graph.point(graph[i][j]);
        typename GeomTraits::Vector_2 vec (previous_cdt->point(), next);
        std::size_t nb_pts = std::max(std::size_t(1), std::size_t(std::sqrt(vec*vec) / epsilon));

        typename SMCDT::Vertex_handle current_cdt = previous_cdt;

        for (std::size_t k = 1; k <= nb_pts; ++ k)
        {
          typename GeomTraits::Point_2 point
            = previous_cdt->point() + (k / double(nb_pts)) * vec;

          typename SMCDT::Vertex_handle next_cdt = output.insert (point);
          output.insert_constraint (current_cdt, next_cdt);
          current_cdt = next_cdt;
        }
        previous_cdt = current_cdt;

        typename SMCDT::Vertex_handle current_map = map.insert (graph.point(graph[i][j]));
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

    filter_closed_polygons_with_no_3D_point_inside<GeomTraits> (output);

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

} // namespace Top_view_surface_reconstruction_3
  
template <typename GeomTraits, typename PointInputIterator, typename PointMap>
void top_view_surface_reconstruction (PointInputIterator begin,
                                      PointInputIterator end,
                                      PointMap point_map,
                                      Surface_mesh_on_cdt<GeomTraits>& output_mesh,
                                      double spacing,
                                      double meshing_factor = 2.,
                                      double quantile = 0.5)
{
//  typedef typename GeomTraits::Point_3 Point_3;
  typedef Surface_mesh_on_cdt<GeomTraits> SMCDT;
  typedef Border_graph<GeomTraits> Graph;
  typedef Line_detector<GeomTraits> Detect;

  SMCDT tmp_mesh;

  TOP_VIEW_CERR << "Inserting filtered points" << std::endl;
  Top_view_surface_reconstruction_3::internal::insert_filtered_points<GeomTraits>
    (begin, end, point_map, tmp_mesh, spacing, quantile);

  TOP_VIEW_CERR << "Filtering faces" << std::endl;
  Top_view_surface_reconstruction_3::internal::filter_faces<GeomTraits>
    (tmp_mesh, spacing * meshing_factor);

  tmp_mesh.DEBUG_dump_off_0();

  TOP_VIEW_CERR << "Cleaning buffer zone" << std::endl;
  Top_view_surface_reconstruction_3::internal::cleanup_buffer_zone<GeomTraits>
    (tmp_mesh);

  tmp_mesh.DEBUG_dump_off_3();

  Graph graph;

  TOP_VIEW_CERR << "Building border graph" << std::endl;
  Top_view_surface_reconstruction_3::internal::build_border_graph<GeomTraits>
    (tmp_mesh, graph);

  graph.filter_small_terminal_borders(10. * spacing * meshing_factor);
  graph.split_into_polylines();
  
  tmp_mesh.DEBUG_dump_poly();
  graph.DEBUG_dump_poly("poly_complex.polylines.txt");

  
  TOP_VIEW_CERR << "Simplifying border graph" << std::endl;
  Top_view_surface_reconstruction_3::internal::simplify_border_graph<GeomTraits>
    (tmp_mesh, graph, spacing * meshing_factor);
  graph.DEBUG_dump_poly("poly_simplified.polylines.txt");

  TOP_VIEW_CERR << "Detecting lines in buffer" << std::endl;
  Detect detect (tmp_mesh);

  detect.sort_candidates_with_graph (graph);
  detect.run (spacing /* * meshing_factor*/);

  detect.DEBUG_dump_ply("buffer_lines.ply");
  detect.DEBUG_dump_polyline("detected.polylines.txt");

  tmp_mesh.DEBUG_dump_off_2();


  TOP_VIEW_CERR << "Creating mesh with borders" << std::endl;
  Top_view_surface_reconstruction_3::internal::create_mesh_with_borders<GeomTraits>
    (tmp_mesh, graph, output_mesh, spacing * meshing_factor);

//  output_mesh.DEBUG_dump_off_0();
  output_mesh.DEBUG_dump_off_4();

  output_mesh.DEBUG_dump_off_1();

  output_mesh.check_structure_integrity();
  TOP_VIEW_CERR << "Generating missing 3D points" << std::endl;
  Top_view_surface_reconstruction_3::internal::generate_missing_3d_points<GeomTraits>
    (output_mesh, spacing * meshing_factor);


  TOP_VIEW_CERR << "Snapping intersecting borders" << std::endl;
  Top_view_surface_reconstruction_3::internal::snap_intersecting_borders<GeomTraits>
    (output_mesh);
  output_mesh.check_structure_integrity();

  TOP_VIEW_CERR << "Generating vertical walls" << std::endl;
  Top_view_surface_reconstruction_3::internal::generate_vertical_walls<GeomTraits>
    (output_mesh, spacing * meshing_factor);

#ifndef TOP_VIEW_FIX_DUPLICATE_VERTICES
  output_mesh.stitch_borders();
#endif

  output_mesh.check_structure_integrity();

  output_mesh.DEBUG_dump_off_5();

}

} // namespace CGAL

#endif // CGAL_TOP_VIEW_SURFACE_RECONSTRUCTION_H
