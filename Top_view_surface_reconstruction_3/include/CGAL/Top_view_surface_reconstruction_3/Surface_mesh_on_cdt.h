#ifndef CGAL_TVSR_SURFACE_MESH_ON_CDT
#define CGAL_TVSR_SURFACE_MESH_ON_CDT

#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/Surface_mesh.h>

#include <fstream>
#include <queue>

namespace CGAL
{

template <typename GeomTraits>
class Surface_mesh_on_cdt
{
public:

  typedef Surface_mesh_on_cdt<GeomTraits> Self;
  
  typedef GeomTraits Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Direction_2 Direction_2;
  typedef Surface_mesh<Point_3> Mesh;
  typedef typename Mesh::Vertex_index Vertex_index;
  typedef typename Mesh::Face_index Face_index;

  typedef Triangulation_vertex_base_with_info_2<std::vector<std::pair<Direction_2, Vertex_index> >, Kernel>  Vbwi;
  typedef Triangulation_face_base_with_info_2<Face_index, Kernel> Fbwi;
  typedef Constrained_triangulation_face_base_2<Kernel, Fbwi> Cfbwi;
  typedef Triangulation_data_structure_2<Vbwi, Cfbwi> TDS;
  typedef Exact_predicates_tag Itag;
  typedef Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> CDT;

  typedef typename CDT::Face_handle Face_handle;
  typedef typename CDT::Vertex_handle Vertex_handle;
  typedef typename CDT::Edge Edge;
  typedef typename CDT::Finite_edges_iterator Finite_edges_iterator;
  typedef typename CDT::Finite_faces_iterator Finite_faces_iterator;
  typedef typename CDT::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename CDT::Edge_circulator Edge_circulator;
  typedef typename CDT::Vertex_circulator Vertex_circulator;
  typedef typename CDT::Face_circulator Face_circulator;
  typedef typename CDT::Locate_type Locate_type;

  typedef typename Mesh::template Property_map<Vertex_index, Vertex_handle> V2V_map;
  typedef typename Mesh::template Property_map<Face_index, Face_handle> F2F_map;

private:

  Mesh m_mesh;
  CDT m_cdt;

  V2V_map m_v2v_map;
  F2F_map m_f2f_map;

public:

  Surface_mesh_on_cdt()
  {
    bool okay = false;
    boost::tie(m_v2v_map, okay) = m_mesh.template add_property_map<Vertex_index, Vertex_handle>("v:cdt_vertex");
    assert (okay);
    boost::tie(m_f2f_map, okay) = m_mesh.template add_property_map<Face_index, Face_handle>("f:cdt_face");
    assert (okay);
  }

  Vertex_handle cdt_vertex (Vertex_index vi) const { return m_v2v_map[vi]; }
  Vertex_index mesh_vertex (Vertex_handle vh, std::size_t idx = 0) const { return vh->info()[idx].second; }
  bool has_mesh_vertex (Vertex_handle vh) const { return !vh->info().empty(); }
  
  bool has_unique_mesh_vertex (Vertex_handle vh) const
  {
    return (vh->info().size() == 1 && vh->info()[0].first == Direction_2(0,1));
  }
  bool is_border_vertex (Vertex_index vi) const { return (degree(cdt_vertex(vi)) > 1); }

  bool has_defined_height (Vertex_index vi) const
  {
    const Point_3& p = point(vi);
    return (p.z() == p.z()); // NaN test
  }


  Face_handle cdt_face (Face_index fi) const { return m_f2f_map[fi]; }
  Face_index mesh_face (Face_handle fh) const { return fh->info(); }

  Face_handle locate (const Point_2& point, Face_handle hint = Face_handle()) const
  { return m_cdt.locate (point, hint); }
  typename GeomTraits::Triangle_2 triangle (Face_handle fh) const { return m_cdt.triangle(fh); }

  Vertex_handle infinite_vertex() { return m_cdt.infinite_vertex(); }
  Finite_faces_iterator finite_faces_begin() { return m_cdt.finite_faces_begin(); }
  Finite_faces_iterator finite_faces_end() { return m_cdt.finite_faces_end(); }
  Finite_vertices_iterator finite_vertices_begin() { return m_cdt.finite_vertices_begin(); }
  Finite_vertices_iterator finite_vertices_end() { return m_cdt.finite_vertices_end(); }
  Finite_edges_iterator finite_edges_begin() { return m_cdt.finite_edges_begin(); }
  Finite_edges_iterator finite_edges_end() { return m_cdt.finite_edges_end(); }
  Face_circulator incident_faces (Vertex_handle vh) { return m_cdt.incident_faces (vh); }
  Edge_circulator incident_edges (Vertex_handle vh) { return m_cdt.incident_edges (vh); }
  Vertex_circulator incident_vertices (Vertex_handle vh) { return m_cdt.incident_vertices (vh); }
  void insert_constraint (Vertex_handle v0, Vertex_handle v1) { m_cdt.insert_constraint (v0, v1); }
  bool are_there_incident_constraints (Vertex_handle vh) { return m_cdt.are_there_incident_constraints (vh); }
  bool is_constrained (Edge e) { return m_cdt.is_constrained (e); }
                                  
  bool is_infinite (Face_handle fh) const { return m_cdt.is_infinite(fh); }

  bool has_mesh_face (Face_handle fh) const { return (int(fh->info()) >= 0); }
  bool is_default (Face_handle fh) const { return (int(fh->info()) == -1); }
  void make_default (Face_handle fh) { fh->info() = Face_index(-1); }
  bool is_buffer (Face_handle fh) const { return (int(fh->info()) == -2); }
  void make_buffer (Face_handle fh) { fh->info() = Face_index(-2); }
  bool is_ignored (Face_handle fh) const { return (int(fh->info()) == -3); }
  void make_ignored (Face_handle fh) { fh->info() = Face_index(-3); }

  int cw (int i) const { return m_cdt.cw(i); }
  int ccw (int i) const { return m_cdt.ccw(i); }

  const Point_3& point (Vertex_index vi) const { return m_mesh.point(vi); }
  const Point_3& point (Vertex_handle vh) const { return m_mesh.point(mesh_vertex(vh)); }

  const Point_3& point_from_vertex_view (Vertex_handle vh, Vertex_handle view) const
  {
    // TODO
  }
  
  Vertex_handle insert (const Point_2& point)
  {
    return m_cdt.insert (point);
  }
  
  Vertex_handle insert (const Point_3& point)
  {
    Vertex_handle vh = m_cdt.insert (Point_2 (point.x(), point.y()));
    Vertex_index vi = m_mesh.add_vertex (point);
    vh->info().push_back(std::make_pair(Direction_2(0,1), vi));
    m_v2v_map[vi] = vh;
    return vh;
  }

  Vertex_handle insert (Vertex_handle vh, const Point_3& point, const Direction_2& direction)
  {
    Vertex_index vi = m_mesh.add_vertex (point);
    vh->info().push_back(std::make_pair(direction, vi));
    m_v2v_map[vi] = vh;
    return vh;
  }

  void remove (Vertex_handle vh)
  {
    m_cdt.remove (vh);
  }

  void remove_incident_constraints (Vertex_handle vh)
  {
    m_cdt.remove_incident_constraints (vh);
  }

  std::size_t degree (Vertex_handle vh)
  {
    std::size_t out = 0;
    Edge_circulator circ = m_cdt.incident_edges(vh), start = circ;
    do
    {
      if (m_cdt.is_constrained(*circ))
        ++ degree;
      ++ circ;
    }
    while (circ != start);
    return degree;
  }

  void add_face (Face_handle fh)
  {
    Face_index fi = m_mesh.add_face (fh->vertex(0)->info()[0].second,
                                     fh->vertex(1)->info()[0].second,
                                     fh->vertex(2)->info()[0].second);
    fh->info() = fi;
    m_f2f_map[fi] = fh;
  }

  void add_face (Face_handle fh, Vertex_index a, Vertex_index b, Vertex_index c)
  {
    Face_index fi = m_mesh.add_face (a, b, c);
    fh->info() = fi;
    m_f2f_map[fi] = fh;
  }

  std::size_t find_section_of_point_from_vertex_view (Vertex_handle vh, const Point_2& point)
  {
    if (vh->info().size() == 1)
      return 0;

    Direction_2 direction (Vector_2 (vh->point(), point));
    for (std::size_t i = 0; i <= vh->info().size(); ++ i)
      if (direction.counterclockwise_in_between (vh->info()[i].first,
                                                 vh->info()[(i+1)%(vh->info().size())].first))
        return i;
    
    return 0;
  }

  void get_neighborhood (Vertex_handle vh, double epsilon,
                         std::vector<std::vector<Vertex_handle> >& neighborhood)
  {
    neighborhood.resize (vh->info().size());
    
    std::queue<Vertex_handle> todo;
    Vertex_circulator circ = m_cdt.incident_vertices (vh), start = circ;
    do
    {
      todo.push (circ);
      ++ circ;
    }
    while (circ != start);
    
    std::set<Vertex_handle> done;
    while (!todo.empty())
    {
      Vertex_handle current = todo.front();
      todo.pop();

      if (!(done.insert(current).second) ||
          m_cdt.is_infinite(current) ||
          !has_unique_mesh_vertex(current) ||
          CGAL::squared_distance (vh->point(), current->point()) > epsilon * epsilon)
        continue;

      std::size_t i = find_section_of_point_from_vertex_view (vh, current->point());
      neighborhood[i].push_back (current);
      
      Vertex_circulator circ = m_cdt.incident_vertices (current), start = circ;
      do
      {
        todo.push (circ);
        ++ circ;
      }
      while (circ != start);
    }
  }

  double estimate_height_with_pca (Vertex_handle vh, std::vector<Vertex_handle>& neighborhood)
  {
    // TODO
    return 0.;
  }

  void estimate_missing_heights (Vertex_handle vh, double epsilon)
  {
    std::ofstream f ("pts.xyz", std::ios_base::app);
    std::vector<std::vector<Vertex_handle> > small_neighborhood;
    get_neighborhood (vh, 3. * epsilon, small_neighborhood);

    std::vector<std::vector<Vertex_handle> > big_neighborhood;
    for (std::size_t i = 0; i < small_neighborhood.size(); ++ i)
      if (small_neighborhood[i].empty())
      {
        get_neighborhood (vh, 15. * epsilon, big_neighborhood);
        break;
      }
    
    for (std::size_t i = 0; i < vh->info().size(); ++ i)
    {
      // Get small neighborhood
      std::vector<Vertex_handle>* neighborhood = &(small_neighborhood[i]);
      
      // If not large enough, get big neighborhood
      if (neighborhood->empty())
        neighborhood = &(big_neighborhood[i]);

      // If not large enough, give up for now
      if (neighborhood->empty())
        continue;

      // Reference height = closest neighbor's height
      Vertex_handle closest = Vertex_handle();
      double dist_min = std::numeric_limits<double>::max();
      for (std::size_t j = 0; j < neighborhood->size(); ++ j)
      {
        double dist = CGAL::squared_distance (vh->point(), (*neighborhood)[j]->point());
        if (dist < dist_min)
        {
          dist_min = dist;
          closest = (*neighborhood)[j];
        }
      }


      Point_3 new_point = m_mesh.point(vh->info()[i].second);

      double ref_z = point(closest).z();

      // Try PCA if not too far from reference
      double pca_z = estimate_height_with_pca (vh, *neighborhood);

      if (CGAL::abs(ref_z - pca_z) < epsilon)
        new_point = Point_3 (new_point.x(), new_point.y(), pca_z);
      else
        new_point = Point_3 (new_point.x(), new_point.y(), ref_z);

      m_mesh.point(vh->info()[i].second) = new_point;
    }

  }

  void DEBUG_dump_off_1() 
  {
    for (typename Mesh::Vertex_range::iterator it = m_mesh.vertices().begin();
         it != m_mesh.vertices().end(); ++ it)
    {
      const Point_3& p = m_mesh.point(*it);
      if (p.z() != p.z())
        m_mesh.point(*it) = Point_3 (p.x(), p.y(), 0);
    }
      
    
    std::ofstream f("debug1.off");
    f.precision(18);
    f << m_mesh;    
  }
  
  void DEBUG_dump_off_0() const
  {
    std::ofstream f("debug0.off");
    f.precision(18);
    
    std::size_t nb_faces = 0;
    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
      if (has_mesh_face(it))
        nb_faces ++;

    f << "OFF" << std::endl << m_cdt.number_of_vertices() << " " << nb_faces << " 0" << std::endl;

    std::map<Vertex_handle, std::size_t> map;
    std::size_t idx = 0;
    for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin(); it != m_cdt.finite_vertices_end(); ++ it)
    {
      f << it->point() << " 0" << std::endl;
      map[it] = idx ++;
    }

    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
      if (has_mesh_face(it))
      {
        f << "3";
        for (std::size_t i = 0; i < 3; ++ i)
          f << " " << map[it->vertex(i)] << std::endl;
        f << std::endl;
      }
  }

  void DEBUG_dump_off_2() const
  {
    std::ofstream f("debug2.off");
    f.precision(18);
    
    std::size_t nb_faces = 0;
    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
      if (!(has_mesh_face(it)))
        nb_faces ++;

    f << "OFF" << std::endl << m_cdt.number_of_vertices() << " " << nb_faces << " 0" << std::endl;

    std::map<Vertex_handle, std::size_t> map;
    std::size_t idx = 0;
    for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin(); it != m_cdt.finite_vertices_end(); ++ it)
    {
      f << it->point() << " 0" << std::endl;
      map[it] = idx ++;
    }

    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
      if (!(has_mesh_face(it)))
      {
        f << "3";
        for (std::size_t i = 0; i < 3; ++ i)
          f << " " << map[it->vertex(i)] << std::endl;
        f << std::endl;
      }
  }


  void DEBUG_dump_off_3() const
  {
    std::ofstream f("debug3.off");
    f.precision(18);
    
    std::size_t nb_faces = 0;
    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
      if (is_ignored(it))
        nb_faces ++;

    f << "OFF" << std::endl << m_cdt.number_of_vertices() << " " << nb_faces << " 0" << std::endl;

    std::map<Vertex_handle, std::size_t> map;
    std::size_t idx = 0;
    for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin(); it != m_cdt.finite_vertices_end(); ++ it)
    {
      f << it->point() << " 0" << std::endl;
      map[it] = idx ++;
    }

    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
      if (is_ignored(it))
      {
        f << "3";
        for (std::size_t i = 0; i < 3; ++ i)
          f << " " << map[it->vertex(i)] << std::endl;
        f << std::endl;
      }
  }

  void DEBUG_dump_off_4() const
  {
    std::ofstream f("debug4.off");
    f.precision(18);
    
    std::size_t nb_faces = 0;
    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
      nb_faces ++;

    f << "OFF" << std::endl << m_cdt.number_of_vertices() << " " << nb_faces << " 0" << std::endl;

    std::map<Vertex_handle, std::size_t> map;
    std::size_t idx = 0;
    for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin(); it != m_cdt.finite_vertices_end(); ++ it)
    {
      f << it->point() << " 0" << std::endl;
      map[it] = idx ++;
    }

    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
    {
      f << "3";
      for (std::size_t i = 0; i < 3; ++ i)
        f << " " << map[it->vertex(i)] << std::endl;
      f << std::endl;
    }
  }

  void DEBUG_dump_poly()
  {
    std::ofstream f("debug_constraints.polylines.txt");
    f.precision(18);

    for (Finite_edges_iterator it = m_cdt.finite_edges_begin(); it != m_cdt.finite_edges_end(); ++ it)
      if (m_cdt.is_constrained (*it))
        f << "2 "
          << it->first->vertex ((it->second + 1)%3)->point() << " 0 "
          << it->first->vertex ((it->second + 2)%3)->point() << " 0" << std::endl;
   }

};


}
  
#endif // CGAL_TVSR_SURFACE_MESH_ON_CDT
