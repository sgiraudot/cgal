#ifndef CGAL_TVSR_SURFACE_MESH_ON_CDT
#define CGAL_TVSR_SURFACE_MESH_ON_CDT

#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Random.h>

#include <CGAL/unordered.h>

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
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Line_3 Line_3;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Direction_2 Direction_2;
  typedef Surface_mesh<Point_3> Mesh;
  typedef typename Mesh::Vertex_index Vertex_index;
  typedef typename Mesh::Face_index Face_index;
  typedef typename Mesh::Halfedge_index Halfedge_index;
  typedef typename Mesh::Edge_index Edge_index;

  struct Face_info
  {
    Face_index index;
    Point_2 endpoint;
    std::vector<std::size_t> incident_lines;
    std::size_t plane_index;
    
    Face_info()
      : index()
      , endpoint (std::numeric_limits<double>::quiet_NaN(),
                  std::numeric_limits<double>::quiet_NaN())
      , plane_index (std::size_t(-1))
    { }

    bool has_endpoint() const { return (endpoint.x() == endpoint.x()); }
    void erase_endpoint() { endpoint = Point_2(std::numeric_limits<double>::quiet_NaN(),
                                               std::numeric_limits<double>::quiet_NaN()); }

    bool has_plane() const { return plane_index != std::size_t(-1); }
    void erase_plane() { plane_index = std::size_t(-1); }
  };

  typedef Triangulation_vertex_base_with_info_2<std::vector<std::pair<Direction_2, Vertex_index> >, Kernel>  Vbwi;
  typedef Triangulation_face_base_with_info_2<Face_info, Kernel> Fbwi;
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
  typedef typename CDT::Line_face_circulator Line_face_circulator;
  typedef typename CDT::Locate_type Locate_type;

  typedef typename Mesh::template Property_map<Vertex_index, Vertex_handle> V2V_map; // Mesh vertex to CDT vertex
  typedef typename Mesh::template Property_map<Vertex_index, Vertex_index> NVN_map; // Mesh vertex to next vertical mesh vertex
  typedef typename Mesh::template Property_map<Face_index, Face_handle> F2F_map; // Mesh face to CDT face

  struct Sort_faces_by_planarity
  {
    Self* mesh;

    Sort_faces_by_planarity (Self* mesh) : mesh (mesh) { }

    bool operator() (const Face_handle& a, const Face_handle& b) const
    {
      double dev_a = deviation_from_plane (a);
      double dev_b = deviation_from_plane (b);
      if (dev_a == dev_b)
        return a < b;
      return deviation_from_plane (a) < deviation_from_plane (b);
    }

    double deviation_from_plane (const Face_handle& fh) const
    {
      Vector_3 normal_seed = mesh->normal_vector (fh);
      Point_3 pt_seed = mesh->midpoint_3(fh);
      Plane_3 optimal_plane (pt_seed, normal_seed);

      double max = 0.;
      std::size_t nb_neigh = 0;
      
      for (std::size_t i = 0; i < 3; ++ i)
      {
        Face_handle fn = fh->neighbor(i);
        if (mesh->has_mesh_face(fn))
        {
          nb_neigh ++;
          Point_3 point = mesh->point (fn->vertex(fn->index(fh)));
          double dist = CGAL::squared_distance (point, optimal_plane);
          if (dist > max)
            max = dist;
        }
      }

      if (nb_neigh == 0)
        return std::numeric_limits<double>::max();
      return max * (4 - nb_neigh);
    }

  };
    

private:

  struct Vertex_handle_to_point : public std::unary_function<Vertex_handle, Point_3>
  {
    Self& mesh;

    Vertex_handle_to_point (Self& mesh) : mesh (mesh) { }

    Point_3 operator()(const Vertex_handle& vh) const
    {
      return mesh.point(vh);
    }
    
  };

  Mesh m_mesh;
  CDT m_cdt;

  V2V_map m_v2v_map;
  NVN_map m_nvv_map;
  F2F_map m_f2f_map;
  std::vector<Plane_3> m_planes;
  std::vector<double> m_plane_weights;

public:

  Surface_mesh_on_cdt()
  {
    bool okay = false;
    boost::tie(m_v2v_map, okay) = m_mesh.template add_property_map<Vertex_index, Vertex_handle>("v:cdt_vertex");
    assert (okay);
    boost::tie(m_nvv_map, okay) = m_mesh.template add_property_map<Vertex_index, Vertex_index>("v:next_vertical_vertex");
    assert (okay);
    boost::tie(m_f2f_map, okay) = m_mesh.template add_property_map<Face_index, Face_handle>("f:cdt_face");
    assert (okay);
  }

  const Mesh& mesh() const { return m_mesh; }
  const CDT& cdt() const { return m_cdt; }
  std::vector<Plane_3>& planes() { return m_planes; }
  std::vector<double>& plane_weights() { return m_plane_weights; }
  
  Vertex_handle cdt_vertex (Vertex_index vi) const { return m_v2v_map[vi]; }
  Vertex_index mesh_vertex (Vertex_handle vh, std::size_t idx = 0) const { return vh->info()[idx].second; }
  bool has_mesh_vertex (Vertex_handle vh) const { return !vh->info().empty(); }
  std::size_t number_of_mesh_vertices (Vertex_handle vh) const { return vh->info().size(); }
  
  bool has_unique_mesh_vertex (Vertex_handle vh) const
  {
    return (vh->info().size() == 1 && vh->info()[0].first == Direction_2(0,0));
  }
  bool is_border_vertex (Vertex_index vi) const { return (degree(cdt_vertex(vi)) > 1); }

  bool has_defined_height (Vertex_index vi) const
  {
    const Point_3& p = point(vi);
    return (p.z() == p.z()); // NaN test
  }

  bool is_valid (Vertex_handle vh) const
  {
    for (std::size_t i = 0; i < vh->info().size(); ++ i)
      if (point (vh, i).z() != point (vh, i).z())
        return false;
    return true;
  }

  Vertex_index next_vertex (Vertex_index vi)
  {
    return m_nvv_map[vi];
  }
  void set_next_vertex (Vertex_index vi, Vertex_index vnext)
  {
    m_nvv_map[vi] = vnext;
  }


  Face_handle cdt_face (Face_index fi) const { return m_f2f_map[fi]; }
  Face_index mesh_face (Face_handle fh) const { return fh->info().index; }

  Face_handle locate (const Point_2& point, Face_handle hint = Face_handle()) const
  { return m_cdt.locate (point, hint); }

  Edge mirror_edge (const Edge& edge) const { return m_cdt.mirror_edge(edge); }
  
  typename GeomTraits::Triangle_2 triangle (Face_handle fh) const { return m_cdt.triangle(fh); }

  typename GeomTraits::Triangle_3 triangle_3 (Face_handle fh)
  {
    CGAL_assertion (has_unique_mesh_vertex(fh->vertex(0))
                    && has_unique_mesh_vertex(fh->vertex(1))
                    && has_unique_mesh_vertex(fh->vertex(2)));

    return typename GeomTraits::Triangle_3 (point(fh->vertex(0)),
                                            point (fh->vertex(1)),
                                            point (fh->vertex(2)));
  }

  bool is_edge (Vertex_handle va, Vertex_handle vb, Face_handle& fh, int& fi)
  {
    return m_cdt.is_edge (va, vb, fh, fi);
  }
  
  Line_face_circulator line_walk (const Point_2& p, const Point_2& q, Face_handle f = Face_handle())
  {
    return m_cdt.line_walk(p,q,f);
  }

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
  void remove_constraint (const Edge& e) { m_cdt.remove_constraint (e.first, e.second); }
  bool are_there_incident_constraints (Vertex_handle vh) { return m_cdt.are_there_incident_constraints (vh); }

  template <typename OutputItEdges>
  OutputItEdges incident_constraints (Vertex_handle vh, OutputItEdges output) const
  {
    return m_cdt.incident_constraints (vh, output);
  }

  template <typename OutputItFaces>
  OutputItFaces get_conflicts (const Point_2& p, OutputItFaces fit, Face_handle start = Face_handle())
  {
    return m_cdt.get_conflicts (p, fit, start);
  }

  bool remove_constraint (Vertex_handle a, Vertex_handle b)
  {
    Face_handle f;
    int idx;
    if(m_cdt.is_edge (a, b, f, idx))
    {
      m_cdt.remove_constraint (f, idx);
      return true;
    }
    
    return false;
  }
  
  bool is_constrained (Edge e) { return m_cdt.is_constrained (e); }
                                  
  bool is_infinite (Face_handle fh) const { return m_cdt.is_infinite(fh); }

  bool has_at_least_one_mesh_vertex (Face_handle f) const
  {
    for (std::size_t i = 0; i < 3; ++ i)
      if (has_mesh_vertex (f->vertex(i)))
        return true;
    return false;
  }
  
  bool has_three_mesh_vertices (Face_handle f) const
  {
    for (std::size_t i = 0; i < 3; ++ i)
      if (!has_mesh_vertex (f->vertex(i)))
        return false;
    return true;
  }

  bool has_cdt_face (Face_index fi) const { return m_f2f_map[fi] != Face_handle(); }
  bool has_mesh_face (Face_handle fh) const { return (int(fh->info().index) >= 0); }
  bool is_default (Face_handle fh) const { return (int(fh->info().index) == -1); }
  void make_default (Face_handle fh) { fh->info().index = Face_index(-1); }
  bool is_buffer (Face_handle fh) const
  { return (int(fh->info().index) == -2
            || int(fh->info().index) == -4
            || int(fh->info().index) == -5);
  }
  void make_buffer (Face_handle fh) { fh->info().index = Face_index(-2); }
  bool is_ignored (Face_handle fh) const { return (int(fh->info().index) == -3); }
  void make_ignored (Face_handle fh) { fh->info().index = Face_index(-3); }
  bool is_default_buffer (Face_handle fh) const { return (int(fh->info().index) == -2); }
  bool is_wall_buffer (Face_handle fh) const { return (int(fh->info().index) == -4); }
  void make_wall_buffer (Face_handle fh) { fh->info().index = Face_index(-4); }
  bool is_ridge_buffer (Face_handle fh) const { return (int(fh->info().index) == -5); }
  void make_ridge_buffer (Face_handle fh) { fh->info().index = Face_index(-5); }
  bool is_ignored_buffer (Face_handle fh) const { return (int(fh->info().index) == -7); }
  void make_ignored_buffer (Face_handle fh) { fh->info().index = Face_index(-7); }
  bool is_pending (Face_handle fh) const { return (int(fh->info().index) == -6); }
  void make_pending (Face_handle fh) { fh->info().index = Face_index(-6); }

  int cw (int i) const { return m_cdt.cw(i); }
  int ccw (int i) const { return m_cdt.ccw(i); }

  const Point_3& point (Vertex_index vi) const { return m_mesh.point(vi); }
  Point_3& point (Vertex_index vi) { return m_mesh.point(vi); }
  const Point_3& point (Vertex_handle vh, std::size_t idx = 0) const { return m_mesh.point(mesh_vertex(vh, idx)); }
  Point_3& point (Vertex_handle vh, std::size_t idx = 0) { return m_mesh.point(mesh_vertex(vh, idx)); }

  Vector_3 normal_vector (Face_handle fh) const
  {
    CGAL_assertion (has_unique_mesh_vertex(fh->vertex(0))
                    && has_unique_mesh_vertex(fh->vertex(1))
                    && has_unique_mesh_vertex(fh->vertex(2)));

    Vector_3 out = CGAL::orthogonal_vector (point (fh->vertex(0)),
                                            point (fh->vertex(1)),
                                            point (fh->vertex(2)));
    out = out / std::sqrt (out*out);
    return out;
  }
  
  Vertex_handle insert (const Point_2& point)
  {
    return m_cdt.insert (point);
  }
  
  Vertex_handle insert (const Point_3& point)
  {
    Vertex_handle vh = m_cdt.insert (Point_2 (point.x(), point.y()));
    Vertex_index vi = m_mesh.add_vertex (point);
    vh->info().push_back(std::make_pair(Direction_2(0,0), vi));
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

  Vertex_index insert (Vertex_handle vh, const Point_3& point)
  {
    CGAL_assertion (vh->info().empty());

    Vertex_index vi = m_mesh.add_vertex (point);
    vh->info().push_back(std::make_pair(Direction_2(0,0), vi));
    m_v2v_map[vi] = vh;
    return vi;
  }

  Vertex_index insert (Vertex_handle vh)
  {
    CGAL_assertion (vh->info().empty());

    Vertex_index vi = m_mesh.add_vertex (Point_3 (vh->point().x(), vh->point().y(),
                                                  std::numeric_limits<double>::quiet_NaN()));
    vh->info().push_back(std::make_pair(Direction_2(0,0), vi));
    m_v2v_map[vi] = vh;
    return vi;
  }

  void remove (Vertex_handle vh)
  {
    remove_mesh_vertex (vh);
    m_cdt.remove (vh);
  }

  void remove (Vertex_index vi)
  {
    Vertex_handle vh = m_v2v_map[vi];
    for (std::size_t i = 0; i < vh->info().size(); ++ i)
      if (vh->info()[i].second == vi)
      {
        for (std::size_t j = i; j < vh->info().size() - 1; ++ j)
          vh->info()[j] = vh->info()[j+1];
        vh->info().resize(vh->info().size() - 1);
        break;
      }
    m_mesh.remove_vertex (vi);
  }

  void remove_mesh_vertex (Vertex_handle vh)
  {
    if (vh->info().size() == 1)
    {
      bool faces = false;

      BOOST_FOREACH (Face_index fi, faces_around_target (halfedge (vh->info()[0].second, m_mesh), m_mesh))
      {
        if (!faces)
          std::cerr << "Remove vertex " << vh->info()[0].second << std::endl;
        faces = true;
        std::cerr << " * " << fi << std::endl;
      }
      if (faces)
        exit(0);
      
      m_mesh.remove_vertex (vh->info()[0].second);
      vh->info().clear();
    }
  }

  void remove_mesh_face (Face_handle fh)
  {
    m_mesh.remove_face(fh->info().index);
    fh->info().index = Face_index();
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
        ++ out;
      ++ circ;
    }
    while (circ != start);
    return out;
  }



  void add_face (Face_handle fh)
  {
    Face_index fi = m_mesh.add_face (fh->vertex(0)->info()[0].second,
                                     fh->vertex(1)->info()[0].second,
                                     fh->vertex(2)->info()[0].second);
    fh->info().index = fi;
    m_f2f_map[fi] = fh;
  }

  void add_face (Face_handle fh, Vertex_index a, Vertex_index b, Vertex_index c)
  {
    Face_index fi = m_mesh.add_face (a, b, c);
    if (fi != Face_index())
    {
      fh->info().index = fi;
      m_f2f_map[fi] = fh;
    }
  }

  bool add_face (Vertex_index a, Vertex_index b, Vertex_index c)
  {
    Face_index fi =
      m_mesh.add_face (a, b, c);

    return (fi != Face_index());
  }

  void merge_vertices (Vertex_index a, Vertex_index b)
  {
    Point_3 new_point (point(a).x(), point(a).y(),
                       0.5 * (point(a).z() + point(b).z()));

    point(a) = new_point;

#ifdef TOP_VIEW_FIX_DUPLICATE_VERTICES
    Halfedge_index h = m_mesh.halfedge (b);

    do
    {
      Halfedge_index current = h;
      h = m_mesh.next_around_target (h);
      m_mesh.set_target (current, a);
    }
    while (!m_mesh.is_border(h));

    Vertex_handle vh = cdt_vertex (b);

    for (std::size_t i = 0; i < number_of_mesh_vertices (vh); ++ i)
      if (mesh_vertex (vh, i) == b)
        vh->info()[i].second = a;
      
//    m_mesh.remove_vertex (b);
#else
    point(b) = new_point;

    return;
#endif
  }

  void stitch_borders()
  {
    CGAL::Polygon_mesh_processing::stitch_borders (m_mesh);
  }

  std::size_t find_section_of_point_from_vertex_view (Vertex_handle vh, const Point_2& point)
  {
    CGAL_assertion (!vh->info().empty());
    
    if (vh->info().size() == 1)
      return 0;

    Direction_2 direction (Vector_2 (vh->point(), point));
    for (std::size_t i = 0; i < vh->info().size(); ++ i)
      if (direction.counterclockwise_in_between (vh->info()[i].first,
                                                 vh->info()[(i+1)%(vh->info().size())].first))
        return i;

    return 0;
  }

  std::size_t find_section_of_point_from_vertex_view_DEBUG (Vertex_handle vh, const Point_2& point)
  {
    if (vh->info().size() == 1)
      return 0;

    static std::size_t occur = 0;
    static Point_2 to_log (72110.151041666672, 162618.42708333334);

    Direction_2 direction (Vector_2 (vh->point(), point));
    for (std::size_t i = 0; i < vh->info().size(); ++ i)
      if (direction.counterclockwise_in_between (vh->info()[i].first,
                                                 vh->info()[(i+1)%(vh->info().size())].first))
      {
        if (CGAL::squared_distance(vh->point(), to_log) < 1e-10)
        {
          occur ++;
          if (occur == 3)
          {
            std::ofstream f1("section.polylines.txt");
            f1.precision(18);
            f1 << "3 " << (vh->point() + vh->info()[i].first.to_vector()) << " 0 "
               << vh->point() << " 0 "
               << (vh->point() + vh->info()[(i+1)%(vh->info().size())].first.to_vector()) << " 0 " << std::endl;
            std::ofstream f2("section.xyz");
            f2.precision(18);
            f2 << point << " 0" << std::endl;
          }
        }
        return i;
      }
    std::cerr << "Warning section" << std::endl;
    return 0;
  }

  Point_2 midpoint (Face_handle f)
  {
    return Point_2
      ((f->vertex(0)->point().x() + f->vertex(1)->point().x() + f->vertex(2)->point().x()) / 3.,
       (f->vertex(0)->point().y() + f->vertex(1)->point().y() + f->vertex(2)->point().y()) / 3.);
  }

  Point_3 midpoint_3 (Face_handle fh)
  {
    CGAL_assertion (has_unique_mesh_vertex(fh->vertex(0))
                    && has_unique_mesh_vertex(fh->vertex(1))
                    && has_unique_mesh_vertex(fh->vertex(2)));

    Point_3 p0 = point(fh->vertex(0));
    Point_3 p1 = point(fh->vertex(1));
    Point_3 p2 = point(fh->vertex(2));

    return Point_3 ((p0.x() + p1.x() + p2.x()) / 3.,
                    (p0.y() + p1.y() + p2.y()) / 3.,
                    (p0.z() + p1.z() + p2.z()) / 3.);
  }


  void get_neighborhood (Vertex_handle vh, double epsilon,
                         std::vector<std::vector<Point_3> >& neighborhood,
                         bool allow_border_propagation)
  {
    neighborhood.resize (vh->info().size());
    
    std::queue<std::pair<Face_handle, std::size_t> > todo;
    Face_circulator circ = m_cdt.incident_faces (vh), start = circ;
    do
    {
      if (!is_ignored(circ))
      {
        Point_2 p = midpoint(circ);
        std::size_t i = find_section_of_point_from_vertex_view (vh, p);
        todo.push (std::make_pair (circ, i));
      }
      ++ circ;
    }
    while (circ != start);
    
    cpp11::unordered_set<Face_handle> done;
    cpp11::unordered_set<Point_3> vertex_done;
    
    while (!todo.empty())
    {
      Face_handle current = todo.front().first;
      std::size_t section = todo.front().second;
      todo.pop();

      if (!(done.insert(current).second) ||
          m_cdt.is_infinite(current) ||
          is_ignored(current))
        continue;

      bool still_in = false;

      for (std::size_t i = 0; i < 3; ++ i)
      {
        Vertex_handle vcurrent = current->vertex(i);

        if (CGAL::squared_distance (vh->point(), vcurrent->point()) < epsilon * epsilon)
          still_in = true;
        
        if (has_unique_mesh_vertex(vcurrent)
            && vertex_done.insert (point(vcurrent)).second)
          neighborhood[section].push_back (point(vcurrent));
        else if (allow_border_propagation && has_mesh_vertex(vcurrent))
        {
          Point_2 ref = midpoint (current);
          std::size_t this_section = find_section_of_point_from_vertex_view (vcurrent, ref);
          if(has_defined_height(mesh_vertex(vcurrent, this_section))
             && vertex_done.insert (point(vcurrent, this_section)).second)
            neighborhood[section].push_back (point(vcurrent, this_section));
        }

      }

      if (!still_in)
        continue;

      for (std::size_t i = 0; i < 3; ++ i)
      {
        Face_handle neighbor = current->neighbor(i);
        if (!(m_cdt.is_constrained (std::make_pair (current, i))))
          todo.push (std::make_pair(neighbor, section));
      }
      
    }
  }

  double estimate_height_with_pca (Vertex_handle vh, std::vector<Point_3>& neighborhood)
  {
    Plane_3 plane;
    Point_3 centroid;

    CGAL::linear_least_squares_fitting_3 (neighborhood.begin(), neighborhood.end(),
                                          plane, centroid, CGAL::Dimension_tag<0>());

    Line_3 line (Point_3 (vh->point().x(), vh->point().y(), 0.),
                 Vector_3 (0., 0., 1.));

    typename CGAL::cpp11::result_of<typename Kernel::Intersect_3(Line_3, Plane_3)>::type
      result = CGAL::intersection(line, plane);
    Point_3* inter;
    if (result && (inter = boost::get<Point_3>(&*result)))
      return inter->z();
    
    return std::numeric_limits<double>::max();
  }

  void estimate_missing_heights (Vertex_handle vh, double epsilon, bool allow_border_propagation = false)
  {
    std::vector<std::vector<Point_3> > small_neighborhood;
    get_neighborhood (vh, 3. * epsilon, small_neighborhood, allow_border_propagation);
    
    for (std::size_t i = 0; i < vh->info().size(); ++ i)
    {
      if (has_defined_height(vh->info()[i].second))
        continue;

      Point_3 new_point = m_mesh.point(vh->info()[i].second);

      std::vector<Point_3>& neighborhood = small_neighborhood[i];


      // If no neighbors, give up for now
      if (neighborhood.empty())
        continue;

      // Reference height = closest neighbor's height
      Point_3 closest(std::numeric_limits<double>::quiet_NaN(),
                      std::numeric_limits<double>::quiet_NaN(),
                      std::numeric_limits<double>::quiet_NaN());
      
      double dist_min = std::numeric_limits<double>::max();

      for (std::size_t j = 0; j < neighborhood.size(); ++ j)
      {
        double dist = CGAL::squared_distance (vh->point(), Point_2(neighborhood[j].x(),
                                                                   neighborhood[j].y()));
        if (dist < dist_min)
        {
          dist_min = dist;
          closest = neighborhood[j];
        }
      }
      dist_min = std::sqrt (dist_min);

      double ref_z = closest.z();

      // Try PCA if not too far from reference
      double pca_z = estimate_height_with_pca (vh, neighborhood);

      if (CGAL::abs(ref_z - pca_z) < (epsilon + dist_min))
        new_point = Point_3 (new_point.x(), new_point.y(), pca_z);
      else
        new_point = Point_3 (new_point.x(), new_point.y(), ref_z);

      m_mesh.point(vh->info()[i].second) = new_point;
    }

  }

  void estimate_missing_heights_from_planes (Vertex_handle vh)
  {
    std::vector<std::vector<std::size_t> > incident_planes (vh->info().size());
    
    Face_circulator circ = m_cdt.incident_faces (vh), start = circ;
    do
    {
      if (!is_ignored(circ))
      {
        Point_2 p = midpoint(circ);
        std::size_t i = find_section_of_point_from_vertex_view (vh, p);
        incident_planes[i].push_back (circ->info().plane_index);
      }
      ++ circ;
    }
    while (circ != start);

    for (std::size_t i = 0; i < incident_planes.size(); ++ i)
    {
//      CGAL_assertion (!incident_planes[i].empty());
      
      double z = 0.;
      for (std::size_t j = 0; j < incident_planes[i].size(); ++ j)
      {
        double d = estimate_height_with_plane(vh, incident_planes[i][j]);
        CGAL_assertion (d != std::numeric_limits<double>::max());
        z += d;
      }
      z /= incident_planes[i].size();

      Point_3 new_point = m_mesh.point(vh->info()[i].second);
      m_mesh.point(vh->info()[i].second) = Point_3 (new_point.x(), new_point.y(), z);
    }
  }

  template <typename PlaneIndexRange>
  double estimate_height_with_planes (const Point_2& pt, const PlaneIndexRange& planes,
                                      double epsilon)
  {
    return weighted_barycenter_of_projections_on_planes (pt, planes).z();
    
    Point_3 ref = barycenter_of_projections_on_planes (pt, planes);
    Point_3 candidate = least_squares_plane_intersection(planes);

    if (CGAL::abs(ref.z() - candidate.z()) < epsilon * epsilon)
      ref = candidate;
    
    return ref.z();
  }
  
  double estimate_height_with_plane (const Point_2& pt, std::size_t plane_index)
  {
    const Plane_3& plane = m_planes[plane_index];
    Line_3 line (Point_3 (pt.x(), pt.y(), 0.),
                 Vector_3 (0., 0., 1.));

    typename CGAL::cpp11::result_of<typename Kernel::Intersect_3(Line_3, Plane_3)>::type
      result = CGAL::intersection(line, plane);
    Point_3* inter;
    if (result && (inter = boost::get<Point_3>(&*result)))
      return inter->z();
    
    return std::numeric_limits<double>::max();
  }
  
  double estimate_height_with_plane (Vertex_handle vh, std::size_t plane_index)
  {
    return estimate_height_with_plane (vh->point(), plane_index);
  }

  void estimate_missing_heights_no_limit (Vertex_handle vh)
  {
    std::vector<bool> undefined (vh->info().size(), true);
    for (std::size_t i = 0; i < vh->info().size(); ++ i)
      if (has_defined_height(vh->info()[i].second))
        undefined[i] = false;

    std::queue<std::pair<Face_handle, std::size_t> > todo;
    Face_circulator circ = m_cdt.incident_faces (vh), start = circ;
    do
    {
      if (!is_ignored(circ))
      {
        Point_2 p = midpoint(circ);
        std::size_t i = find_section_of_point_from_vertex_view (vh, p);
        if (undefined[i])
          todo.push (std::make_pair (circ, i));
      }
      ++ circ;
    }
    while (circ != start);
    
    cpp11::unordered_set<Face_handle> done;
    
    while (!todo.empty())
    {
      Face_handle current = todo.front().first;
      std::size_t section = todo.front().second;
      todo.pop();
      if (!undefined[section])
        continue;

      if (!(done.insert(current).second) ||
          m_cdt.is_infinite(current) ||
          is_ignored(current))
        continue;

      for (std::size_t i = 0; i < 3; ++ i)
      {
        Vertex_handle vcurrent = current->vertex(i);

        if (has_unique_mesh_vertex(vcurrent))
        {
          Point_3 new_point = m_mesh.point(vh->info()[section].second);
          new_point = Point_3 (new_point.x(), new_point.y(), point(vcurrent).z());
          m_mesh.point (vh->info()[section].second) = new_point;
          undefined[section] = false;
          break;
        }
        else if (has_mesh_vertex(vcurrent))
        {
          Point_2 ref = midpoint (current);
          std::size_t this_section = find_section_of_point_from_vertex_view (vcurrent, ref);
          if(has_defined_height(mesh_vertex(vcurrent, this_section)))
          {
            // FOUND
            Point_3 new_point = m_mesh.point(vh->info()[section].second);
            new_point = Point_3 (new_point.x(), new_point.y(), point(vcurrent, this_section).z());
            m_mesh.point (vh->info()[section].second) = new_point;
            undefined[section] = false;
            break;
          }
        }
      }

      for (std::size_t i = 0; i < 3; ++ i)
      {
        Face_handle neighbor = current->neighbor(i);
        if (!(m_cdt.is_constrained (std::make_pair (current, i))))
          todo.push (std::make_pair(neighbor, section));
      }
      
    }
  }

  bool intersection_line_of_2_planes (std::size_t a, std::size_t b, Line_3& line)
  {
    const Plane_3& p0 = m_planes[a];
    const Plane_3& p1 = m_planes[b];
      
    typename cpp11::result_of<typename Kernel::Intersect_3(Plane_3, Plane_3)>::type
      result = CGAL::intersection(p0, p1);
      
    Line_3* inter;
    if (result && (inter = boost::get<Line_3>(&*result)))
    {
      line = *inter;
      return true;
    }

    return false;
  }

  template <typename PlaneIndexRange>
  Point_3 barycenter_of_projections_on_planes (Vertex_handle vh, const PlaneIndexRange& planes)
  {
    return barycenter_of_projections_on_planes (vh->point(), planes);
  }

  template <typename PlaneIndexRange>
  Point_3 barycenter_of_projections_on_planes (const Point_2& point, const PlaneIndexRange& planes)
  {
    Point_3 out (point.x(), point.y(), 0.);
    std::size_t nb = 0;
    
    Line_3 line (out, typename GeomTraits::Vector_3 (0., 0., 1.));
    
    BOOST_FOREACH (std::size_t idx, planes)
    {
      Plane_3& plane = m_planes[idx];
      typename CGAL::cpp11::result_of<typename GeomTraits::Intersect_3
                                      (Line_3,
                                       Plane_3)>::type
        result = CGAL::intersection(line, plane);
      typename GeomTraits::Point_3* inter;
      if (result && (inter = boost::get<typename GeomTraits::Point_3>(&*result)))
      {
        out = CGAL::barycenter (*inter, 1, out, nb);
        ++ nb;
      }
      else
        CGAL_assertion(false);
    }

    return out;
  }

  template <typename PlaneIndexRange>
  Point_3 weighted_barycenter_of_projections_on_planes (const Point_2& point, const PlaneIndexRange& planes)
  {
    Point_3 out (point.x(), point.y(), 0.);
    
    Line_3 line (out, typename GeomTraits::Vector_3 (0., 0., 1.));

    double accu = 0.;
    
    BOOST_FOREACH (std::size_t idx, planes)
    {
      Plane_3& plane = m_planes[idx];
      typename CGAL::cpp11::result_of<typename GeomTraits::Intersect_3
                                      (Line_3,
                                       Plane_3)>::type
        result = CGAL::intersection(line, plane);
      typename GeomTraits::Point_3* inter;
      if (result && (inter = boost::get<typename GeomTraits::Point_3>(&*result)))
      {
        out = CGAL::barycenter (*inter, m_plane_weights[idx], out, accu);
        accu += m_plane_weights[idx];
      }
      else
        CGAL_assertion(false);
    }

    return out;
  }

  template <typename PlaneIndexRange>
  Point_3 least_squares_plane_intersection (const PlaneIndexRange& incident_planes)
  {
    CGAL::Eigen_matrix<double> mat(incident_planes.size(), 3);
    CGAL::Eigen_vector<double> vec(incident_planes.size());

    std::size_t idx = 0;
    BOOST_FOREACH (std::size_t plane_index, incident_planes)
    {
      const Plane_3& plane = m_planes[plane_index];
      mat.set(idx,0, plane.a()); mat.set(idx,1, plane.b()); mat.set(idx,2, plane.c());
      vec.set(idx, -plane.d());
      ++ idx;
    }

    CGAL::Eigen_svd::solve(mat, vec);

    return Point_3 (vec(0), vec(1), vec(2));
  }

  Point_3 barycenter_plane_intersection (Vertex_handle vh,
                                         const std::vector<std::size_t>& incident_planes)
  {
    Point_3 point (vh->point().x(), vh->point().y(), 0.);
    std::size_t idx = 0;
    Line_3 line (point, Vector_3(0., 0., 1.));
    
    BOOST_FOREACH (std::size_t plane_index, incident_planes)
    {
      const Plane_3& plane = m_mesh.planes()[plane_index];
      
      typename cpp11::result_of<typename Kernel::Intersect_3(Line_3, Plane_3)>::type
        result = CGAL::intersection(line, plane);
      
      Point_3* inter;
      if (result && (inter = boost::get<Point_3>(&*result)))
      {
        point = CGAL::barycenter (*inter, 1, point, idx);
        std::cerr << point << std::endl;
      }
      else
        CGAL_assertion (false);

      ++ idx;
    }

    CGAL_assertion (point.x() == point.x());

    return point;
  }


  std::size_t segment_into_superfacets()
  {
    std::size_t current_index = 0;
    
    for (Finite_faces_iterator it = m_cdt.finite_faces_begin();
         it != m_cdt.finite_faces_end(); ++ it)
    {
      if (!is_default(it))
        continue;

      std::queue<Face_handle> todo;
      todo.push (it);

      while (!todo.empty())
      {
        Face_handle current = todo.front();
        todo.pop();

        if (!is_default(current))
          continue;

        current->info().index = Face_index(current_index);

        for (std::size_t i = 0; i < 3; ++ i)
        {
          if (!is_constrained (std::make_pair (current, i))
              && is_default(current->neighbor(i))
              && !m_cdt.is_infinite(current->neighbor(i)))
            todo.push(current->neighbor(i));
        }
      }
      
      ++ current_index;
    }

    return current_index;
  }

  void apply_planes_and_reset_superfacets (const std::vector<std::size_t>& indices)
  {
    m_plane_weights.resize (m_planes.size(), 0.);
    
    for (Finite_faces_iterator it = m_cdt.finite_faces_begin();
         it != m_cdt.finite_faces_end(); ++ it)
      if (!is_default(it))
      {
        it->info().plane_index = indices[std::size_t(it->info().index)];
        it->info().index = Face_index();
        m_plane_weights[it->info().plane_index]
          += CGAL::abs (CGAL::area (it->vertex(0)->point(),
                                    it->vertex(1)->point(),
                                    it->vertex(2)->point()));
      }
  }
  

  void check_structure_integrity()
  {
#ifdef TOP_VIEW_CHECK_STRUCTURE
    static int nb = 0;
    ++ nb;
    std::cerr << "INTEGRITY CHECK " << nb << " BEGIN" << std::endl;
    std::cerr.precision(18);

    m_mesh.is_valid (true);
    
    for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin();
         it != m_cdt.finite_vertices_end(); ++ it)
    {
      Vertex_handle vh = it;
      for (std::size_t i = 0; i < vh->info().size(); ++ i)
        if (m_v2v_map[vh->info()[i].second] != vh)
          std::cerr << "  [Bad structure] Mesh vertex not connected to correct CDT vertex" << std::endl;
    }

    BOOST_FOREACH (Vertex_index vi, m_mesh.vertices())
    {
      Vertex_handle vh = m_v2v_map[vi];
      if (vh == Vertex_handle())
      {
        std::cerr << "  [Bad structure] Mesh vertex not connected to any CDT vertex" << std::endl;
        continue;
      }
      
      if (m_nvv_map[vi] == Vertex_index())
        continue;
      
      cpp11::unordered_set<Vertex_index> seen;
      Vertex_index circ = vi;
      do
      {
        if (!seen.insert (circ).second)
        {
          std::cerr << "  [Bad structure] Internal loop on vertical vertex circulator" << std::endl;
          std::ofstream file ("struct.xyz");
          file.precision(18);
          BOOST_FOREACH (Vertex_index vi, seen)
            file << point(vi) << std::endl;
          abort();
          break;
        }
        if (circ == Vertex_index())
        {
          std::cerr << "  [Bad structure] Uninitialized index on vertical vertex circulator" << std::endl;
          break;
        }
        circ = m_nvv_map[circ];
      }
      while (circ != vi);
    }

    std::ofstream uf ("ufaces.xyz");
    uf.precision(18);
    
    for (Finite_faces_iterator it = m_cdt.finite_faces_begin();
         it != m_cdt.finite_faces_end(); ++ it)
    {
      Face_handle fh = it;
      if (is_ignored (fh))
        continue;
      if (is_default (fh))
      {
//        std::cerr << "  [Bad structure] Uninitialized face" << std::endl;
        for (std::size_t i = 0; i < 3; ++ i)
          uf << it->vertex(i)->point() << " 0" << std::endl;
        continue;
      }

      if (m_f2f_map[fh->info().index] != fh)
        std::cerr << "  [Bad structure] Mesh face not connected to correct CDT face" << std::endl;
    }
           
    BOOST_FOREACH (Face_index fi, m_mesh.faces())
    {
      Face_handle fh = m_f2f_map[fi];
      if (fh == Face_handle())
        continue;

      if (fh->info().index != fi)
        std::cerr << "  [Bad structure] CDT face not connected to correct mesh face" << std::endl;

      BOOST_FOREACH (Vertex_index vi, vertices_around_face (halfedge (fi, m_mesh), m_mesh))
      {
        Vertex_handle vh = m_v2v_map[vi];

        if (vh == Vertex_handle())
        {
          std::cerr << "  [Bad structure] Mesh face connected to mesh vertex without CDT vertex" << std::endl;
          continue;
        }

        bool found = false;
        for (std::size_t i = 0; i < 3; ++ i)
          if (fh->vertex(i) == vh)
          {
            found = true;
            break;
          }

        if (!found)
          std::cerr << "  [Bad structure] Vertices of mesh face do not correpond to vertices of CDT face" << std::endl;
      }
    }
    std::cerr << "INTEGRITY CHECK " << nb << " END" << std::endl;
#endif
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

  void DEBUG_dump_off_5() 
  {
    for (typename Mesh::Vertex_range::iterator it = m_mesh.vertices().begin();
         it != m_mesh.vertices().end(); ++ it)
    {
      const Point_3& p = m_mesh.point(*it);
      if (p.z() != p.z())
        m_mesh.point(*it) = Point_3 (p.x(), p.y(), 0);
    }
      
    
    std::ofstream f("debug5.off");
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
      if (is_wall_buffer(it))
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
      if (is_wall_buffer(it))
      {
        f << "3";
        for (std::size_t i = 0; i < 3; ++ i)
          f << " " << map[it->vertex(i)] << std::endl;
        f << std::endl;
      }
  }


  void DEBUG_dump_ply_3() const
  {
    std::ofstream f("debug3.ply");
    f.precision(18);
    
    f << "ply" << std::endl
      << "format ascii 1.0" << std::endl
      << "element vertex " << m_cdt.number_of_vertices() << std::endl
      << "property double x" << std::endl
      << "property double y" << std::endl
      << "property double z" << std::endl
      << "element face " << m_cdt.number_of_faces() << std::endl
      << "property list uchar int vertex_indices" << std::endl
      << "property uchar red" << std::endl
      << "property uchar green" << std::endl
      << "property uchar blue" << std::endl
      << "end_header" << std::endl;

    std::map<Vertex_handle, std::size_t> map;
    std::size_t idx = 0;
    for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin(); it != m_cdt.finite_vertices_end(); ++ it)
    {
      f << it->point() << " 0" << std::endl;
      map[it] = idx ++;
    }

    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
    {
      int red = 128, green = 128, blue = 128;
      if (is_ignored(it)) { red = 0; green = 0; blue = 0; }
      else if (is_wall_buffer(it)) { red = 230; green = 0; blue = 0; }
      else if (is_buffer(it)) { red = 0; green = 0; blue = 230; }
        
      f << "3";
      for (std::size_t i = 0; i < 3; ++ i)
        f << " " << map[it->vertex(i)];
      f << " " << red << " " << green << " " << blue << std::endl;
    }
  }

  void DEBUG_dump_off_4() const
  {
    std::ofstream f("debug4.off");
    f.precision(18);
    
    std::size_t nb_faces = 0;
    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
      if(has_mesh_face(it))
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
      if(has_mesh_face(it))
      {
        f << "3";
        for (std::size_t i = 0; i < 3; ++ i)
          f << " " << map[it->vertex(i)] << std::endl;
        f << std::endl;
      }
    }
  }

  void DEBUG_dump_ply(const char* filename) const
  {
    std::ofstream f(filename);
    f.precision(18);
    
    f << "ply" << std::endl
      << "format ascii 1.0" << std::endl
      << "element vertex " << m_cdt.number_of_vertices() << std::endl
      << "property double x" << std::endl
      << "property double y" << std::endl
      << "property double z" << std::endl
      << "element face " << m_cdt.number_of_faces() << std::endl
      << "property list uchar int vertex_indices" << std::endl
      << "property uchar red" << std::endl
      << "property uchar green" << std::endl
      << "property uchar blue" << std::endl
      << "end_header" << std::endl;

    std::map<Vertex_handle, std::size_t> map;
    std::size_t idx = 0;
    for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin(); it != m_cdt.finite_vertices_end(); ++ it)
    {
      f << it->point() << " 0" << std::endl;
      map[it] = idx ++;
    }
    
    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
    {
      int red = 0, green = 0, blue = 0;
      
      if (is_default_buffer(it))
      {
        red = 0;
        green = 128;
        blue = 0;
      }
      else if (is_wall_buffer(it))
      {
        red = 0;
        green = 0;
        blue = 128;
      }
      else if (is_ridge_buffer(it))
      {
        red = 128;
        green = 0;
        blue = 0;
      }
      else if (is_ignored_buffer(it))
      {
        red = 128;
        green = 128;
        blue = 128;
      }
      else if (it->info().has_plane())
      {
        srand(it->info().plane_index);
        red = 64 + rand() % 128;
        green = 64 + rand() % 128;
        blue = 64 + rand() % 128;
      }
               
      f << "3";
      for (std::size_t i = 0; i < 3; ++ i)
        f << " " << map[it->vertex(i)];
      f << " " << red << " " << green << " " << blue << std::endl;
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

  void DEBUG_dump_edges()
  {
    Vector_3 vertical (0., 0., 1.);
    
    std::ofstream f("edges.polylines.txt");
    f.precision(18);

    BOOST_FOREACH (Edge_index ei, m_mesh.edges())
    {
      Vertex_index v0 = m_mesh.vertex(ei, 0);
      Vertex_index v1 = m_mesh.vertex(ei, 1);
      const Point_3& p0 = point(v0);
      const Point_3& p1 = point(v1);

      Vertex_handle vh0 = cdt_vertex (v0);
      Vertex_handle vh1 = cdt_vertex (v1);

      CGAL_assertion (vh0 != Vertex_handle() && vh1 != Vertex_handle());

      bool edge = false;
      
      Face_handle fh;
      int idx;
      if (vh0 == vh1)
        edge = true;
      else if (is_edge(vh0, vh1, fh, idx)
               && is_constrained (std::make_pair(fh, idx)))
      {
        Face_index f0 = m_mesh.face (m_mesh.halfedge(ei, 0));
        Face_index f1 = m_mesh.face (m_mesh.halfedge(ei, 1));

        if (f0 != Face_index() && f1 != Face_index() &&
            (has_cdt_face(f0) || has_cdt_face(f1)))
          edge = true;
      }
      
      if (edge)
        f << "2 " << p0 << " " << p1 << std::endl;
    }
  }

  void DEBUG_dump_planes_inter()
  {
    std::ofstream f("edges_inter.polylines.txt");
    f.precision(18);

    typedef std::map<std::pair<std::size_t, std::size_t>, std::vector<Edge> > Intermap;

    Intermap intermap;

    for (Finite_edges_iterator it = m_cdt.finite_edges_begin ();
         it != m_cdt.finite_edges_end(); ++ it)
    {
      Face_handle f0 = it->first;
      Face_handle f1 = it->first->neighbor (it->second);

      if (!f0->info().has_plane() || !f1->info().has_plane() ||
          is_ignored (f0) || is_ignored (f1))
        continue;

      std::size_t p0 = f0->info().plane_index;
      std::size_t p1 = f1->info().plane_index;
      if (p0 == p1)
        continue;
      if (p1 > p0)
        std::swap (p0, p1);

      typename Intermap::iterator found = intermap.insert(std::make_pair (std::make_pair(p0,p1),
                                                                          std::vector<Edge>())).first;
      found->second.push_back (*it);
    }

    std::cerr << intermap.size() << " intersection(s) found" << std::endl;

    for (typename Intermap::iterator it = intermap.begin(); it != intermap.end(); ++ it)
    {
      const Plane_3& p0 = m_planes[it->first.first];
      const Plane_3& p1 = m_planes[it->first.second];

      typename CGAL::cpp11::result_of<typename Kernel::Intersect_3(Plane_3, Plane_3)>::type
        result = CGAL::intersection(p0, p1);
      Line_3* inter;
      if (result && (inter = boost::get<Line_3>(&*result)))
      {
        Point_3 orig = inter->point();
        Vector_3 vec = inter->to_vector();
        vec /= std::sqrt (vec*vec);

        double min = std::numeric_limits<double>::max();
        double max = -std::numeric_limits<double>::max();
        
        for (std::size_t i = 0; i < it->second.size(); ++ i)
        {
          for (std::size_t j = 1; j <= 2; ++ j)
          {
            Vertex_handle v = it->second[i].first->vertex ((it->second[i].second + j)%3);
            if (has_unique_mesh_vertex(v))
            {
              const Point_3& pt = point(v);
              double val = (inter->projection(pt) - orig) * vec;
              min = (std::min)(min, val);
              max = (std::max)(max, val);
            }
          }
        }

        if (min != std::numeric_limits<double>::max() &&
            max != -std::numeric_limits<double>::max())
          f << "2 " << orig + min * vec << " " << orig + max * vec << std::endl;          
      }
    }
    
  }
};


}
  
#endif // CGAL_TVSR_SURFACE_MESH_ON_CDT
