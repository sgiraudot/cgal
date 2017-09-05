#ifndef CGAL_TVSR_SURFACE_MESH_ON_CDT
#define CGAL_TVSR_SURFACE_MESH_ON_CDT

#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Random.h>

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

  struct Face_info
  {
    Face_index index;
    Point_2 endpoint;
    std::vector<std::size_t> incident_lines;
    Face_info()
      : index()
      , endpoint (std::numeric_limits<double>::quiet_NaN(),
                  std::numeric_limits<double>::quiet_NaN())
    { }

    bool has_endpoint() const { return (endpoint.x() == endpoint.x()); }
    void erase_endpoint() { endpoint = Point_2(std::numeric_limits<double>::quiet_NaN(),
                                               std::numeric_limits<double>::quiet_NaN()); }
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
  typename GeomTraits::Triangle_2 triangle (Face_handle fh) const { return m_cdt.triangle(fh); }

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
  bool are_there_incident_constraints (Vertex_handle vh) { return m_cdt.are_there_incident_constraints (vh); }

  template <typename OutputItEdges>
  OutputItEdges incident_constraints (Vertex_handle vh, OutputItEdges output) const
  {
    return m_cdt.incident_constraints (vh, output);
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
  
  bool has_mesh_face (Face_handle fh) const { return (int(fh->info().index) >= 0); }
  bool is_default (Face_handle fh) const { return (int(fh->info().index) == -1); }
  void make_default (Face_handle fh) { fh->info().index = Face_index(-1); }
  bool is_buffer (Face_handle fh) const { return (int(fh->info().index) == -2) || int(fh->info().index) == -4; }
  void make_buffer (Face_handle fh) { fh->info().index = Face_index(-2); }
  bool is_ignored (Face_handle fh) const { return (int(fh->info().index) == -3); }
  void make_ignored (Face_handle fh) { fh->info().index = Face_index(-3); }
  bool is_unhandled_buffer (Face_handle fh) const { return (int(fh->info().index) == -2); }
  bool is_handled_buffer (Face_handle fh) const { return (int(fh->info().index) == -4); }
  void make_handled_buffer (Face_handle fh) { fh->info().index = Face_index(-4); }

  int cw (int i) const { return m_cdt.cw(i); }
  int ccw (int i) const { return m_cdt.ccw(i); }

  const Point_3& point (Vertex_index vi) const { return m_mesh.point(vi); }
  Point_3& point (Vertex_index vi) { return m_mesh.point(vi); }
  const Point_3& point (Vertex_handle vh, std::size_t idx = 0) const { return m_mesh.point(mesh_vertex(vh, idx)); }
  Point_3& point (Vertex_handle vh, std::size_t idx = 0) { return m_mesh.point(mesh_vertex(vh, idx)); }

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
    Vertex_index vi = m_mesh.add_vertex (point);
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
    else
    {
      std::ofstream file ("bad.off");
      file.precision(18);
      file << "OFF\n3 1 0\n"
           << point(a) << std::endl << point(b) << std::endl << point(c) << std::endl
           << "3 0 1 2" << std::endl;

      Vertex_index v[3]; v[0] = a; v[1] = b; v[2] = c;

      for (std::size_t vn = 0; vn < 3; ++ vn)
      {
        std::set<Face_index> faces;
        BOOST_FOREACH (Face_index idx, faces_around_target (halfedge (v[vn], m_mesh), m_mesh))
          if (idx != Face_index())
            faces.insert (idx);

        std::ostringstream oss;
        oss << "bad_" << vn << ".off";
        std::ofstream file2 (oss.str().c_str());
        
        file2.precision(18);

        file2 << "OFF" << std::endl
              << faces.size() * 3 << " " << faces.size() << " 0" << std::endl;
        BOOST_FOREACH (Face_index fi, faces)
        {
          std::cerr << fi;
          BOOST_FOREACH (Vertex_index vi, vertices_around_face (halfedge(fi, m_mesh), m_mesh))
          {
            std::cerr << " " << vi;
            file2 << point(vi) << std::endl;
          }
          std::cerr << std::endl;
        }

        for (std::size_t i = 0; i < faces.size(); ++ i)
          file2 << "3 " << 3*i << " " << 3*i + 1 << " " << 3*i + 2 << std::endl;
      }
      
      exit(0);
    }
  }

  bool add_face (Vertex_index a, Vertex_index b, Vertex_index c)
  {
    Face_index fi =
      m_mesh.add_face (a, b, c);

    return (fi != Face_index());
//    std::cerr << fi << " ";
    // if (fi >= m_mesh.number_of_faces())
    //   std::cerr << "WHAT?!" << std::endl;
    // m_f2f_map[fi] = Face_handle();
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
    if (vh->info().size() == 1)
      return 0;

    Direction_2 direction (Vector_2 (vh->point(), point));
    for (std::size_t i = 0; i < vh->info().size(); ++ i)
      if (direction.counterclockwise_in_between (vh->info()[i].first,
                                                 vh->info()[(i+1)%(vh->info().size())].first))
        return i;

    std::cerr << "Warning section" << std::endl;
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
    
    std::set<Face_handle> done;
    std::set<Point_3> vertex_done;
    
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

      Point_3 new_point = m_mesh.point(vh->info()[i].second);

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
    
    std::set<Face_handle> done;
    
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
      
      std::set<Vertex_index> seen;
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
      if (is_handled_buffer(it))
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
      if (is_handled_buffer(it))
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
