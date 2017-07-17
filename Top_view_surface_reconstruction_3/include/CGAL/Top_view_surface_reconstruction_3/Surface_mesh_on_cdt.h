#ifndef CGAL_TVSR_SURFACE_MESH_ON_CDT
#define CGAL_TVSR_SURFACE_MESH_ON_CDT

#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <CGAL/Surface_mesh.h>

#include <fstream>

namespace CGAL
{

template <typename GeomTraits>
class Surface_mesh_on_cdt
{
public:
  
  typedef GeomTraits Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Point_2 Point_2;
  typedef Surface_mesh<Point_3> Mesh;
  typedef typename Mesh::Vertex_index Vertex_index;
  typedef typename Mesh::Face_index Face_index;

  typedef Triangulation_vertex_base_with_info_2<std::vector<Vertex_index>, Kernel>  Vbwi;
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
  Vertex_index mesh_vertex (Vertex_handle vh) const { return vh->info()[0]; }

  Face_handle cdt_face (Face_index fi) const { return m_f2f_map[fi]; }
  Face_index mesh_face (Face_handle fh) const { return fh->info(); }

  Face_handle locate (const Point_2& point, Face_handle hint = Face_handle()) const
  { return m_cdt.locate (point, hint); }
  typename GeomTraits::Triangle_2 triangle (Face_handle fh) const { return m_cdt.triangle(fh); }

  Finite_faces_iterator finite_faces_begin() { return m_cdt.finite_faces_begin(); }
  Finite_faces_iterator finite_faces_end() { return m_cdt.finite_faces_end(); }
  Finite_vertices_iterator finite_vertices_begin() { return m_cdt.finite_vertices_begin(); }
  Finite_vertices_iterator finite_vertices_end() { return m_cdt.finite_vertices_end(); }
  Finite_edges_iterator finite_edges_begin() { return m_cdt.finite_edges_begin(); }
  Finite_edges_iterator finite_edges_end() { return m_cdt.finite_edges_end(); }
  Face_circulator incident_faces (Vertex_handle vh) { return m_cdt.incident_faces (vh); }
  void insert_constraint (Vertex_handle v0, Vertex_handle v1) { m_cdt.insert_constraint (v0, v1); }
  bool are_there_incident_constraints (Vertex_handle vh) { return m_cdt.are_there_incident_constraints (vh); }
                                  
  bool is_infinite (Face_handle fh) const { return m_cdt.is_infinite(fh); }

  bool has_mesh_face (Face_handle fh) const { return (int(fh->info()) >= 0); }
  bool is_skiped (Face_handle fh) const { return (int(fh->info()) == -2); }
  void skip (Face_handle fh) { fh->info() = Face_index(-2); }

  const Point_3& point (Vertex_index vi) const { return m_mesh.point(vi); }
  
  void insert (const Point_3& point)
  {
    Vertex_handle vh = m_cdt.insert (Point_2 (point.x(), point.y()));
    Vertex_index vi = m_mesh.add_vertex (point);
    vh->info().push_back(vi);
    m_v2v_map[vi] = vh;
  }

  void remove (Vertex_handle vh)
  {
    m_cdt.remove (vh);
  }

  void add_face (Face_handle fh)
  {
    Face_index fi = m_mesh.add_face (fh->vertex(0)->info()[0],
                                     fh->vertex(1)->info()[0],
                                     fh->vertex(2)->info()[0]);
    fh->info() = fi;
    m_f2f_map[fi] = fh;
  }


  void DEBUG_dump_off_1() const
  {
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

    for (Finite_vertices_iterator it = m_cdt.finite_vertices_begin(); it != m_cdt.finite_vertices_end(); ++ it)
      f << it->point() << " 0" << std::endl;

    for (Finite_faces_iterator it = m_cdt.finite_faces_begin(); it != m_cdt.finite_faces_end(); ++ it)
      if (has_mesh_face(it))
      {
        f << "3";
        for (std::size_t i = 0; i < 3; ++ i)
          f << " " << std::size_t(it->vertex(i)->info()[0]);
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


};


}
  
#endif // CGAL_TVSR_SURFACE_MESH_ON_CDT
