#ifndef CGAL_CLASSIFICATION_DEM_H
#define CGAL_CLASSIFICATION_DEM_H

#include <fstream>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/remove_outliers.h>

#if defined(CGAL_DEM_VERBOSE)
#define CGAL_DEM_SILENT false
#else
#define CGAL_DEM_SILENT true
#endif

#define CGAL_DEM_CERR \
  if(CGAL_DEM_SILENT) {} else std::cerr

#ifdef CGAL_LINKED_WITH_TBB
#  include <tbb/parallel_for.h>
#  include <tbb/blocked_range.h>
#  include <tbb/scalable_allocator.h>
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL
{

namespace Classification
{

template <typename Kernel, typename PointRange, typename PointMap>
class Digital_elevation_model
{
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Triangle_3 Triangle_3;
  typedef typename boost::property_traits<PointMap>::key_type Point_item;
  typedef typename boost::property_traits<PointMap>::reference reference;

  struct Face_info
  {
    std::vector<size_t> inliers;
    size_t idx;
    double cos_angle;
    double sq_distance;
  };

  typedef boost::int32_t size_t;
  typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, Kernel> Vbi;
  typedef CGAL::Triangulation_face_base_with_info_2<Face_info, Kernel> Fbi;
  typedef CGAL::Triangulation_data_structure_2<Vbi, Fbi> Tds;
  typedef CGAL::Delaunay_triangulation_2<Kernel,Tds> DT;

  typedef typename DT::Vertex_handle Vertex_handle;
  typedef typename DT::Face_handle Face_handle;
  typedef typename DT::Vertex_circulator Vertex_circulator;
  typedef typename DT::Face_circulator Face_circulator;
  typedef typename DT::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename DT::Finite_faces_iterator Finite_faces_iterator;

  const PointRange& m_points;
  PointMap m_point_map;
  std::vector<bool> m_is_ground;
  std::vector<Point_3> m_bbox_pts;
  
  double m_max_building_size;
  double m_min_terrain_cos_angle;
  double m_min_cos_angle;
  double m_max_sq_distance;
  double m_min_edge_length;
  
  mutable std::vector<double> m_cos_angle_histogram;
  mutable std::vector<double> m_sq_dist_histogram;

  DT m_dt;
  mutable Face_handle m_hint;

  struct Index_to_point_map
  {
    typedef size_t key_type;
    typedef typename boost::property_traits<PointMap>::value_type value_type;
    typedef typename boost::property_traits<PointMap>::reference reference;
    typedef typename boost::property_traits<PointMap>::category category;

    const PointRange* points;
    const PointMap point_map;
    Index_to_point_map () : points (nullptr) { }
    Index_to_point_map (const PointRange& points, const PointMap point_map)
      : points (&points), point_map (point_map) {}

    friend reference get (const Index_to_point_map& map, size_t idx)
    {
      return get (map.point_map, *(map.points->begin() + idx));
    }
  };

  Index_to_point_map index_to_point_map() const { return Index_to_point_map(m_points, m_point_map); }

public:

  Digital_elevation_model (const PointRange& points,
                                 PointMap point_map,
                                 double max_building_size = 50.,
                                 double max_terrain_angle = 88,
                                 double max_angle = 10,
                                 double max_distance = 1.4,
                                 double min_edge_length = 2.0,
                                 unsigned int init_outlier_neighbors = 12,
                                 double init_outlier_percent = 25.,
                                 const std::function<void(const Digital_elevation_model&)>&
                                 callback = std::function<void(const Digital_elevation_model&)>())
    : m_points (points)
    , m_point_map (point_map)
    , m_max_building_size (max_building_size)
    , m_min_terrain_cos_angle (std::cos(CGAL_PI * (max_terrain_angle / 180.)))
    , m_min_cos_angle (std::cos(CGAL_PI * (max_angle / 180.)))
    , m_max_sq_distance (max_distance * max_distance)
    , m_min_edge_length (min_edge_length * min_edge_length)
  {
    init_triangulation(init_outlier_neighbors, init_outlier_percent);
//    return;

    std::vector<size_t> to_add;

    Face_handle fh;

    std::size_t iter = 0;
    std::size_t prev_nb_ground = 0;

    for (Finite_faces_iterator it = m_dt.finite_faces_begin();
         it != m_dt.finite_faces_end(); ++ it)
      for (size_t idx : it->info().inliers)
      {
        double cos_angle;
        double sq_distance;
        m_is_ground[idx] = is_ground(idx, it);
        ++ prev_nb_ground;
      }
    
    while (true)
    {
      ++ iter;
      CGAL_DEM_CERR << "Iteration " << iter << ":" << std::endl;
      
      std::size_t nb_ground = 0;
      for (const bool b : m_is_ground)
        if (b)
          ++ nb_ground;

      CGAL_DEM_CERR << " * " << nb_ground << " ground point(s), "
                << m_dt.number_of_faces() << " face(s) in DT" << std::endl;

      if (iter != 1 && nb_ground == prev_nb_ground) // Algorithm finished
      {
        CGAL_DEM_CERR << "No new ground point" << std::endl;
        break;
      }
      prev_nb_ground = nb_ground;

      for (Finite_faces_iterator it = m_dt.finite_faces_begin();
           it != m_dt.finite_faces_end(); ++ it)
        if (max_edge_length(it) > min_edge_length)
        {
#define STRATEGY_CENTROID // Best
//#define STRATEGY_RANDOM
          
#if defined(STRATEGY_RANDOM) // Close to best
          for (size_t idx : it->info().inliers)
            if (m_is_ground[idx])
            {
              to_add.push_back (idx);
              break;
            }
#elif defined(STRATEGY_CENTROID) // Best
          Point_2 centroid = CGAL::barycenter (it->vertex(0)->point(), 1,
                                               it->vertex(1)->point(), 1,
                                               it->vertex(2)->point(), 1);
          double dist_min = std::numeric_limits<double>::max();
          size_t chosen = no_element();
          for (size_t idx : it->info().inliers)
            if (m_is_ground[idx])
            {
              double dist = CGAL::squared_distance(centroid, point_2(idx));
              if (dist < dist_min)
              {
                dist_min = dist;
                chosen = idx;
              }
            }
          if (chosen != no_element())
            to_add.push_back (chosen);
#elif defined(STRATEGY_CIRCUMCENTER)
          Point_2 circumcenter = CGAL::circumcenter (it->vertex(0)->point(),
                                                     it->vertex(1)->point(),
                                                     it->vertex(2)->point());
          double dist_min = std::numeric_limits<double>::max();
          size_t chosen = no_element();
          for (size_t idx : it->info().inliers)
            if (m_is_ground[idx])
            {
              double dist = CGAL::squared_distance(circumcenter, point_2(idx));
              if (dist < dist_min)
              {
                dist_min = dist;
                chosen = idx;
              }
            }
          if (chosen != no_element())
            to_add.push_back (chosen);
#elif defined(STRATEGY_CLOSEST) // Close to best but really long
          Triangle_3 triangle = triangle_3(it);
          double dist_min = std::numeric_limits<double>::max();
          size_t chosen = no_element();
          for (size_t idx : it->info().inliers)
            if (m_is_ground[idx])
            {
              double dist = CGAL::squared_distance(triangle, point_3(idx));
              if (dist < dist_min)
              {
                dist_min = dist;
                chosen = idx;
              }
            }
          if (chosen != no_element())
            to_add.push_back (chosen);
#elif defined(STRATEGY_FARTHEST) // Really bad
          Triangle_3 triangle = triangle_3(it);
          double dist_max = 0;
          size_t chosen = no_element();
          for (size_t idx : it->info().inliers)
            if (m_is_ground[idx])
            {
              double dist = CGAL::squared_distance(triangle, point_3(idx));
              if (dist > dist_max)
              {
                dist_max = dist;
                chosen = idx;
              }
            }
          if (chosen != no_element())
            to_add.push_back (chosen);
#endif
        }

      // m_min_cos_angle = 1 - ((1 - m_min_cos_angle) * 0.95);
      m_max_sq_distance *= 1;
      
      CGAL_DEM_CERR << " * cos_angle_limit = " << m_min_cos_angle << std::endl;
      CGAL_DEM_CERR << " * sq_distance_limit = " << m_max_sq_distance << std::endl;
      CGAL_DEM_CERR << " * " << to_add.size() << " point(s) to be added" << std::endl;

      // std::nth_element (m_cos_angle_histogram.begin(),
      //                   m_cos_angle_histogram.begin() + (m_cos_angle_histogram.size() / 2),
      //                   m_cos_angle_histogram.end());
      // std::nth_element (m_sq_dist_histogram.begin(),
      //                   m_sq_dist_histogram.begin() + (m_sq_dist_histogram.size() / 2),
      //                   m_sq_dist_histogram.end());

      // CGAL_DEM_CERR << " * median cos_angle = " << m_cos_angle_histogram[m_cos_angle_histogram.size() / 2] << std::endl
      //               << " * median sd_distance = " << m_sq_dist_histogram[m_sq_dist_histogram.size() / 2] << std::endl;
      // m_min_cos_angle = m_cos_angle_histogram[m_cos_angle_histogram.size() / 2];
      // m_max_sq_distance = m_sq_dist_histogram[m_sq_dist_histogram.size() / 2];
      
      // m_cos_angle_histogram.clear();
      // m_sq_dist_histogram.clear();

      callback(*this);
      
      if (to_add.size() == 0)
        break;

      Face_handle fh;
      for (size_t idx : to_add)
        fh = insert (idx, fh);

      to_add.clear();
    }
    //    std::cerr << m_regular << " " << m_mirror << std::endl;
    CGAL_DEM_CERR << "All done" << std::endl;
  }

  friend std::ostream& operator<< (std::ostream& os, const Digital_elevation_model& ptd)
  {
    std::size_t nb_vertices = 0, nb_faces = 0;
    
    for (Finite_vertices_iterator it = ptd.m_dt.finite_vertices_begin();
         it != ptd.m_dt.finite_vertices_end(); ++ it)
      if (it->info() >= 0)
        ++ nb_vertices;
    for (Finite_faces_iterator it = ptd.m_dt.finite_faces_begin();
         it != ptd.m_dt.finite_faces_end(); ++ it)
    {
      bool ok = true;
      for (int i = 0; i < 3; ++ i)
        if (it->vertex(i)->info() < 0)
        {
          ok = false;
          break;
        }
      if (ok)
        ++ nb_faces;
    }
    
    os << "OFF" << std::endl << nb_vertices << " " << nb_faces << " 0" << std::endl;
    std::unordered_map<size_t, size_t> map_v2v;
    size_t idx = 0;
    for (Finite_vertices_iterator it = ptd.m_dt.finite_vertices_begin();
         it != ptd.m_dt.finite_vertices_end(); ++ it)
    {
      if (it->info() >= 0)
      {
        os << ptd.point_3(it->info()) << std::endl;
        map_v2v.insert (std::make_pair (it->info(), idx ++));
      }
    }
    for (Finite_faces_iterator it = ptd.m_dt.finite_faces_begin();
         it != ptd.m_dt.finite_faces_end(); ++ it)
    {
      bool ok = true;
      for (int i = 0; i < 3; ++ i)
        if (it->vertex(i)->info() < 0)
        {
          ok = false;
          break;
        }
      if (!ok)
        continue;
      
      os << "3";
      for (int i = 0; i < 3; ++ i)
        os << " " << map_v2v[it->vertex(i)->info()];
      os << std::endl;
    }

    return os;
  }

  bool is_ground (typename PointRange::const_iterator it) const
  {
    return m_is_ground[it - m_points.begin()];
  }
  
  double distance (typename PointRange::const_iterator it) const
  {
    CGAL_assertion (m_points.begin() <= it && it < m_points.end());
    if (is_ground(it))
      return 0.;

    m_hint = m_dt.locate (point_2 (it - m_points.begin()), m_hint);

    return CGAL::approximate_sqrt(cos_angle_and_sq_distance (point_3 (it - m_points.begin()), m_hint).second)
      - CGAL::approximate_sqrt (m_max_sq_distance);
  }


private:

  inline size_t no_element() const { return std::numeric_limits<size_t>::max(); }

  inline Point_3 point_3 (const size_t& idx) const
  {
    if (idx < 0)
      return m_bbox_pts[-1 - idx];
    return get (m_point_map, *(m_points.begin() + idx));
  }
  inline Point_2 point_2 (const size_t& idx) const
  {
    if (idx < 0)
      return Point_2 (m_bbox_pts[-1 - idx].x(), m_bbox_pts[-1 - idx].y());
    return Point_2 (point_3(idx).x(), point_3(idx).y());
  }

  inline Triangle_3 triangle_3 (Face_handle fh) const
  {
    return Triangle_3 (point_3(fh->vertex(0)->info()),
                       point_3(fh->vertex(1)->info()),
                       point_3(fh->vertex(2)->info()));
  }

  double max_edge_length (Face_handle fh) const
  {
    double out = 0.;
    for (int i = 0; i < 3; ++ i)
    {
      double l = CGAL::squared_distance(fh->vertex(i)->point(), fh->vertex((i+1)%3)->point());
      out = (std::max) (l, out);
    }
    return out;
  }

  Face_handle path_independant_locate (const Point_2& point, Face_handle fh) const
  {
    typename DT::Locate_type lt;
    int li;
    fh = m_dt.locate (point, lt, li,fh);

    if (lt == DT::EDGE)
    {
      Face_handle candidate = fh->neighbor(li);
      if (candidate < fh)
        fh = candidate;
    }
    else if (lt == DT::VERTEX)
    {
      Face_circulator circ = m_dt.incident_faces(fh->vertex(li)),
        start = circ;
      do
      {
        Face_handle candidate = circ;
        if (candidate < fh)
          fh = candidate;
      }
      while (++ circ != start);
    }
    
    return fh;
  }

  bool is_ground (const size_t& idx, Face_handle fh) const
  {
    Vector_3 ortho = triangle_3(fh).supporting_plane().orthogonal_vector();
    ortho = ortho / CGAL::approximate_sqrt(ortho * ortho);
    
    if (CGAL::abs(ortho * Vector_3(0,0,1)) > m_min_terrain_cos_angle)
    {
      double cos_angle;
      double sq_distance;
      std::tie (cos_angle, sq_distance) = cos_angle_and_sq_distance(point_3(idx), fh);

      if (cos_angle > m_min_cos_angle && sq_distance < m_max_sq_distance)
      {
        // m_cos_angle_histogram.push_back(cos_angle);
        // m_sq_dist_histogram.push_back(sq_distance);
        return true;
      }
      // else
      return false;
    }
    // else
    if (mirror_criteria_met (idx, fh))
      return true;
    // else
    return false;
  }

  std::pair<double, double> cos_angle_and_sq_distance (const Point_3& point, Face_handle fh) const
  {
    Triangle_3 triangle = triangle_3 (fh);
    
    double distance = CGAL::squared_distance (point, triangle);

    double cos_angle = 1.;

    Point_3 projection = triangle.supporting_plane().projection(point);

    double min_dist = std::numeric_limits<double>::max();
    
    for (int i = 0; i < 3; ++ i)
    {
      Point_3 pvertex = point_3 (fh->vertex(i)->info());
      double dist = CGAL::squared_distance (point, pvertex);
      if (dist < min_dist)
      {
        min_dist = dist;
        Vector_3 vpoint (pvertex, point);
        Vector_3 vproj (pvertex, projection);

        if (pvertex == point)
          cos_angle = 1.;
        else if (pvertex == projection)
          cos_angle = 0.;
        else
        {
          vpoint = vpoint / CGAL::approximate_sqrt (vpoint * vpoint);
          vproj = vproj / CGAL::approximate_sqrt (vproj * vproj);
          cos_angle = vpoint * vproj;
        }
      }
    }

    return std::make_pair(cos_angle, distance);
  }

  bool mirror_criteria_met (size_t idx, Face_handle fh) const
  {
    Point_3 point = point_3(idx);
    
    double min_dist = std::numeric_limits<double>::max();
    Vertex_handle vh;
    
    for (int i = 0; i < 3; ++ i)
    {
      Point_3 pvertex = point_3 (fh->vertex(i)->info());
      double dist = CGAL::squared_distance (point, pvertex);
      if (dist < min_dist)
      {
        vh = fh->vertex(i);
        min_dist = dist;
      }
    }

    Point_2 p2 = point_2(idx);
    Vector_2 vec (p2, point_2 (vh->info()));
    Point_2 new_p = p2 + 2. * vec;
    Point_3 p3 (new_p.x(), new_p.y(), point.z());

    double cos_angle;
    double sq_distance;

    Face_handle hint = m_dt.locate(new_p, fh);
    if (m_dt.is_infinite(hint))
      return false;
    
    std::tie (cos_angle, sq_distance) = cos_angle_and_sq_distance(p3, hint);

    if (cos_angle > m_min_cos_angle && sq_distance < m_max_sq_distance)
      return true;
    // else
    return false;
  }

void init_triangulation (unsigned int neighbors, double percent)
  {
    std::vector<size_t> indices;
    indices.reserve(m_points.size());
    m_is_ground.resize (m_points.size());
    for (size_t i = 0; i < size_t(m_points.size()); ++ i)
    {
      indices.push_back(i);
      m_is_ground[i] = false;
    }
       
    // Get initial points as lowest points in each cell of a large grid
    CGAL::Bbox_3 bbox
      = CGAL::bbox_3
      (boost::make_transform_iterator
       (m_points.begin(),
        CGAL::Property_map_to_unary_function<PointMap>(m_point_map)),
       boost::make_transform_iterator
       (m_points.end(),
        CGAL::Property_map_to_unary_function<PointMap>(m_point_map)));

    double dx = bbox.xmax() - bbox.xmin();
    double dy = bbox.ymax() - bbox.ymin();

    std::size_t nb_x = std::size_t(dx / m_max_building_size) + 1;
    std::size_t nb_y = std::size_t(dy / m_max_building_size) + 1;

    CGAL_DEM_CERR << "Using a grid " << nb_x << "*" << nb_y << std::endl;

    std::vector<std::vector<size_t> > grid;
    grid.resize(nb_x);
    for (auto& g : grid)
      g.resize(nb_y, no_element());

    typename std::vector<size_t>::iterator
      last = CGAL::remove_outliers (indices, neighbors,
                                    CGAL::parameters::point_map(index_to_point_map())
                                    .threshold_percent(percent));

    CGAL::Iterator_range<typename std::vector<size_t>::iterator>
      filtered_points (indices.begin(), last);

    for (size_t idx : filtered_points)
    {
      std::size_t x = (point_3(idx).x() - bbox.xmin()) / m_max_building_size;
      std::size_t y = (point_3(idx).y() - bbox.ymin()) / m_max_building_size;
      CGAL_assertion (x < nb_x && y < nb_y);

      if (grid[x][y] == no_element())
        grid[x][y] = idx;
      else if (point_3(idx).z() < point_3(grid[x][y]).z())
        grid[x][y] = idx;
    }

    std::size_t nb = 0;
    for (auto& g : grid)
      for (auto& idx : g)
        if (idx != no_element())
        {
          Vertex_handle vh = m_dt.insert (point_2(idx));
          vh->info() = idx;
          ++ nb;
        }
    CGAL_DEM_CERR << "Added " << nb << " points from grid" << std::endl;

    CGAL_DEM_CERR << "Adding corners of bbox" << std::endl;

    std::array<Point_2, 4> bbox_pts
      = { Point_2 (bbox.xmin() - dx * 0.01, bbox.ymin() - dy * 0.01),
          Point_2 (bbox.xmin() - dx * 0.01, bbox.ymax() + dy * 0.01),
          Point_2 (bbox.xmax() + dx * 0.01, bbox.ymin() - dy * 0.01),
          Point_2 (bbox.xmax() + dx * 0.01, bbox.ymax() + dy * 0.01) };
    
    for (const Point_2& p : bbox_pts)
    {
      Vertex_handle vh = m_dt.insert (p);
      double h = 0.;
      double dist_min = std::numeric_limits<double>::max();
      Vertex_circulator start = m_dt.incident_vertices(vh),
        circ = start;
      do
      {
        if (!m_dt.is_infinite(circ) && circ->info() >= 0)
        {
          double dist = CGAL::squared_distance (p, circ->point());
          if (dist < dist_min)
            h = point_3(circ->info()).z();
        }
      }
      while (++ circ != start);

      m_bbox_pts.push_back (Point_3(p.x(), p.y(), h));
      vh->info() = - m_bbox_pts.size();
    }

    std::vector<Vertex_handle> vertices;
    vertices.resize (m_points.size(), Vertex_handle());
    for (Finite_vertices_iterator it = m_dt.finite_vertices_begin();
         it != m_dt.finite_vertices_end(); ++ it)
    {
      if (it->info() >= 0)
      {
        vertices[it->info()] = it;
        m_is_ground[it->info()] = true;
      }
    }

    Face_handle fh;
    
    for (size_t i = 0; i < size_t(m_points.size()); ++ i)
      if (vertices[i] == Vertex_handle())
      {
        typename DT::Locate_type lt;
        int li;
        fh = m_dt.locate (point_2(i), lt, li,fh);

        if (lt == DT::VERTEX)
        {
          CGAL_DEM_CERR << "Ignore same vertex" << std::endl;
          continue;
        }
        
        CGAL_assertion (!m_dt.is_infinite(fh));
        fh->info().inliers.push_back(i);
      }
  }

  Face_handle insert (size_t idx, Face_handle hint)
  {
    Point_2 p2 = point_2(idx);
    
    std::vector<size_t> inliers;

    typename DT::Locate_type lt;
    int li;
    hint = m_dt.locate (p2, lt, li,hint);

    Vertex_handle vh;
    if (lt == DT::VERTEX)
    {
      for (int i = 0; i < 3; ++ i)
        if (p2 == hint->vertex(i)->point())
        {
          vh = hint->vertex(i);
          break;
        }
      CGAL_assertion (vh != Vertex_handle());
      
      Face_circulator start = m_dt.incident_faces(vh),
        circ = start;
      do
      {
        inliers.reserve (inliers.size() + circ->info().inliers.size());
        for (size_t i : circ->info().inliers)
          if (i != idx)
            inliers.push_back(i);
        circ->info().inliers.clear();
      }
      while (++ circ != start);
    }
    else
    {
      std::vector<Face_handle> conflict;
      m_dt.get_conflicts (p2, std::back_inserter (conflict), hint);

      for (Face_handle fh : conflict)
      {
        inliers.reserve (inliers.size() + fh->info().inliers.size());
        for (size_t i : fh->info().inliers)
          if (i != idx)
            inliers.push_back(i);
        fh->info().inliers.clear();
      }

      vh = m_dt.insert (p2, hint);
      vh->info() = idx;
      m_is_ground[vh->info()] = true;
    }
    
    hint = m_dt.incident_faces(vh);
    if (m_dt.is_infinite(hint))
    {
      std::ofstream out ("dbg.off");
      out.precision(18);
      out << *this;
      std::ofstream outx ("dbg.xyz");
      outx.precision(18);
      outx << point_3(idx) << std::endl;
    }
    
    CGAL_assertion (!m_dt.is_infinite(hint));

    Face_circulator start = m_dt.incident_faces(vh),
      circ = start;
    do
    {
      CGAL_assertion (circ->info().inliers.empty());
    }
    while (++ circ != start);

#if 0
    // Test failing -> now succeeding with path_independant locate
    {
      std::vector<Face_handle> parallel_locations;
      parallel_locations.resize (inliers.size(), Face_handle());
      
      tbb::parallel_for (tbb::blocked_range<std::size_t>(0, inliers.size()),
                         [&](const tbb::blocked_range<std::size_t>& r)
                         {
                           Face_handle fh = hint;
                           for (std::size_t s = r.begin(); s != r.end(); ++ s)
                           {
                             size_t i = inliers[s];
                             fh = path_independant_locate (point_2(i), fh);
                             parallel_locations[s] = fh;
                           }
                         });

      std::vector<Face_handle> sequential_locations;
      sequential_locations.resize (inliers.size(), Face_handle());
      for (std::size_t s = 0; s < inliers.size(); ++ s)
      {
        size_t i = inliers[s];
        hint = path_independant_locate (point_2(i), hint);
        sequential_locations[s] = hint;
      }

      std::cerr << parallel_locations.size() << std::endl;
      CGAL_assertion (parallel_locations == sequential_locations);
    }
#endif

#if 1 // Try TBB -> failing so far
    if (inliers.size() > 100)
    {
      std::vector<Face_handle> locations;
      locations.resize (inliers.size(), Face_handle());
      
      tbb::parallel_for (tbb::blocked_range<std::size_t>(0, inliers.size()),
                         [&](const tbb::blocked_range<std::size_t>& r)
                         {
                           Face_handle fh = hint;
                           for (std::size_t s = r.begin(); s != r.end(); ++ s)
                           {
                             size_t i = inliers[s];
                             fh = m_dt.locate (point_2(i), fh);
                             locations[s] = fh;
                             m_is_ground[i] = is_ground(i, fh);
                           }
                         });
      for (std::size_t i = 0; i < locations.size(); ++ i)
      {
        CGAL_assertion (locations[i] != Face_handle());
        locations[i]->info().inliers.push_back(inliers[i]);
      }
      hint = locations.back();
    }
    else
#endif
    for (size_t i : inliers)
    {
      hint = m_dt.locate (point_2(i), hint);
      hint->info().inliers.push_back(i);
      m_is_ground[i] = is_ground(i, hint);
    }


    start = m_dt.incident_faces(vh);
    circ = start;
    do
    {
      CGAL::cpp98::random_shuffle (circ->info().inliers.begin(),
                                   circ->info().inliers.end());
    }
    while (++ circ != start);
    return hint;
  }

  
};

} // namespace Classification

} // namespace CGAL
  

#endif // CGAL_CLASSIFICATION_DEM_H
