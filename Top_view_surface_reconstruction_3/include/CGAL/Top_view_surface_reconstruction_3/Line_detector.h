#ifndef CGAL_TVSR_LINE_DETECTOR
#define CGAL_TVSR_LINE_DETECTOR

#include <CGAL/Top_view_surface_reconstruction_3/Surface_mesh_on_cdt.h>
#include <CGAL/Top_view_surface_reconstruction_3/Border_graph.h>

#include <CGAL/Eigen_svd.h>

#include <CGAL/linear_least_squares_fitting_2.h>

namespace CGAL
{

template <typename GeomTraits>
class Line_detector
{
public:
  typedef GeomTraits Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Line_3 Line_3;
  typedef typename Kernel::Line_2 Line_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Triangle_2 Triangle_2;
  typedef typename Kernel::Vector_2 Vector_2;
  typedef typename Kernel::Direction_2 Direction_2;

  typedef Surface_mesh_on_cdt<Kernel> SMCDT;
  typedef Border_graph<Kernel> Graph;
  typedef typename SMCDT::Face_handle Face_handle;
  typedef typename SMCDT::Face_index Face_index;
  typedef typename SMCDT::Finite_faces_iterator Finite_faces_iterator;
  typedef typename SMCDT::Line_face_circulator Line_face_circulator;

  struct Sort_by_planarity
  {
    Line_detector& detector;
    double epsilon;

    Sort_by_planarity (Line_detector& detector, double epsilon)
      : detector (detector), epsilon (epsilon) { }

    bool operator() (const Face_handle& a, const Face_handle& b) const
    {
      std::size_t deg_a = detector.degree(a);
      std::size_t deg_b = detector.degree(b);

      if (deg_a < 3 && deg_b == 3)
        return true;
      if (deg_a == 3 && deg_b < 3)
        return false;
      if (deg_a == 3 && deg_b == 3)
        return a < b;
         
      double width_a = width(a);
      double width_b = width(b);
      if (width_a < epsilon && width_b > epsilon)
        return true;
      if (width_a > epsilon && width_b < epsilon)
        return false;
      if (width_a > epsilon && width_b > epsilon)
        return width_a < width_b;

      double dev_a, dev_b;
      std::size_t nb_a, nb_b;
      boost::tie (dev_a, nb_a) = deviation(a);
      boost::tie (dev_b, nb_b) = deviation(b);
      if (nb_a != nb_b)
        return nb_a > nb_b;
      
      return dev_a < dev_b;
    }

    std::pair<std::size_t, double> deviation (const Face_handle& fh) const
    {
      std::vector<Triangle_2> support;
      support.push_back (detector.m_mesh.triangle(fh));
      for (std::size_t i = 0; i < 3; ++ i)
        if (detector.m_mesh.is_handled_buffer(fh->neighbor(i)))
        {
          std::vector<Face_handle> neighbor;
          detector.get_n_next_faces (fh->neighbor(i), fh, 5,
                                     std::back_inserter (neighbor));
          for (std::size_t n = 0; n < neighbor.size(); ++ n)
            support.push_back (detector.m_mesh.triangle(neighbor[n]));
        }
      Line_2 line;
      Point_2 centroid;
      return std::make_pair (support.size(),
                             CGAL::linear_least_squares_fitting_2 (support.begin(), support.end(),
                                                                   line, centroid, CGAL::Dimension_tag<2>()));
    }

    
  };

  friend Sort_by_planarity;

  struct Sort_by_map
  {
    std::map<Face_handle, double>& scores;
    Sort_by_map (std::map<Face_handle, double>& scores) : scores (scores) { }

    bool operator() (const Face_handle& a, const Face_handle& b) const
    {
      typename std::map<Face_handle, double>::iterator it_a, it_b;
      it_a = scores.find(a);
      it_b = scores.find(b);

      double score_a = std::numeric_limits<double>::max();
      double score_b = std::numeric_limits<double>::max();

      if (it_a == scores.end() && it_b == scores.end())
        return a < b;

      if (it_a != scores.end())
        score_a = it_a->second;
      if (it_b != scores.end())
        score_b = it_b->second;

      return score_a < score_b;
    }

  };

  
  struct Line
  {
    Line_2 support;
    std::deque<Face_handle> buffer;

    bool is_valid() const
    {
      return (source.x() == source.x() &&
              source.y() == source.y() &&
              target.x() == target.x() &&
              target.y() == target.y());
    }

    Segment_2 segment() const
    {
      return Segment_2 (source, target);
    }
  };

  struct Face_to_grow_on
  {
    double distance_to_support;
    Face_handle current;
    Face_handle previous;
    bool push_back;

    Face_to_grow_on (double distance_to_support,
                     Face_handle current,
                     Face_handle previous,
                     bool push_back)
      : distance_to_support (distance_to_support)
      , current (current)
      , previous (previous)
      , push_back (push_back)
    { }

    bool operator< (const Face_to_grow_on& other) const
    {
      if (this->distance_to_support != other.distance_to_support)
        return this->distance_to_support < other.distance_to_support;
      return this->current < other.current;
    }
  };

  
private:
  
  SMCDT& m_mesh;
  std::vector<Face_handle> m_buffer;
  std::vector<Line> m_lines;
  std::map<Face_handle, std::vector<std::size_t> > m_map_f2l;
  std::map<Face_handle, Point_2> m_endpoints;
  
public:

  Line_detector (SMCDT& mesh) : m_mesh (mesh)
  {
    init();
  }

  Point_2& source (const Line& l) { return m_endpoints[l.buffer.front()]; }
  Point_2& target (const Line& l) { return m_endpoints[l.buffer.back()]; }
  Segment_2 segment (const Line& l) { return Segment_2 (source(l), target(l)); }
  double squared_length (const Line& l)
  {
    return CGAL::squared_distance (source(l), target(l));
  }
  
  void init()
  {
    for (Finite_faces_iterator it = m_mesh.finite_faces_begin();
         it != m_mesh.finite_faces_end(); ++ it)
      if (m_mesh.is_handled_buffer(it))
        m_buffer.push_back (it);
  }


  void sort_candidates_with_graph (Graph& graph)
  {
    std::map<Face_handle, double> scores;
    
    Face_handle fbegin = Face_handle ();
    Face_handle fend = Face_handle ();
    for (std::size_t i = 0; i < graph.size(); ++ i)
      for (std::size_t j = 0; j < graph[i].size() - 1; ++ j)
      {
        Point_2 p = graph.point (graph[i][j]);
        Point_2 q = graph.point (graph[i][j+1]);

//        double length = std::sqrt (CGAL::squared_distance (p, q));
        
        fbegin = m_mesh.locate (p, fend);
        fend = m_mesh.locate (q, fbegin);

        if (m_mesh.is_handled_buffer (fbegin))
          scores.insert (std::make_pair (fbegin, std::numeric_limits<double>::max()));
        if (m_mesh.is_handled_buffer (fend))
          scores.insert (std::make_pair (fend, std::numeric_limits<double>::max()));
        
        Line_face_circulator it = m_mesh.line_walk(p, q, fbegin);
        ++ it;
        for (; it != fend; ++ it)
          if (m_mesh.is_handled_buffer (it))
          {
            // scores[it] = width(it) / length;
            Point_2 m = m_mesh.midpoint (it);
            double d = (std::min) (CGAL::squared_distance (p, m),
                                   CGAL::squared_distance (q, m));
            scores[it] = width(it) / d;
          }
      }

    std::sort (m_buffer.begin(), m_buffer.end(), Sort_by_map(scores));

  }
                                   


  void run (double epsilon)
  {
    // Sorting faces to pick most relevant first
//    std::sort (m_buffer.begin(), m_buffer.end(), Sort_by_planarity (*this, epsilon));
//    std::random_shuffle (m_buffer.begin(), m_buffer.end());

    // Lines
    m_lines.clear();
    m_map_f2l.clear();
    for (std::size_t n = 0; n < m_buffer.size(); ++ n)
      m_map_f2l.insert (std::make_pair (m_buffer[n], std::vector<std::size_t>()));

    // region growing on buffer faces
    for (std::size_t n = 0; n < m_buffer.size(); ++ n)
    {
      Face_handle fh = m_buffer[n];

      Line line;

      grow_region (fh, line.buffer, line.support, epsilon);

      if (line.buffer.size() < 10)
        continue;

      m_endpoints[line.buffer.front()] = m_mesh.midpoint(line.buffer.front());
      m_endpoints[line.buffer.back()] = m_mesh.midpoint(line.buffer.back());

      for (std::size_t i = 0; i < line.buffer.size(); ++ i)
        m_map_f2l[line.buffer[i]].push_back (m_lines.size());
      m_lines.push_back (line);
      
    }

    TOP_VIEW_CERR << "  " << m_lines.size() << " line(s) detected" << std::endl;
    
    // region growing without tolerance on remaining faces
    for (std::size_t n = 0; n < m_buffer.size(); ++ n)
    {
      Face_handle fh = m_buffer[n];

      Line line;

      grow_region (fh, line.buffer, line.support,
                   std::numeric_limits<double>::max());

      if (line.buffer.size() < 2)
        continue;

      m_endpoints[line.buffer.front()] = m_mesh.midpoint(line.buffer.front());
      m_endpoints[line.buffer.back()] = m_mesh.midpoint(line.buffer.back());

      for (std::size_t i = 0; i < line.buffer.size(); ++ i)
        m_map_f2l[line.buffer[i]].push_back (m_lines.size());
      m_lines.push_back (line);
      
    }
    
    TOP_VIEW_CERR << "  " << m_lines.size() << " line(s) detected after second pass" << std::endl;

    fix_connections();

    DEBUG_dump_polyline("detected_no_regularization.polylines.txt", true);
    
    regularize (epsilon);
  }

  void grow_region (Face_handle fh, std::deque<Face_handle>& faces, Line_2& line, double epsilon)
  {
    if (degree(fh) == 3) // do not originate region from 3
      return;
      
    if (!face_available (fh))
      return;

    Point_2 centroid = m_mesh.midpoint (fh);
    line = Line_2 (centroid, support_vector (fh));
      
    std::vector<Triangle_2> triangles;
    faces.push_back (fh);
    triangles.push_back (m_mesh.triangle (fh));

    std::set<Face_to_grow_on> todo;
    for (std::size_t i = 0; i < 3; ++ i)
      if (m_mesh.is_handled_buffer (fh->neighbor(i)))
        todo.insert (Face_to_grow_on (0., fh->neighbor(i), fh, (todo.size() == 0)));

    CGAL_assertion (todo.size() == 1 || todo.size() == 2);
      
    std::size_t iterations = 0;

    std::set<Face_handle> done;
      
    while (!todo.empty())
    {
      Face_handle current = todo.begin()->current;
      Face_handle previous = todo.begin()->previous;
      bool push_back = todo.begin()->push_back;
      
      todo.erase(todo.begin());

      if (!done.insert (current).second)
        continue;

      Triangle_2 triangle = m_mesh.triangle(current);
      double dist = CGAL::squared_distance (triangle, line);
      if (!face_available (current, previous) ||
          dist > epsilon * epsilon)
        continue;

      if (push_back)
        faces.push_back (current);
      else
        faces.push_front (current);

      std::size_t deg = degree (current);
//      double w = width (current);
      
      if (deg < 3)
        triangles.push_back (triangle);

      // if ((iterations < 10)
      //     || (iterations < 50 && iterations%10 == 0)
      //     || (iterations > 50 && iterations%500 = 0)
      CGAL::linear_least_squares_fitting_2
        (triangles.begin(), triangles.end(),
         line, centroid, CGAL::Dimension_tag<2>());


      if (deg == 1) // End of region reached
        continue;

      if (deg == 2) // Standard case
        for (std::size_t i = 0; i < 3; ++ i)
        {
          if (m_mesh.is_handled_buffer (current->neighbor(i)) && current->neighbor(i) != previous)
          {
            double dist = CGAL::squared_distance (m_mesh.triangle(current->neighbor(i)), line);
            todo.insert (Face_to_grow_on (dist, current->neighbor(i), current, push_back));
          }
        }
      else // Degree 3
      {
        int first_buffer = -1;
        int chosen = -1;
        for (std::size_t i = 0; i < 3; ++ i)
        {
          if (m_mesh.is_handled_buffer (current->neighbor(i)) && current->neighbor(i) != previous)
          {
            if (first_buffer ==  -1)
              first_buffer = int(i);
            else
            {
              std::vector<Face_handle> next_path;
              get_n_next_faces (current->neighbor(first_buffer), current, 5,
                                std::back_inserter (next_path));
                                
              Triangle_2 tr_first = m_mesh.triangle (next_path.back());

              next_path.clear();
              get_n_next_faces (current->neighbor(i), current, 5,
                                std::back_inserter (next_path));
              Triangle_2 tr_i = m_mesh.triangle (next_path.back());
                
              if (CGAL::squared_distance (tr_first, line) < CGAL::squared_distance (tr_i, line))
                chosen = first_buffer;
              else
                chosen = i;
              break;
            }
          }
        }
        CGAL_assertion (chosen != -1);
        double dist = CGAL::squared_distance (m_mesh.triangle(current->neighbor(chosen)), line);
        todo.insert (Face_to_grow_on (dist, current->neighbor(chosen), current, push_back));
      }
        
      ++ iterations;
    }
  }

  void fix_connections()
  {
    std::vector<Face_handle> to_add;
    typename std::map<Face_handle, Point_2>::iterator it = m_endpoints.begin();
    while (it != m_endpoints.end())
    {
      typename std::map<Face_handle, Point_2>::iterator current = it ++;
      
      Face_handle fh = current->first;
      std::vector<std::size_t>& incident_lines = m_map_f2l[fh];

      if (incident_lines.size() == 1) // Find if incident to another line
      {
        Line& line = m_lines[incident_lines[0]];
        
        for (std::size_t i = 0; i < 3; ++ i)
          if (m_mesh.is_handled_buffer (fh->neighbor(i)))
          {
            std::vector<std::size_t>& il2 = m_map_f2l[fh->neighbor(i)];
            bool okay = true;
            for (std::size_t j = 0; j < il2.size(); ++ j)
              if (il2[j] == incident_lines[0])
              {
                okay = false;
                break;
              }
            if (!okay)
              continue;

//            m_endpoints[fh->neighbor(i)] = m_mesh.midpoint (fh->neighbor(i));
            to_add.push_back (fh->neighbor(i));

            if (line.buffer.front() == fh)
              line.buffer.push_front(fh->neighbor(i));
            else
              line.buffer.push_back(fh->neighbor(i));
            il2.push_back(incident_lines[0]);

            m_endpoints.erase (current);
            break;
          }
      }
    }

    std::ofstream conn ("connections.xyz");
    conn.precision(18);
    for (std::size_t i = 0; i < to_add.size(); ++ i)
    {
      m_endpoints[to_add[i]] = m_mesh.midpoint (to_add[i]);
      conn << m_mesh.midpoint(to_add[i]) << " 0" << std::endl;
    }
  }

  void regularize (double epsilon)
  {
    std::ofstream deg1 ("deg1.xyz");
    deg1.precision(18);
    std::ofstream deg2_border_ls ("deg2_border_ls.xyz");
    deg2_border_ls.precision(18);
    std::ofstream deg2_border_mid ("deg2_border_mid.xyz");
    deg2_border_mid.precision(18);
    std::ofstream deg3_border_ls ("deg3_border_ls.xyz");
    deg3_border_ls.precision(18);
    std::ofstream deg3_border_mid ("deg3_border_mid.xyz");
    deg3_border_mid.precision(18);

    std::set<Face_handle> done;

    for (typename std::map<Face_handle, Point_2>::iterator it = m_endpoints.begin();
         it != m_endpoints.end(); ++ it)
    {
      Face_handle fh = it->first;
      Point_2& point = it->second;
      
      std::vector<std::size_t>& incident_lines = m_map_f2l[fh];

      if (incident_lines.size() == 1) // Simple case, just project
      {
        Line& line = m_lines[incident_lines[0]];
        point = line.support.projection (point);
        done.insert (fh);
        deg1 << point << " 0" << std::endl;
      }
      else if (incident_lines.size() == 2) // Intersection of two lines
      {
        Line& l0 = m_lines[incident_lines[0]];
        Line& l1 = m_lines[incident_lines[1]];

        if ((l0.buffer.front() == fh || l0.buffer.back() == fh) &&
            (l1.buffer.front() == fh || l1.buffer.back() == fh)) // Border to border case
        {
          // Try intersecting
          typename cpp11::result_of<typename Kernel::Intersect_3(Line_2, Line_2)>::type
            result = CGAL::intersection(l0.support, l1.support);
          Point_2* inter;
          if (result && (inter = boost::get<Point_2>(&*result))
              && CGAL::squared_distance (*inter, m_mesh.triangle (fh)) < epsilon * epsilon)
          {
            point = *inter;
            done.insert (fh);
            deg2_border_ls << point << " 0" << std::endl;
          }
          // If not good, use weighted midpoint
          else
          {
            Point_2 p0 = l0.support.projection (point);
            Point_2 p1 = l1.support.projection (point);
            point = CGAL::barycenter (p0, std::sqrt (squared_length(l0)),
                                      p1, std::sqrt (squared_length(l1)));
            done.insert (fh);
            deg2_border_mid << point << " 0" << std::endl;
          }
        }
      }
      else // Intersection of three lines
      {
        bool border = true;
        for (std::size_t i = 0; i < 3; ++ i)
          if (m_lines[incident_lines[i]].buffer.front() != fh &&
              m_lines[incident_lines[i]].buffer.back() != fh) // not border
          {
            border = false;
            break;
          }

        if (border)
        {
          Line& l0 = m_lines[incident_lines[0]];
          Line& l1 = m_lines[incident_lines[1]];
          Line& l2 = m_lines[incident_lines[2]];

          CGAL::Eigen_matrix<double> mat(3, 2);
          CGAL::Eigen_vector<double> vec(3);
          mat.set(0,0, l0.support.a()); mat.set(0,1, l0.support.b()); vec.set(0, -l0.support.c());
          mat.set(1,0, l1.support.a()); mat.set(1,1, l1.support.b()); vec.set(1, -l1.support.c());
          mat.set(2,0, l2.support.a()); mat.set(2,1, l2.support.b()); vec.set(2, -l2.support.c());

          CGAL::Eigen_svd::solve(mat, vec);

          Point_2 candidate (vec(0), vec(1));

          if (CGAL::squared_distance (point, m_mesh.triangle (fh)) < epsilon * epsilon)
          {
            point = candidate;
            done.insert (fh);
            deg3_border_ls << point << " 0" << std::endl;
          }
          // If not good, use weighted midpoint            
          else
          {
            Point_2 p0 = l0.support.projection (point);
            Point_2 p1 = l1.support.projection (point);
            Point_2 p2 = l2.support.projection (point);

            point = CGAL::barycenter (p0, std::sqrt (squared_length(l0)),
                                      p1, std::sqrt (squared_length(l1)),
                                      p2, std::sqrt (squared_length(l2)));
            done.insert (fh);
            deg3_border_mid << point << " 0" << std::endl;
          }

        }
      }
    }

    for (std::size_t i = 0; i < m_lines.size(); ++ i)
      if (done.find(m_lines[i].buffer.front()) != done.end() &&
          done.find(m_lines[i].buffer.back()) != done.end ())
        m_lines[i].support = Line_2 (source(m_lines[i]), target(m_lines[i]));

    std::ofstream deg2_inside_proj ("deg2_inside_proj.xyz");
    deg2_inside_proj.precision(18);
    std::ofstream deg2_inside_mid ("deg2_inside_mid.xyz");
    deg2_inside_mid.precision(18);

    for (typename std::map<Face_handle, Point_2>::iterator it = m_endpoints.begin();
         it != m_endpoints.end(); ++ it)
    {
      Face_handle fh = it->first;
      Point_2& point = it->second;

      if (done.find(fh) != done.end())
        continue;

      std::vector<std::size_t>& incident_lines = m_map_f2l[fh];

      if (incident_lines.size() == 2) // Intersection of two lines
      {
        Line& l0 = m_lines[incident_lines[0]];
        Line& l1 = m_lines[incident_lines[1]];

        Line* lfix = &l0;
        Line* lvary = &l1;
        if ((l0.buffer.front() == fh || l0.buffer.back() == fh))
        {
          lfix = &l1;
          lvary = &l0;
        }

        // Try intersecting
        typename cpp11::result_of<typename Kernel::Intersect_3(Line_2, Segment_2)>::type
          result = CGAL::intersection(segment(*lfix), lvary->support);
        Point_2* inter;
        if (result && (inter = boost::get<Point_2>(&*result))
            && CGAL::squared_distance (*inter, m_mesh.triangle (fh)) < epsilon * epsilon)
        {
          point = *inter;
          done.insert (fh);
          deg2_inside_proj << point << " 0" << std::endl;
        }
        // If not good, use closest point
        else
        {
          Point_2 pvary = lvary->support.projection (m_mesh.midpoint (fh));
          Point_2 proj = lfix->support.projection(pvary);

          Vector_2 p2s (proj, source(*lfix));
          Vector_2 p2t (proj, target(*lfix));

          if (p2s * p2t < 0) // In segment
            point = proj;
          else
          {
            point = target(*lfix);
            if (CGAL::squared_distance (pvary, source (*lfix))
                < CGAL::squared_distance (pvary, target (*lfix)))
              point = source(*lfix);
          }
          done.insert (fh);
          deg2_inside_mid << point << " 0" << std::endl;
        }
        update_support (*lvary);
      }
    }
  }

  void update_support (Line& line)
  {
    line.support = Line_2 (source(line), target(line));        

    for (typename std::deque<Face_handle>::iterator it = line.buffer.begin();
         it != line.buffer.end(); ++ it)
    {
      if (*it == line.buffer.front() || *it == line.buffer.back())
        continue;

      typename std::map<Face_handle, Point_2>::iterator found
        = m_endpoints.find (*it);
      if (found == m_endpoints.end())
        continue;

      found->second = line.support.projection (found->second);
    }

  }
  
  Vector_2 support_vector (Face_handle fh)
  {
    int first_buffer = -1;
    for (std::size_t i = 0; i < 3; ++ i)
      if (!m_mesh.is_handled_buffer(fh->neighbor(i)))
      {
        if (first_buffer == -1)
          first_buffer = int(i);
        else // If face of degree 1, use vector from common vertex to midpoint of opposite edge
        {
          const Point_2& p = fh->vertex(3 - (std::size_t(first_buffer) + i))->point();
          const Point_2& p0 = fh->vertex(first_buffer)->point();
          const Point_2& p1 = fh->vertex(i)->point();

          return Vector_2 (p, Point_2 ((p0.x() + p1.x()) / 2.,
                                       (p0.y() + p1.y()) / 2.));
        }
      }

    return Vector_2 (fh->vertex((first_buffer+1)%3)->point(),
                     fh->vertex((first_buffer+2)%3)->point());
  }
  
  std::size_t degree (Face_handle fh)
  {
    std::size_t out = 0;
    for (std::size_t i = 0; i < 3; ++ i)
      if (m_mesh.is_handled_buffer(fh->neighbor(i)))
        ++ out;
    return out;
  }

  double width (const Face_handle& fh) const
  {
    int first_buffer = -1;
    for (std::size_t i = 0; i < 3; ++ i)
      if (!m_mesh.is_handled_buffer(fh->neighbor(i)))
      {
        if (first_buffer == -1)
          first_buffer = int(i);
        else
        {
          return std::sqrt (CGAL::squared_distance (fh->vertex(first_buffer)->point(),
                                                    fh->vertex(i)->point()));
        }
      }
    if (first_buffer == -1)
      return 0.;
    
    return std::sqrt (CGAL::squared_distance (fh->vertex(first_buffer)->point(),
                                              Line_2 (fh->vertex((first_buffer+1)%3)->point(),
                                                      fh->vertex((first_buffer+2)%3)->point())));
  }
  
  bool face_available (Face_handle fh)
  {
    return m_map_f2l[fh].empty();
  }

  bool face_available (Face_handle fh, Face_handle previous)
  {
    const std::vector<std::size_t>& list_fh = m_map_f2l[fh];
    if (list_fh.empty())
      return true;
    
    const std::vector<std::size_t>& list_prev = m_map_f2l[previous];

    for (std::size_t i = 0; i < list_fh.size(); ++ i)
      for (std::size_t j = 0; j < list_prev.size(); ++ j)
        if (list_fh[i] == list_prev[j])
          return false;
    return true;
  }

  template <typename OutputIterator>
  void get_n_next_faces (Face_handle fh, Face_handle previous, std::size_t n,
                         OutputIterator output)
  {
    while (n != 0)
    {
      *(output ++) = fh;

      Face_handle next = Face_handle();
      for (std::size_t i = 0; i < 3; ++ i)
        if (m_mesh.is_unhandled_buffer (fh->neighbor(i))
            && fh->neighbor(i) != previous)
        {
          if (next == Face_handle())
            next = fh->neighbor(i);
          else
            return;
        }
      if (next == Face_handle())
        return;
      previous = fh;
      fh = next;
      -- n;
    }
  }

  void DEBUG_dump_ply(const char* filename)
  {
    std::ofstream f(filename);
    f.precision(18);

    f << "ply" << std::endl
      << "format ascii 1.0" << std::endl
      << "element vertex " << 3 * m_buffer.size() << std::endl
      << "property double x" << std::endl
      << "property double y" << std::endl
      << "property double z" << std::endl
      << "element face " << m_buffer.size() << std::endl
      << "property list uchar int vertex_indices" << std::endl
      << "property uchar red" << std::endl
      << "property uchar green" << std::endl
      << "property uchar blue" << std::endl
      << "end_header" << std::endl;

    for (std::size_t n = 0; n < m_buffer.size(); ++ n)
      for (std::size_t i = 0; i < 3; ++ i)
        f << m_buffer[n]->vertex(i)->point() << " 0" << std::endl;

    for (std::size_t n = 0; n < m_buffer.size(); ++ n)
    {
      f << "3";
      for (std::size_t i = 0; i < 3; ++ i)
        f << " " << 3*n + i;

      int red = 0, green = 0, blue = 0;
      if (m_map_f2l[m_buffer[n]].size() == 1)
      {
        srand(m_map_f2l[m_buffer[n]][0]);
        red = 64 + rand() % 128;
        green = 64 + rand() % 128;
        blue = 64 + rand() % 128;
      }
      else if (m_map_f2l[m_buffer[n]].size() > 1)
      {
        red = 128; green = 128; blue = 128;
      }

      f << " " << red << " " << green << " " << blue << std::endl;
    }
      
  }

  void DEBUG_dump_polyline(const char* filename, bool project = false)
  {
    std::ofstream file (filename);
    file.precision(18);

    for (std::size_t i = 0; i < m_lines.size(); ++ i)
    {
      if (project)
        file << "2 " << m_lines[i].support.projection(source(m_lines[i])) << " 0 "
             << m_lines[i].support.projection(target(m_lines[i])) << " 0" << std::endl;
      else
        file << "2 " << source(m_lines[i]) << " 0 " << target(m_lines[i]) << " 0" << std::endl;
    }
    

  }
  
};

} // namespace CGAL


#endif
