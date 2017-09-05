#ifndef CGAL_TVSR_LINE_DETECTOR
#define CGAL_TVSR_LINE_DETECTOR

#include <CGAL/Top_view_surface_reconstruction_3/Surface_mesh_on_cdt.h>
#include <CGAL/Top_view_surface_reconstruction_3/Border_graph.h>

#include <CGAL/Eigen_svd.h>

#include <CGAL/linear_least_squares_fitting_2.h>

#include <boost/iterator/iterator_facade.hpp>

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

  class iterator
    : public boost::iterator_facade< iterator,
                                     Point_2,
                                     std::forward_iterator_tag >
  {
    typedef boost::iterator_facade< iterator,
                                    Point_2,
                                    std::forward_iterator_tag > Facade;
  public:
    iterator() : m_idx(0), m_detector(NULL) { }
                 
    iterator(Line_detector* detector, std::size_t idx)
      : m_idx (idx), m_it (detector->m_lines[idx].buffer.begin()), m_detector (detector)
    {
      m_current = &((*m_it)->info().endpoint);
      CGAL_assertion ((*m_it)->info().has_endpoint());
    }
    iterator(Line_detector* detector, std::size_t idx, bool) // end
      : m_idx (idx), m_it (detector->m_lines[idx].buffer.end()), m_detector (detector)
    {
    }
  private:
    
    friend class boost::iterator_core_access;
    void increment()
    {
      CGAL_assertion(m_detector != NULL);

      ++ m_it;
      
      for (; m_it != m_detector->m_lines[m_idx].buffer.end(); ++ m_it)
      {
        if ((*m_it)->info().has_endpoint()
            && (*m_it)->info().endpoint != *m_current)
        {
          m_current = &((*m_it)->info().endpoint);
          return;
        }
      }

      return;
    }
  
    bool equal(const iterator& other) const
    {
      return m_it == other.m_it;
    }

    Point_2& dereference() const { return const_cast<Point_2&>(*m_current); }

    std::size_t m_idx;
    typename std::deque<Face_handle>::iterator m_it;
    Point_2* m_current;
    Line_detector* m_detector;
  };

  friend iterator;
  
private:
  
  SMCDT& m_mesh;
  std::vector<Face_handle> m_buffer;
  std::vector<Line> m_lines;
  
public:

  Line_detector (SMCDT& mesh) : m_mesh (mesh)
  {
    init();
  }

  std::size_t size() const { return m_lines.size(); }
  iterator begin(std::size_t i) { return iterator(this, i); }
  iterator end(std::size_t i) { return iterator(this, i, true); }

  Point_2& source (const Line& l) { return l.buffer.front()->info().endpoint; }
  Point_2& target (const Line& l) { return l.buffer.back()->info().endpoint; }
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

    // region growing on buffer faces
    for (std::size_t n = 0; n < m_buffer.size(); ++ n)
    {
      Face_handle fh = m_buffer[n];

      Line line;

      grow_region (fh, line.buffer, line.support, epsilon);

      if (line.buffer.size() < 10)
        continue;

      line.buffer.front()->info().endpoint = m_mesh.midpoint(line.buffer.front());
      line.buffer.back()->info().endpoint = m_mesh.midpoint(line.buffer.back());

      for (std::size_t i = 0; i < line.buffer.size(); ++ i)
        line.buffer[i]->info().incident_lines.push_back (m_lines.size());

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

      line.buffer.front()->info().endpoint = m_mesh.midpoint(line.buffer.front());
      line.buffer.back()->info().endpoint = m_mesh.midpoint(line.buffer.back());

      for (std::size_t i = 0; i < line.buffer.size(); ++ i)
        line.buffer[i]->info().incident_lines.push_back (m_lines.size());
      
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

    for (std::size_t n = 0; n < m_buffer.size();)
    {
      if (!m_buffer[n]->info().has_endpoint())
      {
        ++ n;
        continue;
      }

      Face_handle fh = m_buffer[n ++];
      std::vector<std::size_t>& incident_lines = fh->info().incident_lines;

      if (incident_lines.size() == 1) // Find if incident to another line
      {
        Line& line = m_lines[incident_lines[0]];
        
        for (std::size_t i = 0; i < 3; ++ i)
          if (m_mesh.is_handled_buffer (fh->neighbor(i)))
          {
            std::vector<std::size_t>& il2 = fh->neighbor(i)->info().incident_lines;
            bool okay = true;
            for (std::size_t j = 0; j < il2.size(); ++ j)
              if (il2[j] == incident_lines[0])
              {
                okay = false;
                break;
              }
            if (!okay)
              continue;

            to_add.push_back (fh->neighbor(i));

            if (line.buffer.front() == fh)
              line.buffer.push_front(fh->neighbor(i));
            else
              line.buffer.push_back(fh->neighbor(i));
            il2.push_back(incident_lines[0]);

            fh->info().erase_endpoint();
            break;
          }
      }
    }

    std::ofstream conn ("connections.xyz");
    conn.precision(18);
    for (std::size_t i = 0; i < to_add.size(); ++ i)
    {
      to_add[i]->info().endpoint = m_mesh.midpoint (to_add[i]);
      conn << m_mesh.midpoint(to_add[i]) << " 0" << std::endl;
    }
  }

  void regularize (double epsilon)
  {
    for (std::size_t n = 0; n < m_buffer.size(); ++ n)
    {
      if (!m_buffer[n]->info().has_endpoint())
        continue;

      Face_handle fh = m_buffer[n];
      Point_2& point = fh->info().endpoint;
      
      std::vector<std::size_t>& incident_lines = fh->info().incident_lines;

      if (incident_lines.size() == 1) // Simple case, just project
        point = regularized_point_degree_1 (point, m_lines[incident_lines[0]]);
      else if (incident_lines.size() == 2) // Intersection of two lines
      {
        Line& l0 = m_lines[incident_lines[0]];
        Line& l1 = m_lines[incident_lines[1]];

        if ((l0.buffer.front() == fh || l0.buffer.back() == fh) &&
            (l1.buffer.front() == fh || l1.buffer.back() == fh)) // Border to border case
        {
          Point_2 candidate = regularized_point_degree_2_border_intersection (l0, l1);

          if (candidate.x() != candidate.x() ||
              CGAL::squared_distance (candidate, m_mesh.triangle (fh)) > epsilon * epsilon)
            point = regularized_point_degree_2_border_barycenter (point, l0, l1);
          else
            point = candidate;
        }
        else
        {
          Line* lfix = &l0;
          Line* lvary = &l1;
          if ((l0.buffer.front() == fh || l0.buffer.back() == fh))
          {
            lfix = &l1;
            lvary = &l0;
          }

          Point_2 candidate = regularized_point_degree_2_internal_intersection (*lfix, *lvary);

          if (candidate.x() != candidate.x() ||
              CGAL::squared_distance (candidate, m_mesh.triangle (fh)) > epsilon * epsilon)
            point = regularized_point_degree_2_internal_barycenter (point, *lfix, *lvary);
          else
            point = candidate;
        }
      }
      else // Intersection of three lines
      {
        Line& l0 = m_lines[incident_lines[0]];
        Line& l1 = m_lines[incident_lines[1]];
        Line& l2 = m_lines[incident_lines[2]];
        
        Point_2 candidate = regularized_point_degree_3_intersection (l0, l1, l2);

        if (candidate.x() != candidate.x() ||
            CGAL::squared_distance (candidate, m_mesh.triangle (fh)) > epsilon * epsilon)
          point = regularized_point_degree_3_barycenter (point, l0, l1, l2);
        else
          point = candidate;
      }
    }

    // Check and repair self intersection of polyline
    std::ofstream file("intersect.xyz");
    file.precision(18);
    for (std::size_t i = 0; i < m_lines.size(); ++ i)
    {
      Line& line = m_lines[i];
      
      Face_handle source = line.buffer.front();
      for (std::size_t j = 1; j < line.buffer.size(); ++ j)
      {
        if (!line.buffer[j]->info().has_endpoint())
          continue;

        Face_handle target = line.buffer[j];
        if (target == source)
          continue;
 
        if (adjacent_lines_intersect (source, target)) // If intersecting, use barycenters (safer)
        {
          file << source->info().endpoint << " 0" << std::endl
               << target->info().endpoint << " 0" << std::endl;

          std::vector<std::size_t>& ilsource = source->info().incident_lines;
          if (ilsource.size() == 2)
            source->info().endpoint = regularized_point_degree_2_border_barycenter (m_mesh.midpoint(source),
                                                                                    m_lines[ilsource[0]],
                                                                                    m_lines[ilsource[1]]);
          else if (ilsource.size() == 3)
            source->info().endpoint = regularized_point_degree_3_barycenter (m_mesh.midpoint(source),
                                                                             m_lines[ilsource[0]],
                                                                             m_lines[ilsource[1]],
                                                                             m_lines[ilsource[2]]);
          std::vector<std::size_t>& iltarget = target->info().incident_lines;
          if (iltarget.size() == 2)
            target->info().endpoint = regularized_point_degree_2_border_barycenter (m_mesh.midpoint(target),
                                                                                    m_lines[iltarget[0]],
                                                                                    m_lines[iltarget[1]]);
          else if (iltarget.size() == 3)
            target->info().endpoint = regularized_point_degree_3_barycenter (m_mesh.midpoint(target),
                                                                             m_lines[iltarget[0]],
                                                                             m_lines[iltarget[1]],
                                                                             m_lines[iltarget[2]]);
          if (adjacent_lines_intersect (source, target)) // If still intersecting, go back to midpoint (guaranteed)
          {
            source->info().endpoint = m_mesh.midpoint(source);
            target->info().endpoint = m_mesh.midpoint(target);
          }
        }
        source = target;
      }
    }
  }

  Point_2 regularized_point_degree_1 (const Point_2& point, Line& line)
  {
    return line.support.projection (point);
  }

  Point_2 regularized_point_degree_2_border_intersection (Line& l0, Line& l1)
  {
    typename cpp11::result_of<typename Kernel::Intersect_3(Line_2, Line_2)>::type
      result = CGAL::intersection(l0.support, l1.support);
    Point_2* inter;
    if (result && (inter = boost::get<Point_2>(&*result)))
      return *inter;
    return Point_2 (std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN());
  }

  Point_2 regularized_point_degree_2_border_barycenter (const Point_2& point, Line& l0, Line& l1)
  {
    Point_2 p0 = l0.support.projection (point);
    Point_2 p1 = l1.support.projection (point);
    return CGAL::barycenter (p0, std::sqrt (squared_length(l0)),
                             p1, std::sqrt (squared_length(l1)));
  }

  Point_2 regularized_point_degree_2_internal_intersection (Line& lfix, Line& lvary)
  {
    typename cpp11::result_of<typename Kernel::Intersect_3(Line_2, Segment_2)>::type
      result = CGAL::intersection(segment(lfix), lvary.support);
    Point_2* inter;
    if (result && (inter = boost::get<Point_2>(&*result)))
      return *inter;
    return Point_2 (std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN());
  }

  Point_2 regularized_point_degree_2_internal_barycenter (const Point_2& point, Line& lfix, Line& lvary)
  {
    Point_2 pvary = lvary.support.projection (point);
    Point_2 proj = lfix.support.projection(pvary);

    Vector_2 p2s (proj, source(lfix));
    Vector_2 p2t (proj, target(lfix));

    if (p2s * p2t < 0) // In segment
      return CGAL::barycenter (proj, std::sqrt (squared_length(lfix)),
                               pvary, std::sqrt (squared_length(lvary)));

    if (CGAL::squared_distance (pvary, source (lfix))
        < CGAL::squared_distance (pvary, target (lfix)))
      return source(lfix);
    else
      return target(lfix);
  }
  
  Point_2 regularized_point_degree_3_intersection (Line& l0, Line& l1, Line& l2)
  {
    CGAL::Eigen_matrix<double> mat(3, 2);
    CGAL::Eigen_vector<double> vec(3);
    mat.set(0,0, l0.support.a()); mat.set(0,1, l0.support.b()); vec.set(0, -l0.support.c());
    mat.set(1,0, l1.support.a()); mat.set(1,1, l1.support.b()); vec.set(1, -l1.support.c());
    mat.set(2,0, l2.support.a()); mat.set(2,1, l2.support.b()); vec.set(2, -l2.support.c());

    CGAL::Eigen_svd::solve(mat, vec);

    return Point_2 (vec(0), vec(1));
  }

  Point_2 regularized_point_degree_3_barycenter (const Point_2& point, Line& l0, Line& l1, Line& l2)
  {
    Point_2 p0 = l0.support.projection (point);
    Point_2 p1 = l1.support.projection (point);
    Point_2 p2 = l2.support.projection (point);

    return CGAL::barycenter (p0, std::sqrt (squared_length(l0)),
                             p1, std::sqrt (squared_length(l1)),
                             p2, std::sqrt (squared_length(l2)));
  }

  bool adjacent_lines_intersect (Face_handle source, Face_handle target)
  {
    std::vector<std::size_t>& ilsource = source->info().incident_lines;
    std::vector<std::size_t>& iltarget = target->info().incident_lines;
    if (ilsource.size() < 2 || iltarget.size() < 2)
      return false;

    std::vector<Segment_2> segsource;
    for (std::size_t i = 0; i < ilsource.size(); ++ i)
      get_adjacent_segment (m_lines[ilsource[i]], source, std::back_inserter(segsource));
    
    std::vector<Segment_2> segtarget;
    for (std::size_t i = 0; i < iltarget.size(); ++ i)
      get_adjacent_segment (m_lines[iltarget[i]], target, std::back_inserter(segtarget));

    for (std::size_t i = 0; i < segsource.size(); ++ i)
      for (std::size_t j = 0; j < segtarget.size(); ++ j)
        if (CGAL::do_intersect (segsource[i], segtarget[j]))
        {
          // Special case = loop, not considered intersection
          if (segsource[i].source() == segtarget[j].source() ||
              segsource[i].source() == segtarget[j].target() ||
              segsource[i].target() == segtarget[j].source() ||
              segsource[i].target() == segtarget[j].target())
            continue;
          
          return true;
        }
    return false;
  }

  template <typename OutputIterator>
  void get_adjacent_segment (Line& line, Face_handle fh, OutputIterator output)
  {
    std::queue<std::pair<std::size_t, bool> > todo;
    if (fh == line.buffer.front())
      todo.push (std::make_pair (0, true));
    else if (fh == line.buffer.back())
      todo.push (std::make_pair (line.buffer.size() - 1, false));
    else
      for (std::size_t i = 0; i < line.buffer.size(); ++ i)
        if (fh == line.buffer[i])
        {
          todo.push (std::make_pair (i, true));
          todo.push (std::make_pair (i, false));
        }

    while (!todo.empty())
    {
      std::size_t current = todo.front().first;
      bool move_forward = todo.front().second;
      todo.pop();

      if (move_forward)
        ++ current;
      else
        -- current;

      if (line.buffer[current]->info().has_endpoint())
        *(output ++) = Segment_2 (fh->info().endpoint, line.buffer[current]->info().endpoint);
      else
        todo.push (std::make_pair (current, move_forward));
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
    return fh->info().incident_lines.empty();
  }

  bool face_available (Face_handle fh, Face_handle previous)
  {
    const std::vector<std::size_t>& list_fh = fh->info().incident_lines;
    if (list_fh.empty())
      return true;
    
    const std::vector<std::size_t>& list_prev = previous->info().incident_lines;

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
      if (m_buffer[n]->info().incident_lines.size() == 1)
      {
        srand(m_buffer[n]->info().incident_lines[0]);
        red = 64 + rand() % 128;
        green = 64 + rand() % 128;
        blue = 64 + rand() % 128;
      }
      else if (m_buffer[n]->info().incident_lines.size() > 1)
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

    if (project)
      for (std::size_t i = 0; i < m_lines.size(); ++ i)
        file << "2 " << m_lines[i].support.projection(source(m_lines[i])) << " 0 "
             << m_lines[i].support.projection(target(m_lines[i])) << " 0" << std::endl;
    else
    {
      for (std::size_t i = 0; i < m_lines.size(); ++ i)
      {
        iterator previous = begin(i);
        iterator it = previous;
        ++ it;
        for (; it != end(i); ++ it)
        {
          file << "2 " << *previous << " 0 " << *it << " 0" << std::endl;
          previous = it;
        }
      }
    }
  }
  
};

} // namespace CGAL


#endif
