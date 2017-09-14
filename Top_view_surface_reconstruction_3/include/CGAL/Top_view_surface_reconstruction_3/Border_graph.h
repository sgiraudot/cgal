#ifndef CGAL_TVSR_BORDER_GRAPH
#define CGAL_TVSR_BORDER_GRAPH

#include <CGAL/boost/graph/split_graph_into_polylines.h>

#include <boost/graph/adjacency_list.hpp>

namespace CGAL
{

template <typename GeomTraits>
class Border_graph : public boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS,
                                                  boost::property<boost::vertex_point_t, typename GeomTraits::Point_2> >
{
public:
  
  typedef GeomTraits Kernel;
  typedef typename Kernel::Point_2 Point_2;

  typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS,
                                boost::property<boost::vertex_point_t, Point_2> > Graph;

  typedef boost::graph_traits<Graph> Traits;
 
  typedef typename Traits::vertex_descriptor     vertex_descriptor;
  typedef typename Traits::edge_descriptor       edge_descriptor;

  typedef std::vector<vertex_descriptor> Polyline;
  typedef typename std::vector<Polyline*>::iterator iterator;

  
private:

  std::vector<Polyline*> m_polylines;
  
public:

  Border_graph()
  {
  }

  virtual ~Border_graph()
  {
    clear();
  }

  std::size_t size() const { return m_polylines.size(); }
  std::vector<vertex_descriptor>& operator[] (const std::size_t& i) { return *(m_polylines[i]); }
  const std::vector<vertex_descriptor>& operator[] (const std::size_t& i) const { return *(m_polylines[i]); }
  iterator begin() { return m_polylines.begin(); }
  iterator end() { return m_polylines.end(); }

  void clear()
  {
    for (std::size_t i = 0; i < m_polylines.size(); ++ i)
      delete m_polylines[i];
    m_polylines.clear();
  }

  const Point_2& point (vertex_descriptor v) const { return get (get(boost::vertex_point, *this), v); }
  Point_2& point (vertex_descriptor v) { return get (get(boost::vertex_point, *this), v); }
  
  void start_new_polyline()
  {
    m_polylines.push_back (new std::vector<vertex_descriptor>());
  }
  
  void add_node (vertex_descriptor v)
  {
    m_polylines.back()->push_back (v);
  }

  void end_polyline()
  {
  }

  void split_into_polylines()
  {
    CGAL::split_graph_into_polylines (*this, *this);
  }

  bool is_terminal (const std::vector<vertex_descriptor>& border)
  {
    return (degree(border.front(), *this) == 1 ||
            degree(border.back(), *this) == 1);
  }

  template <typename OutputIterator>
  void filter_small_terminal_borders (double epsilon,
                                      OutputIterator output)
  {
    std::vector<std::pair<double, vertex_descriptor> > terminal_vertices;
    
    BOOST_FOREACH (vertex_descriptor vd, vertices(*this))
    {
      if (degree(vd, *this) != 1)
        continue;
      
      double length = 0.;
      
      vertex_descriptor previous = vd;
      vertex_descriptor current = *(boost::adjacent_vertices(previous, *this).first);
      std::vector<vertex_descriptor> visited;
      visited.push_back (previous);
      
      while (true)
      {
        visited.push_back (current);
        length += std::sqrt (CGAL::squared_distance (point(previous),
                                                     point(current)));

        if (length > epsilon)
          break;

        if (degree (current, *this) != 2)
          break;
        
        vertex_descriptor next = current;
        BOOST_FOREACH (vertex_descriptor avd, boost::adjacent_vertices(current, *this))
          if (avd != previous)
          {
            next = avd;
            break;
          }

        previous = current;
        current = next;
      }

      if (length < epsilon)
        terminal_vertices.push_back (std::make_pair (length, vd));
    }

    std::sort (terminal_vertices.begin(), terminal_vertices.end());
    for (std::size_t i = 0; i < terminal_vertices.size(); ++ i)
    {
      vertex_descriptor vd = terminal_vertices[i].second;
      if (degree(vd, *this) != 1)
        continue;
      
      double length = 0.;
      
      vertex_descriptor previous = vd;
      vertex_descriptor current = *(boost::adjacent_vertices(previous, *this).first);
      std::vector<vertex_descriptor> visited;
      visited.push_back (previous);
      
      while (true)
      {
        visited.push_back (current);
        length += std::sqrt (CGAL::squared_distance (point(previous),
                                                     point(current)));

        if (length > epsilon)
          break;

        if (degree (current, *this) != 2)
          break;
        
        vertex_descriptor next = current;
        BOOST_FOREACH (vertex_descriptor avd, boost::adjacent_vertices(current, *this))
          if (avd != previous)
          {
            next = avd;
            break;
          }

        previous = current;
        current = next;
      }

      if (length < epsilon)
        for (std::size_t j = 0; j < visited.size() - 1; ++ j)
          remove_edge (visited[j], visited[j+1], *this);
    }
    

    std::vector<vertex_descriptor> to_remove;
    BOOST_FOREACH (vertex_descriptor vd, vertices(*this))
      if (degree(vd, *this) == 0)
        to_remove.push_back (vd);

    for (std::size_t i = 0; i < to_remove.size(); ++ i)
    {
      *(output ++) = point(to_remove[i]);
      remove_vertex (to_remove[i], *this);
    }
  }


  void DEBUG_dump_poly(const char* filename)
  {
    std::ofstream f(filename);
    f.precision(18);

    for (std::size_t i = 0; i < m_polylines.size(); ++ i)
    {
      f << m_polylines[i]->size();
      for (typename std::vector<vertex_descriptor>::iterator it = m_polylines[i]->begin();
           it != m_polylines[i]->end(); ++ it)
        f << " " << point(*it) << " 0.";
      f << std::endl;
    }
  }
};


}
  
#endif // CGAL_TVSR_BORDER_GRAPH
