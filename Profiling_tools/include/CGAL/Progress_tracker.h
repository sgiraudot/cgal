#ifndef CGAL_PROGRESS_TRACKER_H
#define CGAL_PROGRESS_TRACKER_H


namespace CGAL {


template <typename Observed>
class Abstract_progress_tracker
{
protected:
  std::set<Observed*> m_observed;
  
public:

  Abstract_progress_tracker () { }
  virtual ~Abstract_progress_tracker () { }

  void attach (Observed* o)
  {
    m_observed.insert (o);
  }
  void detach (Observed* o)
  {
    m_observed.erase (o);
  }
                
  
  virtual void notify () = 0;

};


template <typename Observed>
class Simple_progress_tracker : public Abstract_progress_tracker<Observed>
{
private:
  unsigned int m_refresh_time;
  unsigned int m_latest;
  
public:

  Simple_progress_tracker (time_t refresh_time = 1) // Default = update every 1 second
    : m_refresh_time (refresh_time)
  {
    m_latest = time (NULL);
  }
  virtual ~Simple_progress_tracker () { }

  void notify ()
  {
    unsigned int current = time (NULL);
    
    if (current < m_latest + m_refresh_time)
      return;

    for (typename std::set<Observed*>::iterator it = this->m_observed.begin ();
         it != this->m_observed.end (); ++ it)
      {
        std::cerr << "\r" << 100. * (*it)->progress () << "% done";
      }

    m_latest = time (NULL);
  }

};


} // namespace CGAL


#endif // CGAL_PROGRESS_TRACKER_H
