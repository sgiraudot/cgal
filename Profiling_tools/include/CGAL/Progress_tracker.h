#ifndef CGAL_PROGRESS_TRACKER_H
#define CGAL_PROGRESS_TRACKER_H


namespace CGAL {


template <typename Observed>
class Abstract_progress_tracker
{
public:
  virtual void notify (const Observed* obs) = 0;
};


template <typename Observed>
class Simple_progress_tracker : public Abstract_progress_tracker<Observed>
{
private:
  
  time_t m_refresh_time;
  time_t m_latest;
  
public:

  Simple_progress_tracker (time_t refresh_time = 1) // Default = update every 1 second
    : m_refresh_time (refresh_time)
  {
    m_latest = time (NULL);
  }
  virtual ~Simple_progress_tracker () { }

  
  virtual void notify (const Observed* obs)
  {
    time_t current = time (NULL);
    
    if (current < m_latest + m_refresh_time)
      return;
    
    std::cerr << "\r" << (unsigned int)(100. * obs->progress ()) << "% done";

    m_latest = time (NULL);
  }

};

template <typename Observed>
class Ascii_bar_progress_tracker : public Abstract_progress_tracker<Observed>
{
private:

  time_t m_refresh_time;
  time_t m_latest;
  unsigned int m_width;
  
public:

  Ascii_bar_progress_tracker (time_t refresh_time = 1, // Default = update every 1 second
                              unsigned int width = 50)  // Width of progress bar
    : m_refresh_time (refresh_time), m_width (width)
  {
    m_latest = time (NULL);
  }
  virtual ~Ascii_bar_progress_tracker () { }
  
  virtual void notify (const Observed* obs)
  {
    time_t current = time (NULL);
    
    if (current < m_latest + m_refresh_time)
      return;

    double done = obs->progress ();
    unsigned int bars_filled = (unsigned int)(done * m_width);
    
    std::cerr << "\r[";
    unsigned int i = 0;
    for (; i < bars_filled; ++ i)
      std::cerr << "=";
    for (; i < m_width; ++ i)
      std::cerr << " ";
    
    std::cerr << "] " << (unsigned int)(100. * done) << "%";
    

    m_latest = time (NULL);
  }

};


template <typename Observed>
class Abstract_remaining_time_progress_tracker : public Abstract_progress_tracker<Observed>
{
private:
  
public:

  virtual void notify (const Observed* obs) = 0;


  time_t remaining_time (time_t time_done,
                         double fraction_done) const
  {
    return static_cast<time_t>((1. - fraction_done) * time_done / fraction_done);
  }

  template <typename Stream>
  void display_time (Stream& stream, time_t seconds) const
  {
    if (seconds > 3600.)
      {
	int hours = (int)(seconds / 3600.);
	stream << hours << "h ";
	seconds -= hours * 3600.;
      }
    if (seconds > 60.)
      {
	int minutes = (int)(seconds / 60.);
	stream << minutes << "min ";
	seconds -= minutes * 60.;
      }
    stream << seconds << "sec";

  }

};


  
template <typename Observed>
class Simple_remaining_time_progress_tracker
  : public Abstract_remaining_time_progress_tracker<Observed>
{
private:

  time_t m_refresh_time;
  time_t m_latest;
  time_t m_starting_time;
  
public:

  Simple_remaining_time_progress_tracker (time_t refresh_time = 1) // Default = update every 1 second
    : m_refresh_time (refresh_time)
  {
    m_latest = time (NULL);
    m_starting_time = m_latest;
  }
  virtual ~Simple_remaining_time_progress_tracker () { }

  
  virtual void notify (const Observed* obs)
  {
    time_t current = time (NULL);
    
    if (current < m_latest + m_refresh_time)
      return;

    double done = obs->progress ();

    // Erase line
    std::cerr << "\r";
    for (unsigned int i = 0; i < 45; ++ i)
      std::cerr << " ";

    std::cerr << "\r" << (unsigned int)(100. * done) << "% done (";
    this->display_time (std::cerr, this->remaining_time (current - m_starting_time, done));
    std::cerr << " remaining)";

    m_latest = time (NULL);
  }

};


template <typename Observed>
class Ascii_bar_with_remaining_time_progress_tracker
  : public Abstract_remaining_time_progress_tracker<Observed>
{
private:

  time_t m_refresh_time;
  time_t m_latest;
  unsigned int m_width;
  time_t m_starting_time;
  
public:

  Ascii_bar_with_remaining_time_progress_tracker (time_t refresh_time = 1, // Default = update every 1 second
                                                  unsigned int width = 30)  // Width of progress bar
    : m_refresh_time (refresh_time), m_width (width)
  {
    m_latest = time (NULL);
    m_starting_time = m_latest;
  }
  virtual ~Ascii_bar_with_remaining_time_progress_tracker () { }
  
  virtual void notify (const Observed* obs)
  {
    time_t current = time (NULL);
    
    if (current < m_latest + m_refresh_time)
      return;

    double done = obs->progress ();
    unsigned int bars_filled = (unsigned int)(done * m_width);

    // Erase line
    std::cerr << "\r";
    for (unsigned int i = 0; i < 45 + m_width; ++ i)
      std::cerr << " ";
    
    std::cerr << "\r[";
    unsigned int i = 0;
    for (; i < bars_filled; ++ i)
      std::cerr << "=";
    for (; i < m_width; ++ i)
      std::cerr << " ";
    
    std::cerr << "] " << (unsigned int)(100. * done) << "% (";
    this->display_time (std::cerr, this->remaining_time (current - m_starting_time, done));
    std::cerr << " remaining)";
    
    m_latest = time (NULL);
  }

};

  


} // namespace CGAL


#endif // CGAL_PROGRESS_TRACKER_H
