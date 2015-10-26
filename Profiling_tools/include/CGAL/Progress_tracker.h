#ifndef CGAL_PROGRESS_TRACKER_H
#define CGAL_PROGRESS_TRACKER_H


namespace CGAL {


template <typename Observed>
class Abstract_progress_tracker
{
public:
  virtual void notify (const Observed* obs) = 0;
};


template < typename Observed,
           bool ProgressBar = false,
           bool EstimateRemainingTime = false >
class Ascii_progress_tracker : public Abstract_progress_tracker<Observed>
{
private:

  unsigned int m_refresh_iter;
  unsigned int m_current_iter;

  time_t m_refresh_time;
  time_t m_latest;
  time_t m_starting_time;

  unsigned int m_progress_bar_width;
  
public:

  Ascii_progress_tracker
  (unsigned int refresh_iter = 1000, // Check if it is time to update every 1000 iterations
   time_t refresh_time = 1, // Update if at least 1 second has passed since last update
   unsigned int progress_bar_width = 40) // Useless if ProgressBar is false
    : m_refresh_iter (refresh_iter),
      m_refresh_time (refresh_time),
      m_progress_bar_width (progress_bar_width)
  {
    m_current_iter = 0;
    m_latest = time (NULL);
    m_starting_time = m_latest;
  }
  virtual ~Ascii_progress_tracker () { }

  
  virtual void notify (const Observed* obs)
  {
    if (m_current_iter ++ < m_refresh_iter)
      return;
    m_current_iter = 0;

    time_t current = time (NULL);
    if (current < m_latest + m_refresh_time)
      return;

    double done = obs->progress ();

    std::cerr << "\r";
    
    if (ProgressBar)
      display_progress_bar (done);

    std::cerr << (unsigned int)(100. * done) << "%";

    if (EstimateRemainingTime)
      display_remaining_time (done, current);

    m_latest = time (NULL);
  }

protected:

  inline void display_progress_bar (const double& done) const
  {
    unsigned int bars_filled = (unsigned int)(done * m_progress_bar_width);

    // Erase line
    for (unsigned int i = 0; i < 45 + m_progress_bar_width; ++ i)
      std::cerr << " ";
    
    std::cerr << "\r[";
    unsigned int i = 0;
    for (; i < bars_filled; ++ i)
      std::cerr << "=";
    for (; i < m_progress_bar_width; ++ i)
      std::cerr << " ";

    std::cerr << "] ";

  }

  inline void display_remaining_time (const double& done, const time_t& current) const
  {
    std::cerr << " (";
    display_time (std::cerr, remaining_time (current - m_starting_time, done));
    std::cerr << " remaining)";
  }

  inline time_t remaining_time (const time_t& time_done,
                                const double& fraction_done) const
  {
    return static_cast<time_t>((1. - fraction_done) * time_done / fraction_done);
  }

  template <typename Stream>
  inline void display_time (Stream& stream, time_t seconds) const
  {
    if (seconds > 3600)
      {
	int hours = seconds / 3600;
	stream << hours << "h ";
	seconds -= hours * 3600.;
      }
    if (seconds > 60)
      {
	int minutes = seconds / 60;
	stream << minutes << "min ";
	seconds -= minutes * 60.;
      }
    stream << seconds << "sec";

  }
  
};



  


} // namespace CGAL


#endif // CGAL_PROGRESS_TRACKER_H
