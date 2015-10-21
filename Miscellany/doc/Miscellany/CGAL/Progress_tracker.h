
namespace CGAL {

/*!
\ingroup PkgProfilingTools

The class `Abstract_progress_tracker` is an abstract class for keeping
track of the progresses of algorithm. 

\tparam Observed The observed algorithm class. This class must provide
a method `progress()`.

*/


template <typename Observed>
class Abstract_progress_tracker
{
public:

/// \name Abstract callback
/// @{

/*!
  This method is called by the `Observed` object when its progress has
  changed.

  \warning The `Observed` object calls this method whenever its
  progress has changed. It is the responsability of the
  `Progress_tracker` to make sure it does not call the method
  `progress()` too often (which can badly impact the performances if
  `progress()` requires heavy computation).
*/
  virtual void notify (const Observed* obs) = 0;

/// @} 
};


/*!
\ingroup PkgProfilingTools

The class `Simple_progress_tracker` is an implementation of
`Abstract_progress_tracker` which displays the progress of the
observed algorithm on the standard error output.

\tparam Observed The observed algorithm class. This class must provide
a method `progress()`.

*/

template <typename Observed>
class Simple_progress_tracker : public Abstract_progress_tracker<Observed>
{
public:

  Simple_progress_tracker (time_t refresh_time = 1); // Default = update every 1 second
  virtual void notify (const Observed* obs);
};

/*!
\ingroup PkgProfilingTools

The class `Ascii_bar_progress_tracker` is an implementation of
`Abstract_progress_tracker` which displays the progress of the
observed algorithm on the standard error output in the form of an
ASCII progress bar.

\tparam Observed The observed algorithm class. This class must provide
a method `progress()`.

*/

  
template <typename Observed>
class Ascii_bar_progress_tracker : public Abstract_progress_tracker<Observed>
{
public:

  Ascii_bar_progress_tracker (time_t refresh_time = 1, // Default = update every 1 second
                              unsigned int width = 50);  // Width of progress bar
  
  virtual void notify (const Observed* obs);

};

/*!
\ingroup PkgProfilingTools

The class `Abstract_remaining_time_progress_tracker` is a
specialization of `Abstract_progress_tracker` which provides
additional methods to estimate the remaining time and to display a
readable time string.

\tparam Observed The observed algorithm class. This class must provide
a method `progress()`.

*/



template <typename Observed>
class Abstract_remaining_time_progress_tracker : public Abstract_progress_tracker<Observed>
{

public:

  virtual void notify (const Observed* obs) = 0;

  time_t remaining_time (time_t time_done,
                         double fraction_done) const;

  template <typename Stream>
  void display_time (Stream& stream, time_t seconds) const;

};

/*!
\ingroup PkgProfilingTools

The class `Simple_remaining_time_progress_tracker` is an
implementation of `Abstract_remaining_time_progress_tracker` which
displays the progress of the observed algorithm and the estimated
remaining time on the standard error output.

\tparam Observed The observed algorithm class. This class must provide
a method `progress()`.

*/

  
template <typename Observed>
class Simple_remaining_time_progress_tracker
  : public Abstract_remaining_time_progress_tracker<Observed>
{
public:

  Simple_remaining_time_progress_tracker (time_t refresh_time = 1); // Default = update every 1 second
  
  virtual void notify (const Observed* obs);
};

/*!
\ingroup PkgProfilingTools

The class `Ascii_bar_with_remaining_time_progress_tracker` is an
implementation of `Abstract_remaining_time_progress_tracker` which
displays the progress of the observed algorithm and the estimated
remaining time on the standard error output in the form of an ASCII
progress bar.

\tparam Observed The observed algorithm class. This class must provide
a method `progress()`.

*/


template <typename Observed>
class Ascii_bar_with_remaining_time_progress_tracker
  : public Abstract_remaining_time_progress_tracker<Observed>
{
public:

  Ascii_bar_with_remaining_time_progress_tracker (time_t refresh_time = 1, // Default = update every 1 second
                                                  unsigned int width = 30);  // Width of progress bar
  virtual void notify (const Observed* obs);

};

  


} // namespace CGAL

