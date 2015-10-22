
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

The class `Ascii_progress_tracker` is an implementation of
`Abstract_progress_tracker` which displays the percentage of progress
of the observed algorithm on the standard error output.

\tparam Observed The observed algorithm class. This class must provide
a method `progress()`.

\tparam ProgressBar If `true`, the tracker displays an ASCII progress
bar in addition of the percentage of progress.

\tparam EstimateRemainingTime If `true`, the tracker displays an
estimation of the time remaining before the observed algorithm is
expected to finish.

*/

template < typename Observed,
           bool ProgressBar = false,
           bool EstimateRemainingTime = false >
class Ascii_progress_tracker : public Abstract_progress_tracker<Observed>
{
public:

/// \name Constructor
/// @{

/*!

  \param refresh_iter The current time is is only checked if
  `notify()` has been called more than this number since the last
  update.

  \param refresh_time The progress state in only checked and displayed
  if the last display happened more than this time ago.

  \warning Using `refresh_iter=0` along with `refresh_time=0` means
  there is no limitation on the rate of display. This can badly affect
  performances. More generally, decreasing the values of these
  parameters may impact the performances.

  \param progress_bar_width Total number of characters of the ASCII
  progress bar.

*/
  Ascii_progress_tracker
  (unsigned int refresh_iter = 1000, // Check if it is time to update every 1000 iterations
   time_t refresh_time = 1, // Update if at least 1 second has passed since last update
   unsigned int progress_bar_width = 40); // Useless if ProgressBar is false

/// @}

/// \name Callback
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
  virtual void notify (const Observed* obs);

 /// @}
};

  

} // namespace CGAL

