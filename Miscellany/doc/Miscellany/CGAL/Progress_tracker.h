
namespace CGAL {


/*!
\ingroup PkgProfilingTools

The class `Ascii_progress_tracker` is a model for the
`ProgressTracker` concept. It displays the percentage of progress of
the observed algorithm on the standard error output. It can also
estimate the remaining time of the algorithm by assuming a constant
time dependency on the progress percentage.

\cgalModels `ProgressTracker`

\tparam Observed The observed algorithm class. 

\tparam ProgressBar If `true`, the tracker displays an ASCII progress
bar in addition of the percentage of progress.

\tparam EstimateRemainingTime If `true`, the tracker displays an
estimation of the time remaining before the observed algorithm is
expected to finish.

*/

template < typename Observed,
           bool ProgressBar = false,
           bool EstimateRemainingTime = false >
class Ascii_progress_tracker
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

};

/*!
\ingroup PkgProfilingTools

The class `Dummy_progress_tracker` is a model for the
`ProgressTracker` concept. It doesn't do anything and is the default
model when tracking progress is not needed.

\cgalModels `ProgressTracker`

\tparam Observed The observed algorithm class. 


*/

template < typename Observed>
class Dummy_progress_tracker
{
public:
};

  

} // namespace CGAL

