
/*!
\ingroup PkgProfilingToolsConcepts
\cgalConcept

`ProgressTracker` is a concept for keeping track of the progress of an
algorithm.

\cgalHasModel `CGAL::Dummy_progress_tracker`

\cgalHasModel `CGAL::Ascii_progress_tracker<ProgressBar,EstimateRemainingTime>`

*/


template <typename Observed>
class ProgressTracker
{
public:

/// \name Callback
/// @{

/*!
  This method is called by the observed object when its progress has
  changed.

  \tparam Observed The observed algorithm class.

  \warning The `Observed` object calls this method whenever its
  progress has changed. It is the responsability of the
  `Progress_tracker` to make sure it does not call the method
  `progress()` too often (which can badly impact the performances if
  `progress()` requires heavy computation).
*/
  template <typename Observed>
  void notify (const Observed* obs);

/*!

  This variant directly uses the progress value sent by the observed
  object.

*/
  void notify (double progress);

/// @} 
};


