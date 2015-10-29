
/*!
\ingroup PkgProfilingToolsConcepts
\cgalConcept

`ProgressTracker` is a concept for keeping track of the progress of an
algorithm.

\tparam Observed The observed algorithm class.

\cgalHasModel `CGAL::Dummy_progress_tracker<Observed>`

\cgalHasModel `CGAL::Ascii_progress_tracker<Observed,ProgressBar,EstimateRemainingTime>`

*/


template <typename Observed>
class ProgressTracker
{
public:

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
  void notify (const Observed* obs);

/*!

  This variant directly uses the progress value sent by the `Observed`
  object. 

*/
  void notify (double progress);

/// @} 
};


