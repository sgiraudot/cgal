
/*!
\ingroup PkgProfilingToolsConcepts
\cgalConcept

`ProgressTracker` is a concept for keeping track of the progress of an
algorithm.

\cgalHasModel `CGAL::Dummy_progress_tracker`
\cgalHasModel `CGAL::Ascii_progress_tracker<ProgressBar,EstimateRemainingTime>`

*/


class ProgressTracker
{
public:

/// \name Callback
/// @{

/*!
  This method is called by the observed object when its progress has
  changed.

  \tparam Observed The observed algorithm class.

  \note The model of `ProgressTracker` used is responsible for
  checking the progress of the `Observed` object using methods
  provided by this object. This methods are specific to the `Observed`
  class.
*/
  template <typename Observed>
  void notify (const Observed& obs);

/// @} 
};


