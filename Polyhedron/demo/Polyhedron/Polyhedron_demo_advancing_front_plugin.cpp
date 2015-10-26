#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include <Scene_polyhedron_item.h>
#include "Kernel_type.h"
#include "Polyhedron_type.h"
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Progress_tracker.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QProgressDialog>

#include "ui_Polyhedron_demo_advancing_front_plugin.h"


template <typename Observed>
class Qt_progress_tracker :
  public CGAL::Abstract_progress_tracker<Observed>,
  public QProgressDialog
{
private:
  unsigned int m_refresh_iter;
  unsigned int m_current_iter;
  time_t m_refresh_time;
  time_t m_latest;

public:

  Qt_progress_tracker (QWidget* parent = 0,
                       Qt::WindowFlags f = 0,
                       unsigned int refresh_iter = 1000,
                       time_t refresh_time = 1)
    : QProgressDialog (parent, f),
      m_refresh_iter (refresh_iter),
      m_refresh_time (refresh_time)
  {
    m_current_iter = 0;
    m_latest = time (NULL);


    this->setOrientation (Qt::Vertical);
    this->setMinimum(0);
    this->setMaximum(100);
    this->setCancelButton (0);

    this->setMinimumDuration (0);

    this->setValue(0);

    std::cerr << "Created" << std::endl;
    

  }
  virtual ~Qt_progress_tracker () { }

  virtual void notify (const Observed* obs)
  {
    if (m_current_iter ++ < m_refresh_iter)
      return;
    m_current_iter = 0;

    time_t current = time (NULL);
    if (current < m_latest + m_refresh_time)
      return;

    double done = obs->progress ();

    this->setValue ((unsigned int)(100. * done));
    QApplication::processEvents();
    std::cerr << done << std::endl;
    
    m_latest = time (NULL);
  }
  


};



struct Perimeter {

  double bound;

  Perimeter(double bound)
    : bound(bound)
  {}

  bool operator()(const Kernel::Point_3& p, const Kernel::Point_3& q, const Kernel::Point_3& r) const
  {
    if(bound == 0){
      return false;
    }
    double d  = sqrt(squared_distance(p,q));
    if(d>bound) return true;
    d += sqrt(squared_distance(p,r)) ;
    if(d>bound) return true;
    d+= sqrt(squared_distance(q,r));
    return d>bound;
  }
};


class Polyhedron_demo_advancing_front_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionAdvancingFrontReconstruction;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {

    actionAdvancingFrontReconstruction = new QAction(tr("Advancing Front reconstruction"), mainWindow);
    actionAdvancingFrontReconstruction->setObjectName("actionAdvancingFrontReconstruction");
    
    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  //! Applicate for Point_sets with normals.
  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionAdvancingFrontReconstruction;
  }

public Q_SLOTS:
  void on_actionAdvancingFrontReconstruction_triggered();
}; // end class Polyhedron_demo_advancing_front_plugin


class Polyhedron_demo_advancing_front_plugin_dialog : public QDialog, private Ui::AdvancingFrontDialog
{
  Q_OBJECT
  public:
    Polyhedron_demo_advancing_front_plugin_dialog(QWidget* /*parent*/ = 0)
    {
      setupUi(this);
      
    }

    double trianglePerimeter() const { return m_inputPerimeter->value(); }
};

void Polyhedron_demo_advancing_front_plugin::on_actionAdvancingFrontReconstruction_triggered()
{
  typedef CGAL::Advancing_front_surface_reconstruction_vertex_base_3<Kernel> LVb;
  typedef CGAL::Advancing_front_surface_reconstruction_cell_base_3<Kernel> LCb;

  typedef CGAL::Triangulation_data_structure_3<LVb,LCb> Tds;
  typedef CGAL::Delaunay_triangulation_3<Kernel,Tds> Triangulation_3;

  typedef CGAL::Advancing_front_surface_reconstruction<Triangulation_3, Perimeter> Reconstruction;

  
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* point_set_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(point_set_item)
  {
    // Gets point set
    Point_set* points = point_set_item->point_set();
    if(!points) return;

    // Gets options
    Polyhedron_demo_advancing_front_plugin_dialog dialog;
    if(!dialog.exec())
      return;
    const double sm_perimeter     = dialog.trianglePerimeter();


    QApplication::setOverrideCursor(Qt::WaitCursor);

    // Add polyhedron to scene
    
    // Reconstruct point set as a polyhedron
    Scene_polyhedron_item* new_item = new Scene_polyhedron_item(Polyhedron());
    Polyhedron& P = * const_cast<Polyhedron*>(new_item->polyhedron());
    Perimeter filter(sm_perimeter);

    Triangulation_3 dt (points->begin(), points->end());

    Reconstruction reconstruction (dt, filter);

    Qt_progress_tracker<Reconstruction>* tracker = new Qt_progress_tracker<Reconstruction> ();
    
    tracker->show ();
    
    reconstruction.run (5., 0.52, tracker);

    delete tracker;
    
    CGAL::AFSR::construct_polyhedron (P, reconstruction);

    new_item->setName(tr("%1 Advancing Front (%2)")
                      .arg(point_set_item->name())
                      .arg(sm_perimeter));
    new_item->setColor(Qt::red);
    scene->addItem(new_item);
    
    // Hide point set
    point_set_item->setVisible(false);
    scene->itemChanged(index);
    

    QApplication::restoreOverrideCursor();
  
  }
}

#include "Polyhedron_demo_advancing_front_plugin.moc"
