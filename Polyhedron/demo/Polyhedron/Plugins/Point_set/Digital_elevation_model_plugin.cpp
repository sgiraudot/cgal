#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_surface_mesh_item.h"
#include "SMesh_type.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#define CGAL_DEM_VERBOSE
#include <CGAL/Classification/Digital_elevation_model.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Memory_sizer.h>

#include <QMultipleInputDialog.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QSpinBox>
#include "CGAL_double_edit.h"
#include <QMessageBox>

#include "run_with_qprogressdialog.h"

using namespace CGAL::Three;
class Polyhedron_demo_digital_elevation_model_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  QAction* actionDEM;
  
public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    mw = mainWindow;
    actionDEM = new QAction(tr("Digital Elevation Model"), mainWindow);
    actionDEM->setProperty("subMenuName","Point Set Processing");
    actionDEM->setObjectName("actionDEM");
    autoConnectActions();
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionDEM;
  }
  
  //! Applicable if the currently selected item is a
  //! points_with_normal_item.
  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionDEM_triggered();
}; // end Polyhedron_demo_digital_elevation_model_plugin

void Polyhedron_demo_digital_elevation_model_plugin::on_actionDEM_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    // Gets point set
    Point_set* points = item->point_set();
    if(points == NULL)
      return;

    QMultipleInputDialog dialog ("Digital Elevation Model", mw);
    
    DoubleEdit* max_building_size = dialog.add<DoubleEdit> ("Maximum building size:");
    max_building_size->setRange (0, 100000000);
    max_building_size->setValue (50.);
                                 
    DoubleEdit* max_terrain_angle = dialog.add<DoubleEdit> ("Maximum terrain angle:");
    max_terrain_angle->setRange (0, 90);
    max_terrain_angle->setValue (88);
    
    DoubleEdit* max_angle = dialog.add<DoubleEdit> ("Maximum angle:");
    max_angle->setRange (0, 90);
    max_angle->setValue (10.);
    
    DoubleEdit* max_distance = dialog.add<DoubleEdit> ("Maximum distance:");
    max_distance->setRange (0, 100000000);
    max_distance->setValue (1.4);

    DoubleEdit* min_edge_length = dialog.add<DoubleEdit> ("Minimum edge length:");
    min_edge_length->setRange (0, 100000000);
    min_edge_length->setValue (2.0);
    
    QSpinBox* init_outlier_neighbors = dialog.add<QSpinBox> ("Outlier neighbors (initialization):");
    init_outlier_neighbors->setRange (3, 100000);
    init_outlier_neighbors->setValue (12);
    
    DoubleEdit* init_outlier_percent = dialog.add<DoubleEdit> ("Outlier percent (initialization):");
    init_outlier_percent->setRange (0, 100);
    init_outlier_percent->setValue (25.);
    
    Point_set::Property_map<int> label_map;
    bool label_found;
    std::tie (label_map, label_found) = points->property_map<int>("label");
    Point_set::Property_map<unsigned char> classif_map;
    bool classif_found;
    std::tie (classif_map, classif_found) = points->property_map<unsigned char>("classification");

    QSpinBox* label;
    if (label_found || classif_found)
    {
      label = dialog.add<QSpinBox> ("Ground label (for evaluation):");
    }

    if (dialog.exec() != QDialog::Accepted)
      return;


    QApplication::setOverrideCursor(Qt::BusyCursor);
    QApplication::processEvents();
    CGAL::Real_timer task_timer; task_timer.start();

    CGAL::Classification::Digital_elevation_model<Kernel, Point_set, Point_set::Point_map>
      dem (*points, points->point_map(),
           max_building_size->value(),
           max_terrain_angle->value(),
           max_angle->value(),
           max_distance->value(),
           min_edge_length->value(),
           init_outlier_neighbors->value(),
           init_outlier_percent->value(),
           [&](const CGAL::Classification::Digital_elevation_model<Kernel, Point_set, Point_set::Point_map>&)
           { });

    std::stringstream ss;
    ss << dem;
    SMesh *mesh = new SMesh();
    ss >> *mesh;
    Scene_surface_mesh_item* mesh_item = new Scene_surface_mesh_item(mesh);
    mesh_item->setName(tr("%1 (DEM)").arg(item->name()));

    scene->addItem(mesh_item);

    QApplication::restoreOverrideCursor();
  }
}


#include "Digital_elevation_model_plugin.moc"
