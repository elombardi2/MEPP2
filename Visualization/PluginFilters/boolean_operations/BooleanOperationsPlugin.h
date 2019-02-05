// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#if(_MSC_VER >= 1400)
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif

#include "Visualization/Plugins/PluginInterface.h"

#include <QStringList>
#include "Dialogs/DialogBooleanOperations1.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePlugin.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

#include "FEVV/Filters/CGAL/Polyhedron/Boolean_Operations/boolean_operations.hpp"

#include "FEVV/Wrappings/properties.h"

// Polyhedron is used internally
#include "FEVV/Wrappings/properties_polyhedron_3.h"

#ifdef FEVV_USE_CGAL
#include "FEVV/Wrappings/properties_surface_mesh.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"
#endif // FEVV_USE_CGAL
#ifdef FEVV_USE_OPENMESH
#include "FEVV/Wrappings/properties_openmesh.h"
#endif // FEVV_USE_OPENMESH
#ifdef FEVV_USE_AIF
#include "FEVV/Wrappings/properties_aif.h"
#endif // FEVV_USE_AIF
#endif

namespace FEVV {

class BooleanOperationsPlugin : public QObject,
                         public Generic_PluginInterface,
                         public BasePlugin
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "BooleanOperationsPlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  BooleanOperationsPlugin() = default;
  ~BooleanOperationsPlugin() = default;

public:
  void init() override { init(true, 1.0, 1.0, 1.0); }

  void init(bool _forceCompute, double _x, double _y, double _z)
  {
    *value_forceCompute = _forceCompute;

    *value_x = _x;
    *value_y = _y;
    *value_z = _z;
  }

  void reset() override
  {
    init();

    emit resetSignal();
  }

  void addParameters(BaseWindow *_window) override
  {
    window = _window;
    if(window == nullptr || !window->isInit())
    {
      std::cerr << "BaseWindow is null or not initialized." << std::endl;
      return;
    }
  }

  template< typename HalfedgeGraph >
  void process(HalfedgeGraph *_mesh)
  {
    std::cout << "Asking to apply BooleanOperations filter ! " << std::endl;

    auto pm = get(boost::vertex_point, *_mesh);

    //TODO-elo-restore  FEVV::Filters::boolean_union(m1, m2, m_out);
    //TODO-elo  call on 2 meshes
    //TODO-elo  create a mesh for output, how to display?
    //TODO-elo-rm-beg  quick test
    FEVV::Filters::boolean_union(*_mesh, *_mesh, *_mesh);
    //TODO-elo-rm-end  quick test
  }

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    if(*value_forceCompute)
      process(_mesh);

    SimpleViewer< HalfedgeGraph > *viewer =
        dynamic_cast< SimpleViewer< HalfedgeGraph > * >(_adapter->getViewer());
    if(viewer)
      viewer->draw_or_redraw_mesh(_mesh, pmaps_bag, true, false);

    reset();

    viewer->frame();
  }

#ifdef FEVV_USE_OPENMESH
  void apply(BaseAdapterVisu *_adapter,
             MeshOpenMesh *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
#if 0
		//TODO-elo  restore when compil error fixed with OM
    applyHG< MeshOpenMesh >(_adapter, _mesh, pmaps_bag);
#else
    QMessageBox::information(
        0,
        "",
        QObject::tr(
            "Boolean Operations filter is not yet compatible with OpenMesh!"));
#endif
  }
#endif

#ifdef FEVV_USE_CGAL
  void apply(BaseAdapterVisu *_adapter,
             MeshLCC *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshLCC >(_adapter, _mesh, pmaps_bag);
  }

  void apply(BaseAdapterVisu *_adapter,
             MeshSurface *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshSurface >(_adapter, _mesh, pmaps_bag);
  }

  void apply(BaseAdapterVisu *_adapter,
             MeshPolyhedron *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshPolyhedron >(_adapter, _mesh, pmaps_bag);
  }
#endif

#ifdef FEVV_USE_AIF
  void apply(BaseAdapterVisu *_adapter,
             MeshAIF *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
#if 0
		//TODO-elo  restore when num_halfedges() is implemented with AIF
		applyHG<MeshAIF>(_adapter, _mesh, pmaps_bag);
#else
    QMessageBox::information(
        0,
        "",
        QObject::tr(
            "Boolean Operations filter is not yet compatible with AIF!"));
#endif
  }
#endif

  QStringList Generic_plugins() const override
  {
    return QStringList() << "BooleanOperationsPlugin";
  }

  bool Generic_plugin(const QString &plugin) override
  {
    DialogBooleanOperations1 dial1;
    dial1.setParameters(*value_x, *value_y, *value_z);
    if(dial1.exec() == QDialog::Accepted)
    {
      dial1.getParameters(*value_x, *value_y, *value_z);

      SimpleWindow *sw = static_cast< SimpleWindow * >(
          window); // dynamic_cast fails under OS X

      sw->onModificationParam("booleanoperations_qt_p", this);
      sw->onApplyButton();

      return true;
    }

    return false;
  }

signals:
  void resetSignal();

protected:
  bool *value_forceCompute = new bool(false);

  double *value_x = new double(0.0);
  double *value_y = new double(0.0);
  double *value_z = new double(0.0);
};

} // namespace FEVV

