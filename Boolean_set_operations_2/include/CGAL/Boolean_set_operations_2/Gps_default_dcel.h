// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_GPS_DEFAULT_DCEL_H
#define CGAL_GPS_DEFAULT_DCEL_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * This class is the default \dcel{} class used by the General_polygon_set_2
 * and Polygon_set_2} class-templates to represent the undelying internal
 * Arrangement_2 data structure.
 */

#include <CGAL/Arr_default_dcel.h>

#ifdef CGAL_GPS_USE_CC_DCEL
#include <CGAL/Arr_efficient_dcel.h>
#endif

namespace CGAL {

template <class X_monotone_curve_2>
class Gps_halfedge_base
#ifdef CGAL_GPS_USE_CC_DCEL
  : public Arr_eff_halfedge_base<X_monotone_curve_2>
#else
  : public Arr_halfedge_base<X_monotone_curve_2>
#endif
{
  int _flag;
public:

#ifdef CGAL_GPS_USE_CC_DCEL
  typedef Arr_eff_halfedge_base<X_monotone_curve_2> Base;
#else
  typedef Arr_halfedge_base<X_monotone_curve_2> Base;
#endif

  Gps_halfedge_base()
    : Base()
    , _flag(-1)
  {}

  int flag() const {
    return _flag;
  }

  void set_flag(int i) {
    _flag=i;
  }
};

class Gps_face_base
#ifdef CGAL_GPS_USE_CC_DCEL
  : public Arr_eff_face_base
#else
  : public Arr_face_base
#endif
{
protected:
  mutable char m_info;
  enum
  {
    CONTAINED = 1,
    VISITED   = 2
  };
  std::size_t _id;

#ifdef CGAL_GPS_USE_CC_DCEL
  typedef Arr_eff_face_base Base;
#else
  typedef Arr_face_base Base;
#endif


public:
  //Constructor
  Gps_face_base() :
    Base(),
    m_info(0),
    _id(-1)
  {}

   /*! Assign from another face. */
  virtual void assign (const Base& f)
  {
    Base::assign (f);

    const Gps_face_base & ex_f = static_cast<const Gps_face_base&>(f);
    m_info = ex_f.m_info;
  }

  bool contained() const
  {
    return (m_info & CONTAINED) != 0;
  }

  void set_contained(bool b)
  {
    if (b)
      m_info |= CONTAINED;
    else
      m_info &= ~CONTAINED;
  }

  bool visited() const
  {
    return (m_info & VISITED) != 0;
  }

  void set_visited(bool b) const
  {
    if (b)
      m_info |= VISITED;
    else
      m_info &= ~VISITED;
  }

  Base::Outer_ccbs_container&
  _outer_ccbs()
  {
    return this->outer_ccbs;
  }

  Base::Inner_ccbs_container&
  _inner_ccbs()
  {
    return this->inner_ccbs;
  }

  std::size_t id() const {
    return _id;
  }

  bool id_not_set() const {
    return _id==std::size_t(-1);
  }

  void set_id(std::size_t i) {
    _id=i;
  }

  void reset_id()
  {
    _id=std::size_t(-1);
  }
};

template <class Traits_>
class Gps_default_dcel :
#ifdef CGAL_GPS_USE_CC_DCEL
  public Arr_efficient_dcel<Arr_eff_vertex_base<typename Traits_::Point_2>,
                            Gps_halfedge_base<typename Traits_::X_monotone_curve_2>,
                            Gps_face_base >
#else
  public Arr_dcel_base<Arr_vertex_base<typename Traits_::Point_2>,
                       Gps_halfedge_base<typename Traits_::X_monotone_curve_2>,
                       Gps_face_base >
#endif
{
public:
  /*! Default constructor. */
  Gps_default_dcel() {}
};



} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
