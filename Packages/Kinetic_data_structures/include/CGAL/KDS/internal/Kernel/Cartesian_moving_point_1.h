// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_CARTESIAN_MOVING_POINT_1_H_
#define CGAL_KDS_CARTESIAN_MOVING_POINT_1_H_
#include <CGAL/KDS/basic.h>
#include <iostream>

CGAL_KDS_BEGIN_INTERNAL_NAMESPACE;

template <class Coordinate_t>
class Cartesian_moving_point_1
{
    protected:
        typedef Cartesian_moving_point_1<Coordinate_t> This;
    public:
//! The type for coordinate values
        typedef Coordinate_t Coordinate;

//! What should I do for this
        typedef typename Coordinate::NT NT;

//! initialize it from polys
        Cartesian_moving_point_1(const Coordinate &x) {
            _coord=x;
        }

//! initialize it from a still point
        template <class Static_point>
        Cartesian_moving_point_1(const Static_point &pt) {
            _coord=pt.x();
        }

//! null
        Cartesian_moving_point_1(){}

//! homogeneous x
        const Coordinate &hx() const
        {
            return _coord;
        }

//! homogeneous w
        const Coordinate hw() const
        {
            return Coordinate(1);
        }

//! x
        const Coordinate &x() const
        {
            return _coord;
        }

        template <class SK>
            struct Static_traits
        {
            typedef typename SK::RT Static_type;

            static Static_type to_static(const This &o, const NT &t, const SK &) {
                return Static_type(o.x()(t));
            }
        };

        template <class Converter>
            struct Coordinate_converter
        {
            Coordinate_converter(const Converter &c): c_(c){}
            typedef Cartesian_moving_point_1<typename Converter::argument_type> argument_type;
            typedef Cartesian_moving_point_1<typename Converter::result_type> result_type;

            result_type operator()(const argument_type &i) const
            {
                return result_type(c_(i.x()));
            }

            Converter c_;
        };

//! Reverse the motion, time must be negated also
        template <class Negate>
            This transformed_coordinates(const Negate &n) const
        {
            return This(n(_coord));
        }

        void write(std::ostream &out) const
        {
            out <<"(" << x() << ")";
        }
    protected:
        Coordinate _coord;
};

template <class Stream, class Coordinate>
inline Stream &operator<<(Stream &out,
const Cartesian_moving_point_1<Coordinate> &point)
{
    point.write(out);
    return out;
}


CGAL_KDS_END_INTERNAL_NAMESPACE

CGAL_KDS_BEGIN_NAMESPACE

/*template <>
template <class Coord, class SK>
class To_static< internal::Cartesian_moving_point_2<Coord>, SK>:
  public To_static_base<typename Coord::NT,
            typename internal::Cartesian_moving_point_2<Coord>,
            typename SK::Point_2> {
  typedef To_static_base<typename Coord::NT,
             typename internal::Cartesian_moving_point_2<Coord>,
             typename SK::Point_2>  P;
public:
  To_static(){}
typename P::result_type operator()(const typename P::argument_type &arg) const {
return typename P::result_type(arg.x()(P::time()),
arg.y()(P::time()));
}
};*/
CGAL_KDS_END_NAMESPACE
#endif
