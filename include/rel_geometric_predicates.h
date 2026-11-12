#ifndef REL_GEOMETRIC_PREDICATES_H
#define REL_GEOMETRIC_PREDICATES_H

//! Robust geometric predicates.
//! Each predicate is first evaluated in floating point arithmetic together with a bound of its
//! rounding error. Only if the result is smaller than the bound it is evaluated again exactly,
//! with floating point expansions. Sign of a returned value is therefore always correct.
class RGeometricPredicates
{

    private:

        //! Constructor.
        RGeometricPredicates() {}

    public:

        //! Orientation of four points.
        //! Return positive value if point d lies on the side of plane (a,b,c) which normal (b-a)x(c-a) points to,
        //! negative value if it lies on the other side and zero if all four points are coplanar.
        //! Points are arrays of three coordinates.
        static double orient3d(const double *pa, const double *pb, const double *pc, const double *pd);

        //! Position of point e relative to a sphere passing through points a, b, c and d.
        //! Points a, b, c and d must be positively oriented, orient3d(a,b,c,d) > 0.
        //! Return positive value if point e lies inside of the sphere, negative value if it lies
        //! outside and zero if all five points are cospherical.
        static double inSphere(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe);

};

#endif // REL_GEOMETRIC_PREDICATES_H
