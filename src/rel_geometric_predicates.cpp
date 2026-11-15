#include <cmath>
#include <vector>

#include "rel_geometric_predicates.h"

// Exact arithmetic follows J. R. Shewchuk, "Adaptive Precision Floating-Point Arithmetic and
// Fast Robust Geometric Predicates", Discrete & Computational Geometry 18:305-363, 1997.
// A number is held as an expansion, a sum of doubles which do not overlap and are sorted by
// increasing magnitude, so that the last component gives the sign.

namespace
{

    typedef std::vector<double> Expansion;

    //! Half of the distance between 1.0 and the next double, 2^-53.
    const double epsilon = 1.1102230246251565e-16;
    //! Error bound of orient3d evaluated in floating point arithmetic.
    const double orient3dErrorBound = (7.0 + 56.0 * epsilon) * epsilon;
    //! Error bound of inSphere evaluated in floating point arithmetic.
    const double inSphereErrorBound = (16.0 + 224.0 * epsilon) * epsilon;

    inline void fastTwoSum(double a, double b, double &x, double &y)
    {
        x = a + b;
        const double bVirtual = x - a;
        y = b - bVirtual;
    }

    inline void twoSum(double a, double b, double &x, double &y)
    {
        x = a + b;
        const double bVirtual = x - a;
        const double aVirtual = x - bVirtual;
        const double bRoundoff = b - bVirtual;
        const double aRoundoff = a - aVirtual;
        y = aRoundoff + bRoundoff;
    }

    inline void twoDiff(double a, double b, double &x, double &y)
    {
        x = a - b;
        const double bVirtual = a - x;
        const double aVirtual = x + bVirtual;
        const double bRoundoff = bVirtual - b;
        const double aRoundoff = a - aVirtual;
        y = aRoundoff + bRoundoff;
    }

    //! Split a double into two halves of its significand, 2^27 + 1 is the splitter.
    inline void split(double a, double &high, double &low)
    {
        const double c = 134217729.0 * a;
        const double aBig = c - a;
        high = c - aBig;
        low = a - high;
    }

    // Function std::fma() is not used on purpose, as it is not exact with some C libraries.
    inline void twoProduct(double a, double b, double &x, double &y)
    {
        x = a * b;
        double aHigh, aLow, bHigh, bLow;
        split(a, aHigh, aLow);
        split(b, bHigh, bLow);
        const double error1 = x - (aHigh * bHigh);
        const double error2 = error1 - (aLow * bHigh);
        const double error3 = error2 - (aHigh * bLow);
        y = (aLow * bLow) - error3;
    }

    //! Exact difference of two doubles.
    Expansion difference(double a, double b)
    {
        double x, y;
        twoDiff(a, b, x, y);
        if (y != 0.0)
        {
            return Expansion{y, x};
        }
        return Expansion{x};
    }

    //! Exact sum of two expansions.
    Expansion sum(const Expansion &e, const Expansion &f)
    {
        Expansion h;
        h.reserve(e.size() + f.size());

        size_t eIndex = 0;
        size_t fIndex = 0;
        double q, qNew, hh;

        // Components are merged by increasing magnitude.
        auto takeE = [&]() -> bool
        {
            if (eIndex >= e.size())
            {
                return false;
            }
            if (fIndex >= f.size())
            {
                return true;
            }
            const double eNow = e[eIndex];
            const double fNow = f[fIndex];
            return (fNow > eNow) == (fNow > -eNow);
        };

        q = takeE() ? e[eIndex++] : f[fIndex++];

        if (eIndex < e.size() || fIndex < f.size())
        {
            const double next = takeE() ? e[eIndex++] : f[fIndex++];
            fastTwoSum(next, q, qNew, hh);
            q = qNew;
            if (hh != 0.0)
            {
                h.push_back(hh);
            }
            while (eIndex < e.size() || fIndex < f.size())
            {
                const double now = takeE() ? e[eIndex++] : f[fIndex++];
                twoSum(q, now, qNew, hh);
                q = qNew;
                if (hh != 0.0)
                {
                    h.push_back(hh);
                }
            }
        }
        if (q != 0.0 || h.empty())
        {
            h.push_back(q);
        }
        return h;
    }

    //! Exact product of an expansion and a double.
    Expansion scale(const Expansion &e, double b)
    {
        Expansion h;
        h.reserve(2 * e.size());

        double q, hh, product1, product0, s;

        twoProduct(e[0], b, q, hh);
        if (hh != 0.0)
        {
            h.push_back(hh);
        }
        for (size_t i = 1; i < e.size(); i++)
        {
            twoProduct(e[i], b, product1, product0);
            twoSum(q, product0, s, hh);
            if (hh != 0.0)
            {
                h.push_back(hh);
            }
            fastTwoSum(product1, s, q, hh);
            if (hh != 0.0)
            {
                h.push_back(hh);
            }
        }
        if (q != 0.0 || h.empty())
        {
            h.push_back(q);
        }
        return h;
    }

    //! Exact product of two expansions.
    Expansion product(const Expansion &e, const Expansion &f)
    {
        const Expansion &rLong = e.size() >= f.size() ? e : f;
        const Expansion &rShort = e.size() >= f.size() ? f : e;

        Expansion h = scale(rLong, rShort[0]);
        for (size_t i = 1; i < rShort.size(); i++)
        {
            h = sum(h, scale(rLong, rShort[i]));
        }
        return h;
    }

    Expansion negative(const Expansion &e)
    {
        Expansion h(e);
        for (double &rValue : h)
        {
            rValue = -rValue;
        }
        return h;
    }

    inline Expansion subtract(const Expansion &e, const Expansion &f)
    {
        return sum(e, negative(f));
    }

    //! Exact value of a*d - b*c.
    inline Expansion cross(const Expansion &a, const Expansion &b, const Expansion &c, const Expansion &d)
    {
        return subtract(product(a, d), product(b, c));
    }

    //! Exact determinant of rows (a-d), (b-d) and (c-d).
    double orient3dExact(const double *pa, const double *pb, const double *pc, const double *pd)
    {
        const Expansion adx = difference(pa[0], pd[0]);
        const Expansion ady = difference(pa[1], pd[1]);
        const Expansion adz = difference(pa[2], pd[2]);
        const Expansion bdx = difference(pb[0], pd[0]);
        const Expansion bdy = difference(pb[1], pd[1]);
        const Expansion bdz = difference(pb[2], pd[2]);
        const Expansion cdx = difference(pc[0], pd[0]);
        const Expansion cdy = difference(pc[1], pd[1]);
        const Expansion cdz = difference(pc[2], pd[2]);

        const Expansion det = sum(sum(product(adz, cross(bdx, bdy, cdx, cdy)),
                                      product(bdz, cross(cdx, cdy, adx, ady))),
                                  product(cdz, cross(adx, ady, bdx, bdy)));
        return det.back();
    }

    //! Exact determinant of rows (a-e, |a-e|^2), (b-e, |b-e|^2), (c-e, |c-e|^2) and (d-e, |d-e|^2).
    double inSphereExact(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe)
    {
        const Expansion aex = difference(pa[0], pe[0]);
        const Expansion aey = difference(pa[1], pe[1]);
        const Expansion aez = difference(pa[2], pe[2]);
        const Expansion bex = difference(pb[0], pe[0]);
        const Expansion bey = difference(pb[1], pe[1]);
        const Expansion bez = difference(pb[2], pe[2]);
        const Expansion cex = difference(pc[0], pe[0]);
        const Expansion cey = difference(pc[1], pe[1]);
        const Expansion cez = difference(pc[2], pe[2]);
        const Expansion dex = difference(pd[0], pe[0]);
        const Expansion dey = difference(pd[1], pe[1]);
        const Expansion dez = difference(pd[2], pe[2]);

        const Expansion ab = cross(aex, aey, bex, bey);
        const Expansion bc = cross(bex, bey, cex, cey);
        const Expansion cd = cross(cex, cey, dex, dey);
        const Expansion da = cross(dex, dey, aex, aey);
        const Expansion ac = cross(aex, aey, cex, cey);
        const Expansion bd = cross(bex, bey, dex, dey);

        const Expansion abc = sum(subtract(product(aez, bc), product(bez, ac)), product(cez, ab));
        const Expansion bcd = sum(subtract(product(bez, cd), product(cez, bd)), product(dez, bc));
        const Expansion cda = sum(sum(product(cez, da), product(dez, ac)), product(aez, cd));
        const Expansion dab = sum(sum(product(dez, ab), product(aez, bd)), product(bez, da));

        const Expansion aLift = sum(sum(product(aex, aex), product(aey, aey)), product(aez, aez));
        const Expansion bLift = sum(sum(product(bex, bex), product(bey, bey)), product(bez, bez));
        const Expansion cLift = sum(sum(product(cex, cex), product(cey, cey)), product(cez, cez));
        const Expansion dLift = sum(sum(product(dex, dex), product(dey, dey)), product(dez, dez));

        const Expansion det = sum(subtract(product(dLift, abc), product(cLift, dab)),
                                  subtract(product(bLift, cda), product(aLift, bcd)));
        return det.back();
    }

}

double RGeometricPredicates::orient3d(const double *pa, const double *pb, const double *pc, const double *pd)
{
    const double adx = pa[0] - pd[0];
    const double bdx = pb[0] - pd[0];
    const double cdx = pc[0] - pd[0];
    const double ady = pa[1] - pd[1];
    const double bdy = pb[1] - pd[1];
    const double cdy = pc[1] - pd[1];
    const double adz = pa[2] - pd[2];
    const double bdz = pb[2] - pd[2];
    const double cdz = pc[2] - pd[2];

    const double bdxcdy = bdx * cdy;
    const double cdxbdy = cdx * bdy;

    const double cdxady = cdx * ady;
    const double adxcdy = adx * cdy;

    const double adxbdy = adx * bdy;
    const double bdxady = bdx * ady;

    const double det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);

    const double permanent = (std::fabs(bdxcdy) + std::fabs(cdxbdy)) * std::fabs(adz) + (std::fabs(cdxady) + std::fabs(adxcdy)) * std::fabs(bdz) + (std::fabs(adxbdy) + std::fabs(bdxady)) * std::fabs(cdz);
    const double errorBound = orient3dErrorBound * permanent;

    // Determinant of rows (a-d), (b-d) and (c-d) has opposite sign than ((b-a)x(c-a)).(d-a).
    if (det > errorBound || -det > errorBound)
    {
        return -det;
    }
    return -orient3dExact(pa, pb, pc, pd);
}

double RGeometricPredicates::inSphere(const double *pa, const double *pb, const double *pc, const double *pd, const double *pe)
{
    const double aex = pa[0] - pe[0];
    const double bex = pb[0] - pe[0];
    const double cex = pc[0] - pe[0];
    const double dex = pd[0] - pe[0];
    const double aey = pa[1] - pe[1];
    const double bey = pb[1] - pe[1];
    const double cey = pc[1] - pe[1];
    const double dey = pd[1] - pe[1];
    const double aez = pa[2] - pe[2];
    const double bez = pb[2] - pe[2];
    const double cez = pc[2] - pe[2];
    const double dez = pd[2] - pe[2];

    const double aexbey = aex * bey;
    const double bexaey = bex * aey;
    const double ab = aexbey - bexaey;
    const double bexcey = bex * cey;
    const double cexbey = cex * bey;
    const double bc = bexcey - cexbey;
    const double cexdey = cex * dey;
    const double dexcey = dex * cey;
    const double cd = cexdey - dexcey;
    const double dexaey = dex * aey;
    const double aexdey = aex * dey;
    const double da = dexaey - aexdey;

    const double aexcey = aex * cey;
    const double cexaey = cex * aey;
    const double ac = aexcey - cexaey;
    const double bexdey = bex * dey;
    const double dexbey = dex * bey;
    const double bd = bexdey - dexbey;

    const double abc = aez * bc - bez * ac + cez * ab;
    const double bcd = bez * cd - cez * bd + dez * bc;
    const double cda = cez * da + dez * ac + aez * cd;
    const double dab = dez * ab + aez * bd + bez * da;

    const double aLift = aex * aex + aey * aey + aez * aez;
    const double bLift = bex * bex + bey * bey + bez * bez;
    const double cLift = cex * cex + cey * cey + cez * cez;
    const double dLift = dex * dex + dey * dey + dez * dez;

    const double det = (dLift * abc - cLift * dab) + (bLift * cda - aLift * bcd);

    const double aezAbs = std::fabs(aez);
    const double bezAbs = std::fabs(bez);
    const double cezAbs = std::fabs(cez);
    const double dezAbs = std::fabs(dez);
    const double aexbeyAbs = std::fabs(aexbey);
    const double bexaeyAbs = std::fabs(bexaey);
    const double bexceyAbs = std::fabs(bexcey);
    const double cexbeyAbs = std::fabs(cexbey);
    const double cexdeyAbs = std::fabs(cexdey);
    const double dexceyAbs = std::fabs(dexcey);
    const double dexaeyAbs = std::fabs(dexaey);
    const double aexdeyAbs = std::fabs(aexdey);
    const double aexceyAbs = std::fabs(aexcey);
    const double cexaeyAbs = std::fabs(cexaey);
    const double bexdeyAbs = std::fabs(bexdey);
    const double dexbeyAbs = std::fabs(dexbey);

    const double permanent = ((cexdeyAbs + dexceyAbs) * bezAbs + (dexbeyAbs + bexdeyAbs) * cezAbs + (bexceyAbs + cexbeyAbs) * dezAbs) * aLift + ((dexaeyAbs + aexdeyAbs) * cezAbs + (aexceyAbs + cexaeyAbs) * dezAbs + (cexdeyAbs + dexceyAbs) * aezAbs) * bLift + ((aexbeyAbs + bexaeyAbs) * dezAbs + (bexdeyAbs + dexbeyAbs) * aezAbs + (dexaeyAbs + aexdeyAbs) * bezAbs) * cLift + ((bexceyAbs + cexbeyAbs) * aezAbs + (cexaeyAbs + aexceyAbs) * bezAbs + (aexbeyAbs + bexaeyAbs) * cezAbs) * dLift;
    const double errorBound = inSphereErrorBound * permanent;

    // The determinant is positive inside of the sphere for points oriented the other way round.
    if (det > errorBound || -det > errorBound)
    {
        return -det;
    }
    return -inSphereExact(pa, pb, pc, pd, pe);
}
