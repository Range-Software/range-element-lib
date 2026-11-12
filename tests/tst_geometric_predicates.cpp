#include <QTest>

#include <cmath>

#include "rel_geometric_predicates.h"

class TestGeometricPredicates : public QObject
{
    Q_OBJECT

private slots:

    void orient3dSign();
    void orient3dCoplanar();
    void orient3dNearlyCoplanar();
    void orient3dPermutations();
    void inSphereSign();
    void inSphereCospherical();
    void inSphereNearlyCospherical();
    void inSphereRoundedCoordinates();
};

void TestGeometricPredicates::orient3dSign()
{
    const double a[3] = {0.0,0.0,0.0};
    const double b[3] = {1.0,0.0,0.0};
    const double c[3] = {0.0,1.0,0.0};
    const double above[3] = {0.0,0.0,1.0};
    const double below[3] = {0.0,0.0,-1.0};

    // Normal (b-a)x(c-a) points to +Z.
    QVERIFY(RGeometricPredicates::orient3d(a,b,c,above) > 0.0);
    QVERIFY(RGeometricPredicates::orient3d(a,b,c,below) < 0.0);
    QVERIFY(RGeometricPredicates::orient3d(a,c,b,above) < 0.0);
}

void TestGeometricPredicates::orient3dCoplanar()
{
    const double a[3] = {0.1,0.2,0.3};
    const double b[3] = {1.1,0.2,0.3};
    const double c[3] = {0.1,1.2,0.3};
    const double d[3] = {0.7,0.9,0.3};
    QCOMPARE(RGeometricPredicates::orient3d(a,b,c,d),0.0);

    // Points on a plane which is not aligned with axes.
    const double e[3] = {1.0,2.0,3.0};
    const double f[3] = {2.0,4.0,6.0};
    const double g[3] = {3.0,1.0,2.0};
    const double h[3] = {6.0,7.0,11.0};
    QCOMPARE(RGeometricPredicates::orient3d(e,f,g,h),0.0);
}

void TestGeometricPredicates::orient3dNearlyCoplanar()
{
    // Offsets which are lost when coordinates are subtracted in floating point arithmetic.
    const double a[3] = {1.0e10,1.0e10,0.0};
    const double b[3] = {1.0e10 + 1.0,1.0e10,0.0};
    const double c[3] = {1.0e10,1.0e10 + 1.0,0.0};
    const double above[3] = {1.0e10 + 0.25,1.0e10 + 0.25,1.0e-300};
    const double below[3] = {1.0e10 + 0.25,1.0e10 + 0.25,-1.0e-300};
    const double on[3] = {1.0e10 + 0.25,1.0e10 + 0.25,0.0};

    QVERIFY(RGeometricPredicates::orient3d(a,b,c,above) > 0.0);
    QVERIFY(RGeometricPredicates::orient3d(a,b,c,below) < 0.0);
    QCOMPARE(RGeometricPredicates::orient3d(a,b,c,on),0.0);

    // A point moved by the smallest possible step from a tilted plane.
    const double e[3] = {0.0,0.0,0.0};
    const double f[3] = {1.0,0.0,1.0};
    const double g[3] = {0.0,1.0,1.0};
    const double up[3] = {0.5,0.5,std::nextafter(1.0,2.0)};
    const double down[3] = {0.5,0.5,std::nextafter(1.0,0.0)};
    QVERIFY(RGeometricPredicates::orient3d(e,f,g,up) > 0.0);
    QVERIFY(RGeometricPredicates::orient3d(e,f,g,down) < 0.0);
}

void TestGeometricPredicates::orient3dPermutations()
{
    // Sign must change with each swap of two points, however close to a plane the points are.
    const double a[3] = {0.1,0.7,0.3};
    const double b[3] = {0.9,0.2,0.3 + 1.0e-17};
    const double c[3] = {0.4,0.5,0.3};
    const double d[3] = {0.3,0.1,std::nextafter(0.3,1.0)};

    const double abcd = RGeometricPredicates::orient3d(a,b,c,d);
    QVERIFY(abcd != 0.0);
    QVERIFY((RGeometricPredicates::orient3d(b,a,c,d) > 0.0) != (abcd > 0.0));
    QVERIFY((RGeometricPredicates::orient3d(a,b,d,c) > 0.0) != (abcd > 0.0));
    QVERIFY((RGeometricPredicates::orient3d(b,c,a,d) > 0.0) == (abcd > 0.0));
    QVERIFY((RGeometricPredicates::orient3d(d,c,b,a) > 0.0) == (abcd > 0.0));
}

void TestGeometricPredicates::inSphereSign()
{
    const double a[3] = {0.0,0.0,0.0};
    const double b[3] = {1.0,0.0,0.0};
    const double c[3] = {0.0,1.0,0.0};
    const double d[3] = {0.0,0.0,1.0};
    QVERIFY(RGeometricPredicates::orient3d(a,b,c,d) > 0.0);

    const double inside[3] = {0.25,0.25,0.25};
    const double outside[3] = {2.0,2.0,2.0};
    QVERIFY(RGeometricPredicates::inSphere(a,b,c,d,inside) > 0.0);
    QVERIFY(RGeometricPredicates::inSphere(a,b,c,d,outside) < 0.0);
}

void TestGeometricPredicates::inSphereCospherical()
{
    // Corners of a cube lie on one sphere.
    const double a[3] = {0.0,0.0,0.0};
    const double b[3] = {1.0,0.0,0.0};
    const double c[3] = {0.0,1.0,0.0};
    const double d[3] = {0.0,0.0,1.0};
    const double e[3] = {1.0,1.0,1.0};
    const double f[3] = {1.0,1.0,0.0};
    QCOMPARE(RGeometricPredicates::inSphere(a,b,c,d,e),0.0);
    QCOMPARE(RGeometricPredicates::inSphere(a,b,c,d,f),0.0);

    const double shift[3] = {0.1,0.2,0.3};
    double as[3], bs[3], cs[3], ds[3], es[3];
    for (int i=0;i<3;i++)
    {
        as[i] = a[i] * 0.5 + shift[i];
        bs[i] = b[i] * 0.5 + shift[i];
        cs[i] = c[i] * 0.5 + shift[i];
        ds[i] = d[i] * 0.5 + shift[i];
        es[i] = e[i] * 0.5 + shift[i];
    }
    // Shifted coordinates are rounded, but all corners are still built from the same six numbers.
    QCOMPARE(RGeometricPredicates::inSphere(as,bs,cs,ds,es),0.0);
}

void TestGeometricPredicates::inSphereNearlyCospherical()
{
    const double a[3] = {0.0,0.0,0.0};
    const double b[3] = {1.0,0.0,0.0};
    const double c[3] = {0.0,1.0,0.0};
    const double d[3] = {0.0,0.0,1.0};

    // Corner of a cube moved by the smallest possible step to the center of the sphere and from it.
    const double inside[3] = {1.0,1.0,std::nextafter(1.0,0.0)};
    const double outside[3] = {1.0,1.0,std::nextafter(1.0,2.0)};
    QVERIFY(RGeometricPredicates::inSphere(a,b,c,d,inside) > 0.0);
    QVERIFY(RGeometricPredicates::inSphere(a,b,c,d,outside) < 0.0);
}

void TestGeometricPredicates::inSphereRoundedCoordinates()
{
    // Corners of a rectangle and a point out of its plane lie on one sphere, also if their coordinates
    // are rounded, as long as the corners are built from the same numbers. Products of such coordinates
    // are not exact in floating point arithmetic, so this fails if any step of the exact evaluation is rounded.
    const double z2 = 0.1 * 2;
    const double z3 = 0.1 * 3;
    const double a[3] = {0.0,0.0,z2};
    const double b[3] = {0.1 * 4,0.0,0.0};
    const double c[3] = {0.0,0.1,z2};
    const double d[3] = {0.0,0.1,z3};
    const double e[3] = {0.0,0.0,z3};
    QVERIFY(RGeometricPredicates::orient3d(a,b,c,d) != 0.0);
    if (RGeometricPredicates::orient3d(a,b,c,d) > 0.0)
    {
        QCOMPARE(RGeometricPredicates::inSphere(a,b,c,d,e),0.0);
    }
    else
    {
        QCOMPARE(RGeometricPredicates::inSphere(a,c,b,d,e),0.0);
    }
}

QTEST_APPLESS_MAIN(TestGeometricPredicates)
#include "tst_geometric_predicates.moc"
