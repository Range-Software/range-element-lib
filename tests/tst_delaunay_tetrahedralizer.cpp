#include <QTest>

#include <cmath>
#include <cstdint>
#include <vector>

#include <rbl_error.h>

#include "rel_delaunay_tetrahedralizer.h"
#include "rel_geometric_predicates.h"
#include "rel_tetra_mesh.h"

class TestDelaunayTetrahedralizer : public QObject
{
    Q_OBJECT

    private:

        //! Repeatable pseudo random number between 0 and 1.
        static double findRandom(uint64_t &seed)
        {
            seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
            return double(seed >> 11) / double(1ULL << 53);
        }

        //! Check orientation of tetrahedra and that neighbors refer to each other.
        static void verifyTopology(const RTetraMesh &mesh)
        {
            for (int i=0;i<mesh.getNTetraSlots();i++)
            {
                if (!mesh.isAlive(i))
                {
                    continue;
                }
                const RTetraMesh::Tetra &rTetra = mesh.getTetra(i);
                if (!mesh.isGhost(i))
                {
                    QVERIFY(RGeometricPredicates::orient3d(mesh.getPoint(rTetra.v[0]),
                                                           mesh.getPoint(rTetra.v[1]),
                                                           mesh.getPoint(rTetra.v[2]),
                                                           mesh.getPoint(rTetra.v[3])) > 0.0);
                }
                for (int j=0;j<4;j++)
                {
                    const int neighborID = rTetra.n[j];
                    QVERIFY(neighborID != RTetraMesh::noTetra);
                    QVERIFY(mesh.isAlive(neighborID));
                    const RTetraMesh::Tetra &rNeighbor = mesh.getTetra(neighborID);
                    int nBack = 0;
                    for (int k=0;k<4;k++)
                    {
                        if (rNeighbor.n[k] == i)
                        {
                            nBack++;
                            // The shared face is made of the same vertices.
                            for (int l=0;l<4;l++)
                            {
                                if (l != k)
                                {
                                    bool found = false;
                                    for (int m=0;m<4;m++)
                                    {
                                        found = found || (m != j && rTetra.v[m] == rNeighbor.v[l]);
                                    }
                                    QVERIFY(found);
                                }
                            }
                        }
                    }
                    QCOMPARE(nBack,1);
                }
            }
        }

        //! Check that no point lies inside of a circumscribed sphere of any tetrahedron.
        static void verifyDelaunay(const RTetraMesh &mesh)
        {
            for (int i=0;i<mesh.getNTetraSlots();i++)
            {
                if (!mesh.isAlive(i) || mesh.isGhost(i))
                {
                    continue;
                }
                const RTetraMesh::Tetra &rTetra = mesh.getTetra(i);
                for (int j=0;j<mesh.getNPoints();j++)
                {
                    QVERIFY(RGeometricPredicates::inSphere(mesh.getPoint(rTetra.v[0]),
                                                           mesh.getPoint(rTetra.v[1]),
                                                           mesh.getPoint(rTetra.v[2]),
                                                           mesh.getPoint(rTetra.v[3]),
                                                           mesh.getPoint(j)) <= 0.0);
                }
            }
        }

        static double findVolume(const RTetraMesh &mesh)
        {
            double volume = 0.0;
            for (int i=0;i<mesh.getNTetraSlots();i++)
            {
                if (!mesh.isAlive(i) || mesh.isGhost(i))
                {
                    continue;
                }
                const RTetraMesh::Tetra &rTetra = mesh.getTetra(i);
                volume += RGeometricPredicates::orient3d(mesh.getPoint(rTetra.v[0]),
                                                         mesh.getPoint(rTetra.v[1]),
                                                         mesh.getPoint(rTetra.v[2]),
                                                         mesh.getPoint(rTetra.v[3])) / 6.0;
            }
            return volume;
        }

private slots:

    void randomPoints();
    void gridPoints();
    void duplicatePoints();
    void coplanarPoints();
    void findEdgesAndFaces();
};

void TestDelaunayTetrahedralizer::randomPoints()
{
    RTetraMesh mesh;
    uint64_t seed = 12345;
    for (int i=0;i<150;i++)
    {
        const double x = findRandom(seed);
        const double y = findRandom(seed);
        const double z = findRandom(seed);
        mesh.addPoint(x,y,z);
    }

    std::vector<int> duplicates;
    RDelaunayTetrahedralizer::tetrahedralize(mesh,duplicates);

    for (int duplicateID : duplicates)
    {
        QCOMPARE(duplicateID,-1);
    }
    for (int i=0;i<mesh.getNPoints();i++)
    {
        QVERIFY(mesh.getVertexTetra(i) != RTetraMesh::noTetra);
    }
    verifyTopology(mesh);
    verifyDelaunay(mesh);
}

void TestDelaunayTetrahedralizer::gridPoints()
{
    // Points of a regular grid lie on common spheres and planes, which is the worst case for predicates.
    RTetraMesh mesh;
    for (int i=0;i<5;i++)
    {
        for (int j=0;j<5;j++)
        {
            for (int k=0;k<5;k++)
            {
                mesh.addPoint(0.1*i,0.1*j,0.1*k);
            }
        }
    }

    std::vector<int> duplicates;
    RDelaunayTetrahedralizer::tetrahedralize(mesh,duplicates);

    verifyTopology(mesh);
    verifyDelaunay(mesh);

    // Tetrahedra fill the whole grid.
    QVERIFY(std::fabs(findVolume(mesh) - 0.4*0.4*0.4) < 1.0e-12);
}

void TestDelaunayTetrahedralizer::duplicatePoints()
{
    RTetraMesh mesh;
    mesh.addPoint(0.0,0.0,0.0);
    mesh.addPoint(1.0,0.0,0.0);
    mesh.addPoint(0.0,1.0,0.0);
    mesh.addPoint(0.0,0.0,1.0);
    mesh.addPoint(0.3,0.3,0.3);
    mesh.addPoint(0.3,0.3,0.3);
    mesh.addPoint(1.0,0.0,0.0);

    std::vector<int> duplicates;
    RDelaunayTetrahedralizer::tetrahedralize(mesh,duplicates);

    int nDuplicates = 0;
    for (int i=0;i<mesh.getNPoints();i++)
    {
        const int duplicateID = duplicates[size_t(i)];
        if (duplicateID < 0)
        {
            continue;
        }
        nDuplicates++;
        // Duplicate refers to an inserted point with the same coordinates.
        QCOMPARE(duplicates[size_t(duplicateID)],-1);
        for (int j=0;j<3;j++)
        {
            QCOMPARE(mesh.getPoint(i)[j],mesh.getPoint(duplicateID)[j]);
        }
        QCOMPARE(mesh.getVertexTetra(i),int(RTetraMesh::noTetra));
    }
    QCOMPARE(nDuplicates,2);

    verifyTopology(mesh);
    verifyDelaunay(mesh);
}

void TestDelaunayTetrahedralizer::coplanarPoints()
{
    RTetraMesh mesh;
    mesh.addPoint(0.0,0.0,0.5);
    mesh.addPoint(1.0,0.0,0.5);
    mesh.addPoint(0.0,1.0,0.5);
    mesh.addPoint(1.0,1.0,0.5);
    mesh.addPoint(0.3,0.6,0.5);

    std::vector<int> duplicates;
    QVERIFY_THROWS_EXCEPTION(RError,RDelaunayTetrahedralizer::tetrahedralize(mesh,duplicates));
}

void TestDelaunayTetrahedralizer::findEdgesAndFaces()
{
    RTetraMesh mesh;
    mesh.addPoint(0.0,0.0,0.0);
    mesh.addPoint(1.0,0.0,0.0);
    mesh.addPoint(0.0,1.0,0.0);
    mesh.addPoint(0.0,0.0,1.0);
    mesh.addPoint(0.2,0.2,-1.0);

    std::vector<int> duplicates;
    RDelaunayTetrahedralizer::tetrahedralize(mesh,duplicates);

    // Two tetrahedra sharing face 0-1-2.
    QVERIFY(mesh.findEdge(0,1) != RTetraMesh::noTetra);
    QVERIFY(mesh.findEdge(3,1) != RTetraMesh::noTetra);
    QCOMPARE(mesh.findEdge(3,4),int(RTetraMesh::noTetra));

    RTetraMesh::Face face;
    QVERIFY(mesh.findFace(0,1,2,face));
    QVERIFY(mesh.findFace(2,3,0,face));
    QVERIFY(!mesh.findFace(3,4,0,face));

    std::vector<int> star;
    mesh.findStar(3,star);
    for (int tetraID : star)
    {
        const RTetraMesh::Tetra &rTetra = mesh.getTetra(tetraID);
        QVERIFY(rTetra.v[0] == 3 || rTetra.v[1] == 3 || rTetra.v[2] == 3 || rTetra.v[3] == 3);
    }
    // One tetrahedron and three ghost ones.
    QCOMPARE(int(star.size()),4);
}

QTEST_APPLESS_MAIN(TestDelaunayTetrahedralizer)
#include "tst_delaunay_tetrahedralizer.moc"
