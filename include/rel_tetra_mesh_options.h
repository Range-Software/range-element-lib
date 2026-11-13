#ifndef REL_TETRA_MESH_OPTIONS_H
#define REL_TETRA_MESH_OPTIONS_H

#include <string>

//! Options of the tetrahedral mesh generator.
class RTetraMeshOptions
{

    protected:

        //! Print details of mesh generation.
        bool verbose;
        //! Input tetrahedra define regions and their interfaces (existing mesh is reconstructed).
        bool reconstruct;
        //! Refine mesh to meet quality and size criteria.
        bool qualityMesh;
        //! Maximum allowed ratio of circumscribed sphere radius to the shortest edge of a tetrahedron.
        double radiusEdgeRatio;
        //! Minimum allowed dihedral angle of a tetrahedron in degrees, zero if it is not limited.
        double minDihedralAngle;
        //! Maximum allowed volume of a tetrahedron, zero if it is not limited.
        double volumeConstraint;
        //! Use mesh sizes prescribed at input nodes.
        bool useSizeFunction;
        //! Maximum sine of an angle between normals of two neighboring triangles
        //! for them to be considered as lying in one plane.
        double tolerance;
        //! Maximum number of nodes the generator may add, zero if it is not limited.
        unsigned int maxSteinerNodes;
        std::string qualityReportPath;

    public:

        //! Constructor.
        RTetraMeshOptions();

        //! Return verbose.
        bool getVerbose() const;

        //! Set verbose.
        void setVerbose(bool verbose);

        //! Return whether an existing mesh is reconstructed.
        bool getReconstruct() const;

        //! Set whether an existing mesh is reconstructed.
        void setReconstruct(bool reconstruct);

        //! Return whether to generate quality mesh.
        bool getQualityMesh() const;

        //! Set whether to generate quality mesh.
        void setQualityMesh(bool qualityMesh);

        //! Return radius/edge ratio.
        double getRadiusEdgeRatio() const;

        //! Set radius/edge ratio.
        void setRadiusEdgeRatio(double radiusEdgeRatio);

        //! Return minimum dihedral angle in degrees.
        double getMinDihedralAngle() const;

        //! Set minimum dihedral angle in degrees.
        void setMinDihedralAngle(double minDihedralAngle);

        //! Return volume constraint.
        double getVolumeConstraint() const;

        //! Set volume constraint.
        void setVolumeConstraint(double volumeConstraint);

        //! Return whether to use mesh sizes prescribed at input nodes.
        bool getUseSizeFunction() const;

        //! Set whether to use mesh sizes prescribed at input nodes.
        void setUseSizeFunction(bool useSizeFunction);

        //! Return tolerance.
        double getTolerance() const;

        //! Set tolerance.
        void setTolerance(double tolerance);

        //! Return maximum number of added nodes.
        unsigned int getMaxSteinerNodes() const;

        //! Set maximum number of added nodes.
        void setMaxSteinerNodes(unsigned int maxSteinerNodes);

        //! Optional JSON diagnostic export. Empty means no file is written.
        const std::string &getQualityReportPath() const { return qualityReportPath; }
        void setQualityReportPath(const std::string &path) { qualityReportPath=path; }

};

#endif // REL_TETRA_MESH_OPTIONS_H
