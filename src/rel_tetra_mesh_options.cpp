#include "rel_tetra_mesh_options.h"

RTetraMeshOptions::RTetraMeshOptions()
    : verbose(false)
    , reconstruct(false)
    , qualityMesh(true)
    , radiusEdgeRatio(2.0)
    , minDihedralAngle(0.0)
    , volumeConstraint(0.0)
    , useSizeFunction(false)
    , tolerance(1.0e-8)
    , maxSteinerNodes(0)
    , elementOrder(1)
{
}

bool RTetraMeshOptions::getVerbose() const
{
    return this->verbose;
}

void RTetraMeshOptions::setVerbose(bool verbose)
{
    this->verbose = verbose;
}

bool RTetraMeshOptions::getReconstruct() const
{
    return this->reconstruct;
}

void RTetraMeshOptions::setReconstruct(bool reconstruct)
{
    this->reconstruct = reconstruct;
}

bool RTetraMeshOptions::getQualityMesh() const
{
    return this->qualityMesh;
}

void RTetraMeshOptions::setQualityMesh(bool qualityMesh)
{
    this->qualityMesh = qualityMesh;
}

double RTetraMeshOptions::getRadiusEdgeRatio() const
{
    return this->radiusEdgeRatio;
}

void RTetraMeshOptions::setRadiusEdgeRatio(double radiusEdgeRatio)
{
    this->radiusEdgeRatio = radiusEdgeRatio;
}

double RTetraMeshOptions::getMinDihedralAngle() const
{
    return this->minDihedralAngle;
}

void RTetraMeshOptions::setMinDihedralAngle(double minDihedralAngle)
{
    this->minDihedralAngle = minDihedralAngle;
}

double RTetraMeshOptions::getVolumeConstraint() const
{
    return this->volumeConstraint;
}

void RTetraMeshOptions::setVolumeConstraint(double volumeConstraint)
{
    this->volumeConstraint = volumeConstraint;
}

bool RTetraMeshOptions::getUseSizeFunction() const
{
    return this->useSizeFunction;
}

void RTetraMeshOptions::setUseSizeFunction(bool useSizeFunction)
{
    this->useSizeFunction = useSizeFunction;
}

double RTetraMeshOptions::getTolerance() const
{
    return this->tolerance;
}

void RTetraMeshOptions::setTolerance(double tolerance)
{
    this->tolerance = tolerance;
}

unsigned int RTetraMeshOptions::getMaxSteinerNodes() const
{
    return this->maxSteinerNodes;
}

void RTetraMeshOptions::setMaxSteinerNodes(unsigned int maxSteinerNodes)
{
    this->maxSteinerNodes = maxSteinerNodes;
}

int RTetraMeshOptions::getElementOrder() const
{
    return this->elementOrder;
}

void RTetraMeshOptions::setElementOrder(int elementOrder)
{
    this->elementOrder = elementOrder;
}
