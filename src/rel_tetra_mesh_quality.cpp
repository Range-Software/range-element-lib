#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <map>
#include <future>
#include <thread>
#include <system_error>
#include <QJsonArray>
#include <QJsonObject>
#include <QJsonDocument>
#include <QSaveFile>
#include <rbl_error.h>
#include "rel_tetra_mesh_quality.h"
#include "rel_tetra_mesh.h"
#include "rel_geometric_predicates.h"

namespace
{
using V=std::array<double,3>;
V sub(const double *a,const double *b) { return {a[0]-b[0],a[1]-b[1],a[2]-b[2]}; }
V cross(const V &a,const V &b) { return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}; }
double dot(const V &a,const V &b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
}

RTetraMeshQuality::Element RTetraMeshQuality::measure(const double *a,const double *b,const double *c,const double *d,const bool *prescribedFaces)
{
    Element result;
    const double *p[4]={a,b,c,d};
    double longest=0.0, shortest=std::numeric_limits<double>::max();
    for (int i=0;i<4;i++) for (int j=i+1;j<4;j++)
    {
        const double length=std::sqrt(dot(sub(p[i],p[j]),sub(p[i],p[j])));
        longest=std::max(longest,length); shortest=std::min(shortest,length);
    }
    result.radiusEdge=std::numeric_limits<double>::infinity();
    if (!(shortest>0.0) || !std::isfinite(longest)) return result;
    V q[4]={V{0,0,0},sub(b,a),sub(c,a),sub(d,a)};
    for (int i=1;i<4;i++) for (double &v : q[i]) v/=longest;
    const double det=RGeometricPredicates::orient3d(q[0].data(),q[1].data(),q[2].data(),q[3].data());
    result.volume=RGeometricPredicates::orient3d(a,b,c,d)/6.0;
    if (!(det>0.0)) return result;
    result.scaledVolume=det;
    double edges2=0;
    for (int i=0;i<4;i++) for (int j=i+1;j<4;j++) edges2+=dot(sub(q[i].data(),q[j].data()),sub(q[i].data(),q[j].data()));
    result.meanRatio=12.0*std::pow(det/2.0,2.0/3.0)/edges2;
    const V vw=cross(q[2],q[3]),wu=cross(q[3],q[1]),uv=cross(q[1],q[2]);
    V radius;
    for (int i=0;i<3;i++) radius[i]=(dot(q[1],q[1])*vw[i]+dot(q[2],q[2])*wu[i]+dot(q[3],q[3])*uv[i])/(2.0*det);
    for (int i=0;i<3;i++) result.circumcenter[size_t(i)]=a[i]+longest*radius[i];
    result.circumradius=std::sqrt(dot(radius,radius))*longest;
    result.radiusEdge=result.circumradius/shortest;
    V normals[4];
    for (int i=0;i<4;i++)
    {
        const int *f=RTetraMesh::faceVertices[i];
        normals[i]=cross(sub(q[f[1]].data(),q[f[0]].data()),sub(q[f[2]].data(),q[f[0]].data()));
    }
    result.minAngle=180.0; result.minFreeAngle=180.0; result.maxAngle=0.0;
    for (int i=0;i<4;i++) for (int j=i+1;j<4;j++)
    {
        // atan2 retains tiny nonzero angles which acos rounds to zero.
        const V product=cross(normals[i],normals[j]);
        const double angle=std::atan2(std::sqrt(dot(product,product)),-dot(normals[i],normals[j]))*180.0/std::acos(-1.0);
        result.minAngle=std::min(result.minAngle,angle); result.maxAngle=std::max(result.maxAngle,angle);
        if (!prescribedFaces || !prescribedFaces[i] || !prescribedFaces[j]) result.minFreeAngle=std::min(result.minFreeAngle,angle);
    }
    return result;
}

RTetraMeshQuality::Summary RTetraMeshQuality::summarize(const RTetraMeshData &data,const RTetraMeshOptions &options,unsigned int workers)
{
    std::set<std::array<int,3>> constraints;
    for (const auto &t : data.getTriangles())
    {
        std::array<int,3> key={t.nodeIDs[0],t.nodeIDs[1],t.nodeIDs[2]};
        std::sort(key.begin(),key.end()); constraints.insert(key);
    }
    struct Chunk { Summary summary;std::vector<double> ratios; };
    const auto &cells=data.getTetrahedra();
    const auto &fixedFaces=constraints;
    auto evaluate=[&](size_t begin,size_t end) {
        Chunk result;auto &summary=result.summary;auto &ratios=result.ratios;
        ratios.reserve(end-begin);
        for (size_t index=begin;index<end;index++)
        {
            const auto &t=cells[index];
            summary.elements++;
            bool valid=true;
            for (int id : t.nodeIDs) valid=valid && id>=0 && size_t(id)<data.getNodes().size();
            if (!valid)
            {
                summary.invalid++;summary.nearFlat++;
                summary.minAngle=0.0;summary.maxRadiusEdge=std::numeric_limits<double>::infinity();
                ratios.push_back(summary.maxRadiusEdge);continue;
            }
            double p[4][3];
            bool fixed[4];
            for (int i=0;i<4;i++)
            {
                const auto &v=data.getNodes()[size_t(t.nodeIDs[i])]; p[i][0]=v.x;p[i][1]=v.y;p[i][2]=v.z;
                const int *f=RTetraMesh::faceVertices[i];
                std::array<int,3> key={t.nodeIDs[f[0]],t.nodeIDs[f[1]],t.nodeIDs[f[2]]};
                std::sort(key.begin(),key.end());fixed[i]=fixedFaces.count(key)>0;
            }
            const Element e=measure(p[0],p[1],p[2],p[3]);
            summary.invalid+=!(e.volume>0.0) || !std::isfinite(e.volume);
            summary.nearFlat+=e.scaledVolume<1e-12;
            summary.minAngle=std::min(summary.minAngle,e.minAngle);
            summary.maxRadiusEdge=std::max(summary.maxRadiusEdge,e.radiusEdge);
            ratios.push_back(e.radiusEdge);
            summary.badRadius+=options.getRadiusEdgeRatio()>0.0 && e.radiusEdge>options.getRadiusEdgeRatio()*(1.0+1e-10);
            bool tooLarge=options.getVolumeConstraint()>0.0 && e.volume>options.getVolumeConstraint()*(1.0+1e-10);
            if (options.getUseSizeFunction())
            {
                double size=0;int count=0;
                for (int id : t.nodeIDs) if (data.getNodes()[size_t(id)].size>0) { size+=data.getNodes()[size_t(id)].size;count++; }
                tooLarge=tooLarge || (count>0 && e.circumradius>.75*size/count*(1.0+1e-10));
            }
            summary.tooLarge+=tooLarge;
            bool bad=false, forced=false, freeOne=false, prescribedOne=false;
            for (int i=0;i<4;i++) for (int j=i+1;j<4;j++)
            {
                const int *fi=RTetraMesh::faceVertices[i], *fj=RTetraMesh::faceVertices[j];
                const V n=cross(sub(p[fi[1]],p[fi[0]]),sub(p[fi[2]],p[fi[0]]));
                const V m=cross(sub(p[fj[1]],p[fj[0]]),sub(p[fj[2]],p[fj[0]]));
                const V product=cross(n,m);
                const double angle=std::atan2(std::sqrt(dot(product,product)),-dot(n,m))*180.0/std::acos(-1.0);
                if (angle<1.0) { if (fixed[i] && fixed[j]) prescribedOne=true;else freeOne=true; }
                if (angle+1e-8<options.getMinDihedralAngle())
                {
                    if (fixed[i] && fixed[j]) forced=true; else bad=true;
                }
            }
            summary.badAngle+=bad; summary.constrainedAngle+=forced;
            summary.freeBelowOne+=freeOne;summary.prescribedBelowOne+=prescribedOne;
        }
        return result;
    };
    const size_t count=cells.size();
    if (workers==0) workers=count>=32768?std::max(1u,std::min(4u,std::thread::hardware_concurrency())):1;
    workers=std::min(4u,std::max(1u,workers));
    workers=unsigned(std::min(size_t(workers),std::max(size_t(1),count)));
    std::vector<std::future<Chunk>> tasks;
    try
    {
        for (unsigned int i=1;i<workers;i++) tasks.push_back(std::async(std::launch::async,evaluate,count*i/workers,count*(i+1)/workers));
    }
    catch (const std::system_error &)
    {
        // Reporting remains available if the platform cannot start a worker.
        for (auto &task:tasks) task.wait();tasks.clear();workers=1;
    }
    auto result=evaluate(0,count/workers);
    auto &summary=result.summary;auto &ratios=result.ratios;
    for (auto &task:tasks)
    {
        auto chunk=task.get();const auto &s=chunk.summary;
        summary.elements+=s.elements;
        summary.invalid+=s.invalid;
        summary.nearFlat+=s.nearFlat;
        summary.badRadius+=s.badRadius;
        summary.badAngle+=s.badAngle;
        summary.constrainedAngle+=s.constrainedAngle;
        summary.tooLarge+=s.tooLarge;
        summary.freeBelowOne+=s.freeBelowOne;
        summary.prescribedBelowOne+=s.prescribedBelowOne;
        summary.minAngle=std::min(summary.minAngle,s.minAngle);
        summary.maxRadiusEdge=std::max(summary.maxRadiusEdge,s.maxRadiusEdge);
        ratios.insert(ratios.end(),chunk.ratios.begin(),chunk.ratios.end());
    }
    std::sort(ratios.begin(),ratios.end());
    if (!ratios.empty()) summary.radiusEdge95=ratios[size_t(.95*double(ratios.size()-1))];
    return summary;
}

void RTetraMeshQuality::writeReport(const RTetraMeshData &data,const RTetraMeshOptions &options,
                                   const std::vector<PatchAttempts> &attempts,const std::string &path,const Work *work)
{
    using Key=std::array<int,3>;
    std::set<Key> constraints;
    for (const auto &face : data.getTriangles())
    {
        Key key={face.nodeIDs[0],face.nodeIDs[1],face.nodeIDs[2]};std::sort(key.begin(),key.end());constraints.insert(key);
    }
    std::set<int> surfaceVertices;
    for (const auto &face : data.getTriangles()) for (int v : face.nodeIDs) surfaceVertices.insert(v);
    std::map<std::array<int,4>,PatchAttempts> observed;
    for (const auto &record : attempts) observed[record.vertices]=record;
    std::vector<std::pair<double,size_t>> ranked;
    for (size_t id=0;id<data.getTetrahedra().size();id++)
    {
        const auto &t=data.getTetrahedra()[id];double p[4][3];bool valid=true;
        for (int i=0;i<4;i++)
        {
            if (t.nodeIDs[i]<0 || size_t(t.nodeIDs[i])>=data.getNodes().size()) { valid=false;break; }
            const auto &v=data.getNodes()[size_t(t.nodeIDs[i])];p[i][0]=v.x;p[i][1]=v.y;p[i][2]=v.z;
        }
        const auto q=valid?measure(p[0],p[1],p[2],p[3]):Element{};
        const double score=valid?std::min({q.meanRatio,1.0/(1.0+q.radiusEdge*q.radiusEdge),q.minAngle/180.0}):-1;
        ranked.push_back({score,id});
    }
    std::sort(ranked.begin(),ranked.end());
    QJsonArray patches;
    auto number=[](double n)->QJsonValue { return std::isfinite(n)?QJsonValue(n):QJsonValue(); };
    for (size_t rank=0;rank<std::min(size_t(50),ranked.size());rank++)
    {
        const size_t id=ranked[rank].second;const auto &t=data.getTetrahedra()[id];
        QJsonArray vertices,coordinates,neighbors,faces,neighborCells,onSurface;
        double p[4][3]={};bool fixedFaces[4]={};bool valid=true;int constrained=0;
        for (int i=0;i<4;i++)
        {
            vertices.append(t.nodeIDs[i]);neighbors.append(t.neighborIDs[i]);onSurface.append(surfaceVertices.count(t.nodeIDs[i])>0);
            if (t.nodeIDs[i]>=0 && size_t(t.nodeIDs[i])<data.getNodes().size())
            {
                const auto &v=data.getNodes()[size_t(t.nodeIDs[i])];p[i][0]=v.x;p[i][1]=v.y;p[i][2]=v.z;
                coordinates.append(QJsonArray{v.x,v.y,v.z});
            }
            else valid=false;
            const int *f=RTetraMesh::faceVertices[i];Key key={t.nodeIDs[f[0]],t.nodeIDs[f[1]],t.nodeIDs[f[2]]};std::sort(key.begin(),key.end());
            const bool fixed=constraints.count(key)>0;faces.append(fixed);fixedFaces[i]=fixed;constrained+=fixed;
            const int n=t.neighborIDs[i];
            if (n>=0 && size_t(n)<data.getTetrahedra().size())
            {
                const auto &neighbor=data.getTetrahedra()[size_t(n)];QJsonArray points,positions;
                for (int v : neighbor.nodeIDs)
                {
                    points.append(v);
                    if (v>=0 && size_t(v)<data.getNodes().size())
                    { const auto &node=data.getNodes()[size_t(v)];positions.append(QJsonArray{node.x,node.y,node.z}); }
                }
                neighborCells.append(QJsonObject{{"element",n},{"vertices",points},{"coordinates",positions},{"region",neighbor.marker}});
            }
        }
        const auto q=valid?measure(p[0],p[1],p[2],p[3],fixedFaces):Element{};
        auto key=std::array<int,4>{t.nodeIDs[0],t.nodeIDs[1],t.nodeIDs[2],t.nodeIDs[3]};std::sort(key.begin(),key.end());
        auto record=observed.find(key);
        QJsonObject rejections;QStringList causes;
        if (record!=observed.end())
        {
            const auto &r=record->second;
            rejections={{"visits",double(r.visits)},{"accepted",double(r.accepted)},{"boundary",double(r.boundary)},{"region",double(r.region)},
                        {"topology",double(r.topology)},{"quality",double(r.quality)},{"spacing",double(r.spacing)},{"placement",double(r.placement)}};
            if (r.boundary) causes.append("prescribed boundary encountered");
            if (r.region) causes.append("region interface encountered");
            if (r.topology) causes.append("tested connectivity was not admissible");
            if (r.spacing) causes.append("insertion protection rejected candidates");
            if (r.quality || r.placement) causes.append("tested candidates did not improve the patch safely");
        }
        else causes.append("no repair attempt recorded for this final vertex set");
        if (causes.empty()) causes.append("bounded local search found no further improving candidate");
        // A protected face or edge is evidence of a constraint, not proof that
        // every defect in this cell is forced by the input geometry.
        QJsonArray explanation;for (const auto &cause : causes) explanation.append(cause);
        patches.append(QJsonObject{{"element",double(id)},{"vertices",vertices},{"coordinates",coordinates},{"verticesOnPrescribedSurfaces",onSurface},{"neighbors",neighbors},
                                   {"neighborCells",neighborCells},{"region",t.marker},{"constrainedFaces",faces},
                                   {"constrainedFaceCount",constrained},{"minimumDihedral",number(q.minAngle)},
                                   {"minimumFreeDihedral",number(q.minFreeAngle)},
                                   {"minimumOnPrescribedFaces",q.minAngle+1e-8<q.minFreeAngle},
                                   {"radiusEdge",number(q.radiusEdge)},{"scaledVolume",number(q.scaledVolume)},
                                   {"meanRatio",number(q.meanRatio)},{"attempts",rejections},{"observations",explanation}});
    }
    const auto summary=summarize(data,options);
    QJsonObject report{{"schemaVersion",1},{"scope","worst cells of the final mesh; rejection counts describe bounded attempted repairs, not impossibility proofs"},
                       {"elements",double(summary.elements)},{"targetsMet",summary.targetsMet()},
                       {"requestedRadiusEdge",options.getRadiusEdgeRatio()},{"requestedMinimumDihedral",options.getMinDihedralAngle()},
                       {"minimumDihedral",number(summary.minAngle)},{"radiusEdge95",number(summary.radiusEdge95)},
                       {"badRadius",double(summary.badRadius)},{"badAngle",double(summary.badAngle)},
                       {"freeAnglesBelowOne",double(summary.freeBelowOne)},{"prescribedAnglesBelowOne",double(summary.prescribedBelowOne)},
                       {"patches",patches}};
    if (work) report["optimizer"]=QJsonObject{{"queueVisits",double(work->visits)},{"queuesStoppedAtBudget",double(work->budgetStops)},{"acceptedOffCenters",double(work->offCenters)},{"orderedCandidates",double(work->orderedCandidates)},{"acceptedBoundaryPoints",double(work->boundaryPoints)}};
    if (work && work->profile)
    {
        const auto &p=*work->profile;
        report["profile"]=QJsonObject{
            {"total",double(p.total)},
            {"tetrahedralization",double(p.tetrahedralization)},
            {"boundary",double(p.boundary)},
            {"initialRepair",double(p.initialRepair)},
            {"refinement",double(p.refinement)},
            {"optimization",double(p.optimization)},
            {"diagnostics",double(p.diagnostics)},
            {"orderedRepair",double(p.orderedRepair)},
            {"localRepair",double(p.localRepair)},
            {"relocation",double(p.relocation)},
            {"compound",double(p.compound)},
            {"boundaryTrials",double(p.boundaryTrials)},
            {"qualityMeasurements",double(p.qualityMeasurements)},
            {"qualityCacheHits",double(p.qualityCacheHits)},
            {"starQueries",double(p.starQueries)},
            {"starCacheHits",double(p.starCacheHits)},
            {"trialCopies",double(p.trialCopies)},
            {"trialCopyBytes",double(p.trialCopyBytes)},
            {"journalBytes",double(p.journalBytes)},
            {"trialReplacements",double(p.trialReplacements)},
            {"queueSkips",double(p.queueSkips)},
            {"acceptedChanges",double(p.acceptedChanges)}
        };
    }
    QSaveFile file(QString::fromStdString(path));
    if (!file.open(QIODevice::WriteOnly)) throw RError(RError::Type::Application,R_ERROR_REF,"Cannot open mesh quality report '%s'.",path.c_str());
    const QByteArray bytes=QJsonDocument(report).toJson(QJsonDocument::Indented);
    if (file.write(bytes)!=bytes.size() || !file.commit()) throw RError(RError::Type::Application,R_ERROR_REF,"Cannot write mesh quality report '%s'.",path.c_str());
}
