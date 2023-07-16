//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#ifndef SCEPCalSegmentationHandle_h
#define SCEPCalSegmentationHandle_h 1
#include "SCEPCalSegmentation.h"
#include "DD4hep/Segmentations.h"
#include "DD4hep/detail/SegmentationsInterna.h"

namespace dd4hep {
    
class Segmentation;
template <typename T>
class SegmentationWrapper;

typedef Handle<SegmentationWrapper<DDSegmentation::SCEPCalSegmentation>> SCEPCalSegmentationHandle;

class SCEPCalSegmentation : public SCEPCalSegmentationHandle {
public:
    typedef SCEPCalSegmentationHandle::Object Object;

public:
    SCEPCalSegmentation() = default;
    SCEPCalSegmentation(const SCEPCalSegmentation& e) = default;
    SCEPCalSegmentation(const Segmentation& e) : Handle<Object>(e) {}
    SCEPCalSegmentation(const Handle<Object>& e) : Handle<Object>(e) {}
    template <typename Q>
    SCEPCalSegmentation(const Handle<Q>& e) : Handle<Object>(e) {}
    SCEPCalSegmentation& operator=(const SCEPCalSegmentation& seg) = default;
    bool operator==(const SCEPCalSegmentation& seg) const { return m_element == seg.m_element; }

    inline Position position(const CellID& id) const {
        return Position(access()->implementation->position(id));
    }

    inline dd4hep::CellID cellID(const Position& local, const Position& global, const VolumeID& volID) const {
        return access()->implementation->cellID(local, global, volID);
    }

    inline VolumeID setVolumeID(int System, int Eta, int Phi, int Depth) const {
        return access()->implementation->setVolumeID(System, Eta,Phi,Depth);
    }

    inline CellID setCellID(int System, int Eta, int Phi, int Depth) const {
        return access()->implementation->setCellID(System, Eta, Phi, Depth);
    }

    inline int System(const CellID& aCellID) const { return access()->implementation->System(aCellID); }
    inline int Eta(const CellID& aCellID) const { return access()->implementation->Eta(aCellID); }
    inline int Phi(const CellID& aCellID) const { return access()->implementation->Phi(aCellID); }
    inline int Depth(const CellID& aCellID) const { return access()->implementation->Depth(aCellID); }

    inline int getFirst32bits(const CellID& aCellID) const { return access()->implementation->getFirst32bits(aCellID); }
    inline int getLast32bits(const CellID& aCellID) const { return access()->implementation->getLast32bits(aCellID); }

    inline CellID convertFirst32to64(const int aId32) const { return access()->implementation->convertFirst32to64(aId32); }
    inline CellID convertLast32to64(const int aId32) const { return access()->implementation->convertLast32to64(aId32); }

    inline int System(const int& aId32) const { return access()->implementation->System(aId32); }
    inline int Eta(const int& aId32) const { return access()->implementation->Eta(aId32); }
    inline int Phi(const int& aId32) const { return access()->implementation->Phi(aId32); }
    inline int Depth(const int& aId32) const { return access()->implementation->Depth(aId32); }

    inline void setGeomParams(double Fdz, double Rdz,
                              double nomfw, double nomth, 
                              double EBz, double Rin,
                              double sipmth,
                              int phiSegments,
                              int NProjectiveFill) const { 
        return access()->implementation->setGeomParams(Fdz,Rdz,nomfw,nomth,EBz,Rin,sipmth,phiSegments,NProjectiveFill);
    }

    inline double getFdz() const{ return access()->implementation->getFdz(); }
    inline double getRdz() const{ return access()->implementation->getRdz(); }
    inline double getnomfw() const{ return access()->implementation->getnomfw(); }
    inline double getnomth() const{ return access()->implementation->getnomth(); }
    inline double getEBz() const{ return access()->implementation->getEBz(); }
    inline double getRin() const{ return access()->implementation->getRin(); }
    inline double getSipmth() const{ return access()->implementation->getSipmth(); }
    inline int    getphiSegments() const{ return access()->implementation->getphiSegments(); }
    inline int    getNProjectiveFill() const{ return access()->implementation->getNProjectiveFill(); }

};

}
#endif