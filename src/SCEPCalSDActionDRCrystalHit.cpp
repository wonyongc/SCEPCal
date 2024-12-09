//===============================
// Author: Wonyong Chung
//         Princeton University
//===============================
#include "DRCrystalHit.h"
#include "SCEPCalSegmentation.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"

namespace SCEPCal {
  G4double convertEvtoNm(G4double energy)
  {
    return 1239.84187/energy*1000.;
  }

  class SegmentedCrystalCalorimeterSD_DRHit {
    public:
      typedef DRCrystalHit Hit;
  };

  class SegmentedCrystalCalorimeterSD_DRHitSimple {
    public:
      typedef DRCrystalHitSimple Hit;
  };
}

namespace dd4hep {
  namespace sim {
    using namespace SCEPCal;
    
    template <> void Geant4SensitiveAction<SegmentedCrystalCalorimeterSD_DRHit>::defineCollections()    {
      m_collectionID = declareReadoutFilteredCollection<SegmentedCrystalCalorimeterSD_DRHit::Hit>();
    }

    template <> void Geant4SensitiveAction<SegmentedCrystalCalorimeterSD_DRHitSimple>::defineCollections()    {
      m_collectionID = declareReadoutFilteredCollection<SegmentedCrystalCalorimeterSD_DRHitSimple::Hit>();
    }

    template <> bool 
    Geant4SensitiveAction<SegmentedCrystalCalorimeterSD_DRHit>::process(const G4Step* step,G4TouchableHistory* /*hist*/ ) {
      G4double edep = step->GetTotalEnergyDeposit();
      G4StepPoint *thePrePoint = step->GetPreStepPoint();
      G4StepPoint *thePostPoint = step->GetPostStepPoint();
      G4TouchableHandle  thePreStepTouchable = thePrePoint->GetTouchableHandle();
      G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
      G4double pretime = thePrePoint->GetGlobalTime();
      G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();
      G4double posttime = thePostPoint->GetGlobalTime();
      G4String thePrePVName = "";
      if (thePrePV) thePrePVName = thePrePV->GetName();
      G4String thePostPVName = "";
      if (thePostPV) thePostPVName = thePostPV->GetName();

      Geant4StepHandler    h(step);
      Geant4HitData::MonteCarloContrib contrib = Geant4HitData::extractContribution(step);
      Geant4HitCollection* coll    = collection(m_collectionID);

      dd4hep::Segmentation* _geoSeg = &m_segmentation;
      auto segmentation=dynamic_cast<dd4hep::DDSegmentation::SCEPCalSegmentation *>(_geoSeg->segmentation());
      auto copyNum64 = segmentation->convertFirst32to64(thePreStepTouchable->GetCopyNumber(0));
      int cellID = (int)copyNum64;

      SegmentedCrystalCalorimeterSD_DRHit::Hit* hit = coll->findByKey<SegmentedCrystalCalorimeterSD_DRHit::Hit>(cellID);
      if(!hit) {    
        DDSegmentation::Vector3D pos = segmentation->myPosition(copyNum64);    
        
        
        Position global(pos.x(),pos.y(),pos.z());
    
        hit = new SegmentedCrystalCalorimeterSD_DRHit::Hit(global);
        hit->cellID = cellID;
        hit->system = segmentation->System(copyNum64);
        hit->eta = segmentation->Eta(copyNum64);
        hit->phi = segmentation->Phi(copyNum64);
        hit->depth = segmentation->Depth(copyNum64);
        coll->add(cellID, hit);

      }
      G4Track * track =  step->GetTrack();


      if(track->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()) {

        float wavelength=convertEvtoNm(track->GetTotalEnergy()/eV);

        int ibin=-1;
        float binsize=(hit->wavelen_max-hit->wavelen_min)/hit->nbins;
        ibin = (wavelength-hit->wavelen_min)/binsize;

        float avgarrival=(pretime+posttime)/2.;

        int jbin=-1;
        float tbinsize=(hit->time_max-hit->time_min)/hit->nbins;
        jbin = (avgarrival-hit->time_min)/tbinsize;

        int phstep = track->GetCurrentStepNumber();

        // if outside silicon, count 1st and allow to continue
        // if inside silicon, always kill, but count if not 1st step (came in from other material)
        // doesn't count photons generated inside silicon

        if (track->GetCreatorProcess()->G4VProcess::GetProcessName()=="CerenkovPhys") {
          
          std::string amedia = ((track->GetMaterial())->GetName());
          
          if(amedia.find("Silicon")!=std::string::npos) {
            if(phstep>1) {
              hit->ncerenkov+=1;
              if(ibin>-1 && ibin<hit->nbins) ((hit->nwavelen_cer).at(ibin)) +=1;
              if(jbin>-1 && jbin<hit->nbins) ((hit->ntime_cer).at(jbin)) +=1;
            }
            track->SetTrackStatus(fStopAndKill);
          }
          else {
            if(phstep==1) {
              hit->ncerenkov+=1;
              if(ibin>-1 && ibin<hit->nbins) ((hit->nwavelen_cer).at(ibin)) +=1;
              if(jbin>-1 && jbin<hit->nbins) ((hit->ntime_cer).at(jbin)) +=1;
            }
          }
        } 

        else if (track->GetCreatorProcess()->G4VProcess::GetProcessName()=="ScintillationPhys") {
          
          std::string amedia = ((track->GetMaterial())->GetName());
          
          if(amedia.find("Silicon")!=std::string::npos) {
            if(phstep>1) {
              hit->nscintillation+=1;
              if((ibin>-1)&&(ibin<hit->nbins)) ((hit->nwavelen_scint).at(ibin))+=1;
              if(jbin>-1&&jbin<hit->nbins) ((hit->ntime_scint).at(jbin))+=1;
            }
            track->SetTrackStatus(fStopAndKill);
          }
          else {
            if((track->GetCurrentStepNumber()==1)) {
              hit->nscintillation+=1; 
              if((ibin>-1)&&(ibin<hit->nbins)) ((hit->nwavelen_scint).at(ibin))+=1;
              if(jbin>-1&&jbin<hit->nbins) ((hit->ntime_scint).at(jbin))+=1;
            }
          }
        }


      }

      hit->truth.emplace_back(contrib);
      hit->energyDeposit+=edep;

      mark(h.track);

      return true;
    }

    template <> bool 
    Geant4SensitiveAction<SegmentedCrystalCalorimeterSD_DRHitSimple>::process(const G4Step* step,G4TouchableHistory* /*hist*/ ) {
      G4double edep = step->GetTotalEnergyDeposit();
      G4StepPoint *thePrePoint = step->GetPreStepPoint();
      G4StepPoint *thePostPoint = step->GetPostStepPoint();
      G4TouchableHandle  thePreStepTouchable = thePrePoint->GetTouchableHandle();
      G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
      G4double pretime = thePrePoint->GetGlobalTime();
      G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();
      G4double posttime = thePostPoint->GetGlobalTime();
      G4String thePrePVName = "";
      if (thePrePV) thePrePVName = thePrePV->GetName();
      G4String thePostPVName = "";
      if (thePostPV) thePostPVName = thePostPV->GetName();

      Geant4StepHandler    h(step);
      Geant4HitData::MonteCarloContrib contrib = Geant4HitData::extractContribution(step);
      Geant4HitCollection* coll    = collection(m_collectionID);

      dd4hep::Segmentation* _geoSeg = &m_segmentation;
      auto segmentation=dynamic_cast<dd4hep::DDSegmentation::SCEPCalSegmentation *>(_geoSeg->segmentation());
      auto copyNum64 = segmentation->convertFirst32to64(thePreStepTouchable->GetCopyNumber(0));
      int cellID = (int)copyNum64;

      SegmentedCrystalCalorimeterSD_DRHitSimple::Hit* hit = coll->findByKey<SegmentedCrystalCalorimeterSD_DRHitSimple::Hit>(cellID);
      if(!hit) {    
        DDSegmentation::Vector3D pos = segmentation->myPosition(copyNum64);    
        
        Position global(pos.x(),pos.y(),pos.z());
    
        hit = new SegmentedCrystalCalorimeterSD_DRHitSimple::Hit(global);
        hit->cellID = cellID;
        coll->add(cellID, hit);

      }
      G4Track * track =  step->GetTrack();

      if(track->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()) {

        float avgarrival=(pretime+posttime)/2.;

        // count 1st and kill
        // apply scale factor and poisson smearing

        int phstep = track->GetCurrentStepNumber();

        if (track->GetCreatorProcess()->G4VProcess::GetProcessName()=="CerenkovPhys") {
          if(phstep==1) {
            float tAvgC_new = (((hit->tAvgC)*(hit->nCerenkovProd)) +avgarrival)/(hit->nCerenkovProd+1);
            hit->nCerenkovProd+=1;
            hit->tAvgC = tAvgC_new;
          }
          track->SetTrackStatus(fStopAndKill);
        } 

        else if (track->GetCreatorProcess()->G4VProcess::GetProcessName()=="ScintillationPhys") {
          if(phstep==1) {
            float tAvgS_new = (((hit->tAvgS)*(hit->nScintillationProd)) +avgarrival)/(hit->nScintillationProd+1);
            hit->nScintillationProd+=1;
            hit->tAvgS = tAvgS_new;
          }
          track->SetTrackStatus(fStopAndKill);
        }

      }

      hit->truth.emplace_back(contrib);
      hit->energyDeposit+=edep;

      mark(h.track);

      return true;
    }

  }
}


namespace dd4hep { namespace sim {
    typedef Geant4SensitiveAction<SegmentedCrystalCalorimeterSD_DRHit> SCEPCalSDAction_DRHit;
    typedef Geant4SensitiveAction<SegmentedCrystalCalorimeterSD_DRHitSimple> SCEPCalSDAction_DRHitSimple;
  }}
DECLARE_GEANT4SENSITIVE(SCEPCalSDAction_DRHit)
DECLARE_GEANT4SENSITIVE(SCEPCalSDAction_DRHitSimple)