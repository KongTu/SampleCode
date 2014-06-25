// -*- C++ -*-
//
// Package:    V0Validator
// Class:      V0Validator
// 
/**\class V0Validator V0Validator.cc Validation/RecoVertex/interface/V0Validator.h

 Description: Creates validation histograms for RecoVertex/V0Producer

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Wed Feb 18 17:21:04 MST 2009
// $Id: V0Validator.h,v 1.5 2010/02/24 14:22:50 werdmann Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/Records/interface/VertexAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/VertexAssociation/interface/VertexAssociatorBase.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateClosestToBeamLineBuilder.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
//#include "RecoVertex/V0Producer/interface/V0Producer.h"

#include "SimTracker/TrackHistory/interface/TrackClassifier.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"

#include "TROOT.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"

class V0Validator : public edm::EDAnalyzer {

public:
  explicit V0Validator(const edm::ParameterSet&);
  ~V0Validator();


private:
  //virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void endJob() ;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&);

  //Quantities that are to be histogrammed
  float K0sGenEta, LamGenEta, K0sGenpT, LamGenpT, K0sGenR, LamGenR;
  float LamGenX, LamGenY, LamGenZ, KsGenX, KsGenY, KsGenZ;
  float K0sCandEta, LamCandEta, K0sCandpT, LamCandpT, K0sCandR, LamCandR;
  unsigned int K0sGenStatus, LamGenStatus, K0sCandStatus, LamCandStatus;
  unsigned int K0sPiCandStatus[2], LamPiCandStatus[2], K0sPiEff[2], LamPiEff[2];

  //Bookkeeping quantities
  int genLam, genK0s, realLamFound, realK0sFound, realLamFoundEff, realK0sFoundEff;
  int lamTracksFound, k0sTracksFound, lamCandFound, k0sCandFound, noTPforK0sCand, noTPforLamCand;

  std::vector<double> etaBins, ptBins;

  edm::Service<TFileService> theDQMstore;

  // kshort eff numerator
  TH1D* ksEffVsR_num;
  TH1D* ksEffVsEta_num;
  TH1D* ksEffVsPt_num;
  TH2D* ksEffVsEtaPt_num;
  TH1D* ksTkEffVsR_num;
  TH1D* ksTkEffVsEta_num;
  TH1D* ksTkEffVsPt_num;
  TH2D* ksTkEffVsEtaPt_num;

  // kshort fake numerator
  TH1D* ksFakeVsR_num;
  TH1D* ksFakeVsEta_num;
  TH1D* ksFakeVsPt_num;
  TH2D* ksFakeVsEtaPt_num;
  TH1D* ksTkFakeVsR_num;
  TH1D* ksTkFakeVsEta_num;
  TH1D* ksTkFakeVsPt_num;
  TH2D* ksTkFakeVsEtaPt_num;

  // kshort fake and eff denominator 
  TH1D* ksFakeVsR_denom;
  TH1D* ksFakeVsEta_denom;
  TH1D* ksFakeVsPt_denom;
  TH2D* ksFakeVsEtaPt_denom;
  TH1D* ksEffVsR_denom;
  TH1D* ksEffVsEta_denom;
  TH1D* ksEffVsPt_denom;
  TH2D* ksEffVsEtaPt_denom;

  // lambda fake and eff denominator 
  TH1D* lamFakeVsR_denom;
  TH1D* lamFakeVsEta_denom;
  TH1D* lamFakeVsPt_denom;
  TH2D* lamFakeVsEtaPt_denom;
  TH1D* lamEffVsR_denom;
  TH1D* lamEffVsEta_denom;
  TH1D* lamEffVsPt_denom;
  TH2D* lamEffVsEtaPt_denom;

  // lambda eff numerator
  TH1D* lamEffVsR_num;
  TH1D* lamEffVsEta_num;
  TH1D* lamEffVsPt_num;
  TH2D* lamEffVsEtaPt_num;
  TH1D* lamTkEffVsR_num;
  TH1D* lamTkEffVsEta_num;
  TH1D* lamTkEffVsPt_num;
  TH2D* lamTkEffVsEtaPt_num;

  // lambda fake numerator
  TH1D* lamFakeVsR_num;
  TH1D* lamFakeVsEta_num;
  TH1D* lamFakeVsPt_num;
  TH2D* lamFakeVsEtaPt_num;
  TH1D* lamTkFakeVsR_num;
  TH1D* lamTkFakeVsEta_num;
  TH1D* lamTkFakeVsPt_num;
  TH2D* lamTkFakeVsEtaPt_num;

  TH1D* ksXResolution;
  TH1D* ksYResolution;
  TH1D* ksZResolution;
  TH1D* ksAbsoluteDistResolution;
  TH1D* lamXResolution;
  TH1D* lamYResolution;
  TH1D* lamZResolution;
  TH1D* lamAbsoluteDistResolution;

  TH1D* nKs;
  TH1D* nLam;

  TH1D* ksCandStatus;
  TH1D* lamCandStatus;

  TH1D* fakeKsMass;
  TH1D* goodKsMass;
  TH1D* fakeLamMass;
  TH1D* goodLamMass;
  TH2D* fakeKsMassPt;
  TH2D* goodKsMassPt;
  TH2D* fakeLamMassPt;
  TH2D* goodLamMassPt;
  TH1D* ksMassAll;
  TH1D* lamMassAll;
  TH2D* ksMassPtAll;
  TH2D* lamMassPtAll;

  TH1D* ksFakeDauRadDist;
  TH1D* lamFakeDauRadDist;

  edm::InputTag trackCollectionTag;
  edm::InputTag vertexCollectionTag;
  edm::InputTag k0sCollectionTag;
  edm::InputTag lamCollectionTag;

  bool isMatchByHitsOrChi2_;
};
