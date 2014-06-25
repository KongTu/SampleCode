// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TNtuple.h>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


//
// class decleration
//

#define PI 3.1416
using namespace std;

class V0para : public edm::EDAnalyzer {
public:
  explicit V0para(const edm::ParameterSet&);
  ~V0para();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
    
    TNtuple* V0para_ks;
    TNtuple* V0para_la;
    TNtuple* V0para_Xi;
    TNtuple* V0para_Omg;
    
    double etaMin_trg_;
    double etaMax_trg_;
    double etaMin_ass_;
    double etaMax_ass_;
    double multMax_;
    double multMin_;
    int bkgFactor_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

V0para::V0para(const edm::ParameterSet& iConfig)
{

  //now do what ever initialization is needed
    etaMin_trg_ = iConfig.getUntrackedParameter<double>("etaMin_trg", -2.4);
    etaMax_trg_ = iConfig.getUntrackedParameter<double>("etaMax_trg", 2.4);
    etaMin_ass_ = iConfig.getUntrackedParameter<double>("etaMin_ass", -2.4);
    etaMax_ass_ = iConfig.getUntrackedParameter<double>("etaMax_ass", 2.4);
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", 220);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", 185);
    bkgFactor_ = iConfig.getUntrackedParameter<int>("bkgFactor", 10);

    
}


V0para::~V0para()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
V0para::analyze(const edm::Event& iEvent, const edm::EventSetup& 
iSetup)
{
    using std::vector;
    using namespace edm;
    
    // select on requirement of valid vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel("offlinePrimaryVertices",vertices);
    double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    if(bestvz < -15.0 || bestvz>15.0) return;

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_ks;
    iEvent.getByLabel("generalV0CandidatesNew","Kshort",v0candidates_ks);
    if(!v0candidates_ks.isValid()) return;
    
    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_la;
    iEvent.getByLabel("generalV0CandidatesNew","Lambda",v0candidates_la);
    if(!v0candidates_la.isValid()) return;
    
    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_Xi;
    iEvent.getByLabel("generalV0CandidatesNew","Xi",v0candidates_Xi);
    if(!v0candidates_Xi.isValid()) return;
    
    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_Omg;
    iEvent.getByLabel("generalV0CandidatesNew","Omega",v0candidates_Omg);
    if(!v0candidates_Omg.isValid()) return;
    
    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle;
    iEvent.getByLabel("dedxTruncated40", dEdxHandle);
    
    for(unsigned it=0; it<v0candidates_ks->size(); ++it){
        
        const reco::VertexCompositeCandidate & trk = (*v0candidates_ks)[it];
        
        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
                    
        const reco::Candidate * d1 = trk.daughter(0);
        const reco::Candidate * d2 = trk.daughter(1);
        
        auto dau1 = d1->get<reco::TrackRef>();
        auto dau2 = d2->get<reco::TrackRef>();
        
        double dedx1 = -999.9;
        double dedx2 = -999.9;
        
        //dEdx,daughter charge, momentum, need p vs dEdx 2D
        if(dEdxHandle->size()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle.product();
            dedx1 = dEdxTrack[dau1].dEdx();
            dedx2 = dEdxTrack[dau2].dEdx();
        }
        
        double p1 = d1->p()*d1->charge();
        double p2 = d2->p()*d2->charge();
        
        //pt,mass
        double eta = trk.eta();
        double pt = trk.pt();
        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        double mass = trk.mass();
        
        if(eta>etaMax_trg_ || eta<etaMin_trg_) continue;
        
        //vertexCovariance 00-xError 11-y 22-z
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
        //secvzError = sqrt(trk.vertexCovariance(2,2)); secvxError = sqrt(trk.vertexCovariance(0,0)); secvyError = sqrt(trk.vertexCovariance(1,1));
        
        //cout<<trk.vertexCovariance(0,1)<<endl;
        
        //trkNHits
        double nhit1 = dau1->numberOfValidHits();
        double nhit2 = dau2->numberOfValidHits();
        
        //DCA
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzbest1 = dau1->dz(bestvtx);
        double dxybest1 = dau1->dxy(bestvtx);
        double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
        double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);
    
        double dzos1 = dzbest1/dzerror1;
        double dxyos1 = dxybest1/dxyerror1;
        
        double dzbest2 = dau2->dz(bestvtx);
        double dxybest2 = dau2->dxy(bestvtx);
        double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
        double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
        
        double dzos2 = dzbest2/dzerror2;
        double dxyos2 = dxybest2/dxyerror2;
        
        //vtxChi2
        double vtxChi2 = trk.vertexChi2();
        
        //PAngle
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);
            
        double agl = cos(secvec.Angle(ptosvec));
        
        //Decay length
        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        
        SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
        SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        
        double dl = ROOT::Math::Mag(distanceVector);
        double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
        
        /*double dl = sqrt((secvx-bestvx)*(secvx-bestvx)+(secvy-bestvy)*(secvy-bestvy)+(secvz-bestvz)*(secvz-bestvz));
        //double dlerror = sqrt((secvx-bestvx)*(secvx-bestvx)*bestvxError*bestvxError+(secvy-bestvy)*(secvy-bestvy)*bestvyError*bestvyError+(secvz-bestvz)*(secvz-bestvz)*bestvzError*bestvzError)/dl;
        double dlerror = sqrt((secvx-bestvx)*(secvx-bestvx)*(secvxError*secvxError+bestvxError*bestvxError)+2.0*(secvy-bestvy)*(secvx-bestvx)*(sxysec+sxybest)+(secvy-bestvy)*(secvy-bestvy)*(secvyError*secvyError+bestvyError*bestvyError)+(secvz-bestvz)*(secvz-bestvz)*(secvzError*secvzError+bestvzError*bestvzError))/dl;*/
        
        double dlos = dl/dlerror;
        
        //Fill
        V0para_ks->Fill(pt,mass,dzbest1,dzerror1,dzbest2,dzerror2,dxybest1,dxyerror1,dxybest2,dxyerror2,nhit1,nhit2,dl,dlerror,agl);
            
    }
    
    
    for(unsigned it=0; it<v0candidates_la->size(); ++it){
        
        const reco::VertexCompositeCandidate & trk = (*v0candidates_la)[it];
        
        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
        
        const reco::Candidate * d1 = trk.daughter(0);
        const reco::Candidate * d2 = trk.daughter(1);
        
        auto dau1 = d1->get<reco::TrackRef>();
        auto dau2 = d2->get<reco::TrackRef>();
        
        double dedx1 = -999.9;
        double dedx2 = -999.9;
        
        //dEdx
        if(dEdxHandle->size()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle.product();
            dedx1 = dEdxTrack[dau1].dEdx();
            dedx2 = dEdxTrack[dau2].dEdx();
        }
        
        double p1 = d1->p()*d1->charge();
        double p2 = d2->p()*d2->charge();
        
        //pt,mass
        double eta = trk.eta();
        double pt = trk.pt();
        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        double mass = trk.mass();
        
        if(eta>etaMax_trg_ || eta<etaMin_trg_) continue;
        
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
        //secvzError = sqrt(trk.vertexCovariance(2,2)); secvxError = sqrt(trk.vertexCovariance(0,0)); secvyError = sqrt(trk.vertexCovariance(1,1));
        
        //trkNHits
        double nhit1 = dau1->numberOfValidHits();
        double nhit2 = dau2->numberOfValidHits();
        
        //DCA
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzbest1 = dau1->dz(bestvtx);
        double dxybest1 = dau1->dxy(bestvtx);
        double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
        double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);
        
        double dzos1 = dzbest1/dzerror1;
        double dxyos1 = dxybest1/dxyerror1;
        
        double dzbest2 = dau2->dz(bestvtx);
        double dxybest2 = dau2->dxy(bestvtx);
        double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
        double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
        
        double dzos2 = dzbest2/dzerror2;
        double dxyos2 = dxybest2/dxyerror2;
        
        //vtxChi2
        double vtxChi2 = trk.vertexChi2();
        
        //PAngle
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);
        
        double agl = cos(secvec.Angle(ptosvec));
        
        //Decay length
        //3D
        /*double dl = sqrt((secvx-bestvx)*(secvx-bestvx)+(secvy-bestvy)*(secvy-bestvy)+(secvz-bestvz)*(secvz-bestvz));
        //double dlerror = sqrt((secvx-bestvx)*(secvx-bestvx)*bestvxError*bestvxError+(secvy-bestvy)*(secvy-bestvy)*bestvyError*bestvyError+(secvz-bestvz)*(secvz-bestvz)*bestvzError*bestvzError)/dl;
        
        double dlerror = sqrt((secvx-bestvx)*(secvx-bestvx)*(secvxError*secvxError+bestvxError*bestvxError)+2.0*(secvy-bestvy)*(secvx-bestvx)*(sxysec+sxybest)+(secvy-bestvy)*(secvy-bestvy)*(secvyError*secvyError+bestvyError*bestvyError)+(secvz-bestvz)*(secvz-bestvz)*(secvzError*secvzError+bestvzError*bestvzError))/dl;*/
        
        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        
        SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
        SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        
        double dl = ROOT::Math::Mag(distanceVector);
        double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
        
        double dlos = dl/dlerror;
        
        //Fill
        V0para_la->Fill(pt,mass,dzbest1,dzerror1,dzbest2,dzerror2,dxybest1,dxyerror1,dxybest2,dxyerror2,nhit1,nhit2,dl,dlerror,agl);
        
    }

    for(unsigned it=0; it<v0candidates_Xi->size(); ++it){
        
        const reco::VertexCompositeCandidate & trk = (*v0candidates_Xi)[it];
        
        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
        
        const reco::Candidate * d2 = trk.daughter(1);
        
        auto dau2 = d2->get<reco::TrackRef>();
        
        double dedx2 = -999.9;
        
        //dEdx,daughter charge, momentum, need p vs dEdx 2D
        if(dEdxHandle->size()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle.product();
            dedx2 = dEdxTrack[dau2].dEdx();
        }
        
        double p2 = d2->p()*d2->charge();
        
        //pt,mass
        double eta = trk.eta();
        double pt = trk.pt();
        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        double mass = trk.mass();
        
        if(eta>etaMax_trg_ || eta<etaMin_trg_) continue;
        
        //vertexCovariance 00-xError 11-y 22-z
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
        //secvzError = sqrt(trk.vertexCovariance(2,2)); secvxError = sqrt(trk.vertexCovariance(0,0)); secvyError = sqrt(trk.vertexCovariance(1,1));
        
        //cout<<trk.vertexCovariance(0,1)<<endl;
        
        //trkNHits
        double nhit2 = dau2->numberOfValidHits();
        
        //DCA
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzbest2 = dau2->dz(bestvtx);
        double dxybest2 = dau2->dxy(bestvtx);
        double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
        double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
        
        double dzos2 = dzbest2/dzerror2;
        double dxyos2 = dxybest2/dxyerror2;
        
        //vtxChi2
        double vtxChi2 = trk.vertexChi2();
        
        //PAngle
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);
        
        double agl = cos(secvec.Angle(ptosvec));
        
        //Decay length
        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        
        SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
        SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        
        double dl = ROOT::Math::Mag(distanceVector);
        double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
        
        /*double dl = sqrt((secvx-bestvx)*(secvx-bestvx)+(secvy-bestvy)*(secvy-bestvy)+(secvz-bestvz)*(secvz-bestvz));
         //double dlerror = sqrt((secvx-bestvx)*(secvx-bestvx)*bestvxError*bestvxError+(secvy-bestvy)*(secvy-bestvy)*bestvyError*bestvyError+(secvz-bestvz)*(secvz-bestvz)*bestvzError*bestvzError)/dl;
         double dlerror = sqrt((secvx-bestvx)*(secvx-bestvx)*(secvxError*secvxError+bestvxError*bestvxError)+2.0*(secvy-bestvy)*(secvx-bestvx)*(sxysec+sxybest)+(secvy-bestvy)*(secvy-bestvy)*(secvyError*secvyError+bestvyError*bestvyError)+(secvz-bestvz)*(secvz-bestvz)*(secvzError*secvzError+bestvzError*bestvzError))/dl;*/
        
        double dlos = dl/dlerror;
        
        //Fill
        V0para_Xi->Fill(pt,mass,eta,dzbest2,dzerror2,dxybest2,dxyerror2,vtxChi2,nhit2,dl,dlerror,agl);
        
    }
    
    for(unsigned it=0; it<v0candidates_Omg->size(); ++it){
        
        const reco::VertexCompositeCandidate & trk = (*v0candidates_Omg)[it];
        
        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
        
        const reco::Candidate * d2 = trk.daughter(1);
        
        auto dau2 = d2->get<reco::TrackRef>();
        
        double dedx2 = -999.9;
        
        //dEdx,daughter charge, momentum, need p vs dEdx 2D
        if(dEdxHandle->size()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle.product();
            dedx2 = dEdxTrack[dau2].dEdx();
        }
        
        double p2 = d2->p()*d2->charge();
        
        //pt,mass
        double eta = trk.eta();
        double pt = trk.pt();
        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        double mass = trk.mass();
        
        if(eta>etaMax_trg_ || eta<etaMin_trg_) continue;
        
        //vertexCovariance 00-xError 11-y 22-z
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
        //secvzError = sqrt(trk.vertexCovariance(2,2)); secvxError = sqrt(trk.vertexCovariance(0,0)); secvyError = sqrt(trk.vertexCovariance(1,1));
        
        //cout<<trk.vertexCovariance(0,1)<<endl;
        
        //trkNHits
        double nhit2 = dau2->numberOfValidHits();
        
        //DCA
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzbest2 = dau2->dz(bestvtx);
        double dxybest2 = dau2->dxy(bestvtx);
        double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
        double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
        
        double dzos2 = dzbest2/dzerror2;
        double dxyos2 = dxybest2/dxyerror2;
        
        //vtxChi2
        double vtxChi2 = trk.vertexChi2();
        
        //PAngle
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);
        
        double agl = cos(secvec.Angle(ptosvec));
        
        //Decay length
        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        
        SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
        SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        
        double dl = ROOT::Math::Mag(distanceVector);
        double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
        
        /*double dl = sqrt((secvx-bestvx)*(secvx-bestvx)+(secvy-bestvy)*(secvy-bestvy)+(secvz-bestvz)*(secvz-bestvz));
         //double dlerror = sqrt((secvx-bestvx)*(secvx-bestvx)*bestvxError*bestvxError+(secvy-bestvy)*(secvy-bestvy)*bestvyError*bestvyError+(secvz-bestvz)*(secvz-bestvz)*bestvzError*bestvzError)/dl;
         double dlerror = sqrt((secvx-bestvx)*(secvx-bestvx)*(secvxError*secvxError+bestvxError*bestvxError)+2.0*(secvy-bestvy)*(secvx-bestvx)*(sxysec+sxybest)+(secvy-bestvy)*(secvy-bestvy)*(secvyError*secvyError+bestvyError*bestvyError)+(secvz-bestvz)*(secvz-bestvz)*(secvzError*secvzError+bestvzError*bestvzError))/dl;*/
        
        double dlos = dl/dlerror;
        
        //Fill
        V0para_Omg->Fill(pt,mass,eta,dzbest2,dzerror2,dxybest2,dxyerror2,vtxChi2,nhit2,dl,dlerror,agl);
        
    }

}


// ------------ method called once each job just before starting event
//loop  ------------
void 
V0para::beginJob()
{
    edm::Service<TFileService> fs;
        
    TH1D::SetDefaultSumw2();
    
    V0para_ks = fs->make< TNtuple>("V0para_ks","V0para_ks","pt:mass:trkDCA1z:trkDCAsigma1z:trkDCA2z:trkDCAsigma2z:trkDCA1xy:trkDCAsigma1xy:trkDCA2xy:trkDCAsigma2xy:trkNHits1:trkNHits2:L:Lsigma:PAngle");
    V0para_la = fs->make< TNtuple>("V0para_la","V0para_la","pt:mass:trkDCA1z:trkDCAsigma1z:trkDCA2z:trkDCAsigma2z:trkDCA1xy:trkDCAsigma1xy:trkDCA2xy:trkDCAsigma2xy:trkNHits1:trkNHits2:L:Lsigma:PAngle");
    V0para_Xi = fs->make< TNtuple>("V0para_Xi","V0para_Xi","pt:mass:eta:trkDCAz:trkDCAsigmaz:trkDCAxy:trkDCAsigmaxy:vtxChi2:trkNHits:L:Lsigma:PAngle");
    V0para_Omg = fs->make< TNtuple>("V0para_Omg","V0para_Omg","pt:mass:eta:trkDCAz:trkDCAsigmaz:trkDCAxy:trkDCAsigmaxy:vtxChi2:trkNHits:L:Lsigma:PAngle");
}

// ------------ method called once each job just after ending the event
//loop  ------------
void 
V0para::endJob() {
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(V0para);


