// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>

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

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//
// class decleration
//

#define PI 3.1416

/*#define effksl_0(x) ((9.23396e+06*TMath::Power(x,1.00849e+01)*TMath::Exp(3.74519*x*x-1.89234e+01*x+9.53588))/(7.04699e+07*TMath::Power(x,2.69212)*TMath::Exp(1.76962*x*x-9.10509*x+2.04396)))

#define effksh_0(x) ((1.44634e+01*TMath::Power(x,4.03640e-01)*TMath::Exp(1.23220e-01*x*x-2.14135*x+9.74838))/(2.84992e+05*TMath::Power(x,-2.05794)*TMath::Exp(3.43241e-02*x*x-8.38940e-01*x+1.00439)))

#define effksl_1(x) ((6.26825e+06*TMath::Power(x,2.48359e+01)*TMath::Exp(1.41014e+01*x*x-5.44703e+01*x+3.51496e+01))/(7.89456e+07*TMath::Power(x,3.19695)*TMath::Exp(1.96620*x*x-9.91973*x+2.58793)))

#define effksh_1(x) ((2.24818e+02*TMath::Power(x,-1.37187)*TMath::Exp(-9.02356e-03*x*x-7.32575e-01*x+6.03275))/(5.17714e+05*TMath::Power(x,-1.78945)*TMath::Exp(5.02927e-02*x*x-1.01589*x+6.06828e-01)))

#define effksl_2(x) ((7.93665e+07*TMath::Power(x,2.87601e+01)*TMath::Exp(1.53928e+01*x*x-6.07704e+01*x+3.74924e+01))/(7.32916e+05*TMath::Power(x,3.84557)*TMath::Exp(2.30522*x*x-1.11115e+01*x+8.09969)))

#define effksh_2(x) ((7.11878e+01*TMath::Power(x,1.07130)*TMath::Exp(1.13515e-01*x*x-2.24369*x+8.21807))/(1.13902e+02*TMath::Power(x,-1.63500)*TMath::Exp(5.05541e-02*x*x-1.05330*x+9.04806)))

#define effksl_3(x) ((7.54402e+08*TMath::Power(x,2.71092e+01)*TMath::Exp(1.39218e+01*x*x-5.60112e+01*x+3.19798e+01))/(6.24416e+07*TMath::Power(x,3.81176)*TMath::Exp(2.23869*x*x-1.09532e+01*x+3.59543)))

#define effksh_3(x) ((9.11234e-04*TMath::Power(x,8.13692e-01)*TMath::Exp(9.62341e-02*x*x-2.05879*x+1.93640e+01))/(2.97772e+02*TMath::Power(x,-1.44430)*TMath::Exp(6.78137e-02*x*x-1.23135*x+8.26981)))

#define effksl_4(x) ((1.21686e+06*TMath::Power(x,2.20006e+01)*TMath::Exp(1.10768e+01*x*x-4.56693e+01*x+3.11482e+01))/(1.07805e+08*TMath::Power(x,3.43171)*TMath::Exp(2.17143*x*x-1.05289e+01*x+2.77827)))

#define effksh_4(x) ((6.17963*TMath::Power(x,-1.24846)*TMath::Exp(3.59753e-02*x*x-1.05332*x+1.00256e+01))/(2.84028e+03*TMath::Power(x,-1.97138)*TMath::Exp(3.49249e-02*x*x-8.56336e-01*x+5.77800)))

#define effksl_5(x) ((6.58219e+05*TMath::Power(x,1.13154e+01)*TMath::Exp(4.85838*x*x-2.22038e+01*x+1.45414e+01))/(1.36541e+08*TMath::Power(x,2.77573)*TMath::Exp(1.77919*x*x-9.15881*x+1.61922)))

#define effksh_5(x) ((1.12556e+05*TMath::Power(x,-2.68907)*TMath::Exp(-7.81469e-02*x*x+2.53778e-01*x-8.74728e-01))/(1.14665e+03*TMath::Power(x,-1.95576)*TMath::Exp(4.44103e-02*x*x-9.41693e-01*x+6.81057)))

#define efflal_0(x) ((4.52697e+03*TMath::Power(x,1.39899e+01)*TMath::Exp(1.51525*x*x-1.45547e+01*x+1.29369e+01))/(7.87202e+03*TMath::Power(x,3.01668)*TMath::Exp(1.23042*x*x-7.54580*x+9.70611)))

#define efflah_0(x) ((1.34124*TMath::Power(x,4.40016)*TMath::Exp(1.74940e-01*x*x-3.57612*x+1.11304e+01))/(3.65140*TMath::Power(x,-1.16995)*TMath::Exp(6.49969e-02*x*x-1.23244*x+1.22890e+01)))

#define efflal_1(x) ((5.11070e+04*TMath::Power(x,2.23592e+01)*TMath::Exp(4.28951*x*x-2.85957e+01*x+2.15102e+01))/(2.43927e+04*TMath::Power(x,3.51057)*TMath::Exp(1.29094*x*x-7.93633*x+8.82477)))

#define efflah_1(x) ((1.01462e-01*TMath::Power(x,2.36488)*TMath::Exp(3.33604e-02*x*x-1.97685*x+1.22692e+01))/(3.60202e-01*TMath::Power(x,-1.39632e-01)*TMath::Exp(9.56626e-02*x*x-1.79088*x+1.49352e+01)))

#define efflal_2(x) ((2.52993*TMath::Power(x,2.82304e+01)*TMath::Exp(5.69436*x*x-3.67388e+01*x+3.79121e+01))/(3.04280e+06*TMath::Power(x,3.65330)*TMath::Exp(1.21409*x*x-7.66759*x+3.64783)))

#define efflah_2(x) ((5.09173e-03*TMath::Power(x,-9.40931e-01)*TMath::Exp(-5.34324e-02*x*x-1.42480e-01*x+1.41216e+01))/(7.73082e+01*TMath::Power(x,-5.44934e-01)*TMath::Exp(5.25787e-02*x*x-1.36246*x+9.14696)))

#define efflal_3(x) ((2.93101e+01*TMath::Power(x,2.53514e+01)*TMath::Exp(5.14362*x*x-3.30728e+01*x+3.23405e+01))/(6.57499*TMath::Power(x,3.64908)*TMath::Exp(1.20900*x*x-7.64879*x+1.66944e+01)))

#define efflah_3(x) ((3.80077e+08*TMath::Power(x,-1.85514)*TMath::Exp(-8.05000e-02*x*x+3.20202e-01*x-1.10763e+01))/(1.12322e+02*TMath::Power(x,9.66739e-01)*TMath::Exp(1.33976e-01*x*x-2.38345*x+9.45469)))

#define efflal_4(x) ((7.12785*TMath::Power(x,2.96440e+01)*TMath::Exp(6.30591*x*x-3.95837e+01*x+3.93720e+01))/(9.77250e+06*TMath::Power(x,3.55832)*TMath::Exp(1.26341*x*x-7.89410*x+2.89594)))

#define efflah_4(x) ((1.02750e+04*TMath::Power(x,1.24480)*TMath::Exp(4.78918e-02*x*x-1.77726*x+1.19110))/(5.67745e+01*TMath::Power(x,-7.61973e-02)*TMath::Exp(8.45505e-02*x*x-1.75756*x+9.87205)))

#define efflal_5(x) ((1.44598e+01*TMath::Power(x,9.52920)*TMath::Exp(5.61404e-01*x*x-8.42608*x+1.36199e+01))/(1.83734e+07*TMath::Power(x,3.04202)*TMath::Exp(1.17260*x*x-7.38097*x+1.99693)))

#define efflah_5(x) ((2.18290e+01*TMath::Power(x,2.20048)*TMath::Exp(1.00320e-01*x*x-2.50334*x+8.37433))/(1.11728e+02*TMath::Power(x,9.53453e-02)*TMath::Exp(1.06586e-01*x*x-2.00067*x+9.52131)))*/

#define effkss(x) ((1.09925e+11*TMath::Power(x,1.50511e+01)*TMath::Exp(6.58074*x*x-2.94487e+01*x+9.72194))/(6.71504e+11*TMath::Power(x,3.19081)*TMath::Exp(1.97717*x*x-9.90447*x-4.65781)))

#define effksb(x) ((6.47559e+02*TMath::Power(x,2.95369e-01)*TMath::Exp(9.18237e-02*x*x-1.89678*x+7.63891))/(1.20601e+01*TMath::Power(x,-1.86165)*TMath::Exp(4.17033e-02*x*x-9.39659e-01*x+1.30386e+01)))

#define efflas(x) ((4.76443e+06*TMath::Power(x,1.62889e+01)*TMath::Exp(2.74004*x*x-1.97581e+01*x+1.16309e+01))/(5.30771e+08*TMath::Power(x,3.59273)*TMath::Exp(1.50703*x*x-8.45701*x+9.43797e-01)))

#define efflab(x) ((3.86297e-01*TMath::Power(x,1.91207)*TMath::Exp(8.37588e-02*x*x-2.15583*x+1.33689e+01))/(7.01220*TMath::Power(x,-4.80662e-01)*TMath::Exp(7.33837e-02*x*x-1.53854*x+1.35978e+01)))

using namespace std;

class DemoAnalyzerAll220 : public edm::EDAnalyzer {
public:
    explicit DemoAnalyzerAll220(const edm::ParameterSet&);
    ~DemoAnalyzerAll220();
    
    
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    // ----------member data ---------------------------
    
    TH1D* hMult;
    
    TH1D* hMass_ks[13];
    TH1D* hMass_la[13];
    
    TH1D* hKET_ks[13];
    TH1D* hKET_la[13];
    
    TH1D* hPt_ks[13];
    TH1D* hPt_la[13];
    
    TH1D* hMult_ks[13];
    TH2D* hSignal_ks[13];
    TH2D* hBackground_ks[13];
    
    TH1D* hKET_ks_bkg[13];
    TH1D* hKET_la_bkg[13];
    
    TH1D* hPt_ks_bkg[13];
    TH1D* hPt_la_bkg[13];
    
    TH1D* hMult_ks_bkg[13];
    TH2D* hSignal_ks_bkg[13];
    TH2D* hBackground_ks_bkg[13];
    
    TH1D* hMult_ass;
    
    TH1D* hMult_la[13];
    TH2D* hSignal_la[13];
    TH2D* hBackground_la[13];
    
    vector<TVector3> *pVect_trg_ks[13];
    vector< vector<TVector3> > *pVectVect_trg_ks[13];
    vector<TVector3> *pVect_trg_la[13];
    vector< vector<TVector3> > *pVectVect_trg_la[13];
    
    vector<TVector3> *pVect_dau_ks[13];
    vector<TVector3> *pVect_dau_la[13];
    
    TH1D* hMult_la_bkg[13];
    TH2D* hSignal_la_bkg[13];
    TH2D* hBackground_la_bkg[13];
    
    vector<TVector3> *pVect_trg_ks_bkg[13];
    vector< vector<TVector3> > *pVectVect_trg_ks_bkg[13];
    vector<TVector3> *pVect_trg_la_bkg[13];
    vector< vector<TVector3> > *pVectVect_trg_la_bkg[13];
    
    vector<TVector3> *pVect_dau_ks_bkg[13];
    vector<TVector3> *pVect_dau_la_bkg[13];
    
    vector<TVector3> *pVect_ass;
    vector< vector<TVector3> > *pVectVect_ass;
    vector<double> *zvtxVect;
    
    double ptcut[14];
    double mass_ks_max[13];
    double mass_ks_min[13];
    double mass_la_max[13];
    double mass_la_min[13];
    
    double etaMin_trg_;
    double etaMax_trg_;
    double etaMin_ass_;
    double etaMax_ass_;
    double ptMin_ass_;
    double ptMax_ass_;
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

DemoAnalyzerAll220::DemoAnalyzerAll220(const edm::ParameterSet& iConfig)
{
    
    //now do what ever initialization is needed
    etaMin_trg_ = iConfig.getUntrackedParameter<double>("etaMin_trg", -2.4);
    etaMax_trg_ = iConfig.getUntrackedParameter<double>("etaMax_trg", 2.4);
    etaMin_ass_ = iConfig.getUntrackedParameter<double>("etaMin_ass", -2.4);
    etaMax_ass_ = iConfig.getUntrackedParameter<double>("etaMax_ass", 2.4);
    ptMin_ass_ = iConfig.getUntrackedParameter<double>("ptMin_ass", 0.3);
    ptMax_ass_ = iConfig.getUntrackedParameter<double>("ptMax_ass", 3.0);
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", 220);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", 185);
    bkgFactor_ = iConfig.getUntrackedParameter<int>("bkgFactor", 10);
    
}


DemoAnalyzerAll220::~DemoAnalyzerAll220()
{
    
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
DemoAnalyzerAll220::analyze(const edm::Event& iEvent, const edm::EventSetup&
                         iSetup)
{
    using namespace edm;
    
    double ptcut[14] = {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,9.0,12.0};
    /*double etacut[7] = {-2.4,-1.6,-0.8,0,0.8,1.6,2.4};
    double kseff[6] = {1.2,1.2,1.2,1.2,1.2,1.2};
    double laeff[6] = {2,2,2,2,2,1.8};*/
    
    double mass_ks_max[13] = {0.516999,0.515842,0.514229,0.513062,0.511695,0.510907,0.510989,0.510853,0.510397,0.511036,0.517516,0.518819,0.519159};
    double mass_ks_min[13] = {0.477297,0.478353,0.480659,0.482186,0.483628,0.484326,0.484132,0.484217,0.484634,0.484022,0.477581,0.476156,0.476372};
    double mass_la_max[13] = {1.1228,1.1228,1.1228,1.12282,1.12212,1.12202,1.1221,1.12184,1.12151,1.12083,1.12287,1.12342,1.12335};
    double mass_la_min[13] = {1.10979,1.10979,1.10979,1.10995,1.10986,1.10981,1.10974,1.10994,1.11021,1.1108,1.10867,1.10814,1.10749};
    
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
    
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel("generalTracks", tracks);
    
    for(int i=0;i<13;i++)
    {
        pVect_trg_ks[i] = new vector<TVector3>;
        pVect_trg_la[i] = new vector<TVector3>;
        pVect_dau_ks[i] = new vector<TVector3>;
        pVect_dau_la[i] = new vector<TVector3>;
        pVect_trg_ks_bkg[i] = new vector<TVector3>;
        pVect_trg_la_bkg[i] = new vector<TVector3>;
        pVect_dau_ks_bkg[i] = new vector<TVector3>;
        pVect_dau_la_bkg[i] = new vector<TVector3>;
    }
    
    pVect_ass = new vector<TVector3>;
    

    //track selection
    int nMult_ass_good = 0;
    for(unsigned it=0; it<tracks->size(); ++it){
        
        const reco::Track & trk = (*tracks)[it];
        
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;
        
        double eta = trk.eta();
        double pt  = trk.pt();
        
        if(fabs(eta)>2.4) continue;
        if(pt<=0.4) continue;
        nMult_ass_good++;
    }
    hMult->Fill(nMult_ass_good);
    
    if(nMult_ass_good<multMax_ && nMult_ass_good>=multMin_){
        //loop over tracks
        for(unsigned it=0; it<v0candidates_ks->size(); ++it){
            
            const reco::VertexCompositeCandidate & trk = (*v0candidates_ks)[it];
            
            double secvz=-999.9, secvx=-999.9, secvy=-999.9;
            //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
            
            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
            //secvzError = trk.vertexCovariance(2,2); secvxError = trk.vertexCovariance(0,0); secvyError = trk.vertexCovariance(1,1);
            
            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();
            double mass = trk.mass();
            double px = trk.px();
            double py = trk.py();
            double pz = trk.pz();
            
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(px,py,pz);
            
            double agl = cos(secvec.Angle(ptosvec));
            
            if(agl<=0.999) continue;
            
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            
            SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            double dl = ROOT::Math::Mag(distanceVector);
            double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
            
            double dlos = dl/dlerror;
            
            if(dlos<=5) continue;
            
            //efficiency
            double effks = 1.0;
            
            if(pt<1.2)
            {
                effks = effkss(pt);
            }
            else
            {
                effks = effksb(pt);
            }
            
            /*if(eta>etacut[0] && eta<etacut[1])
            {
                if(pt<kseff[0])
                {
                    effks = effksl_0(pt);
                }
                else
                {
                    effks = effksh_0(pt);
                }
            }
            if(eta>etacut[1] && eta<etacut[2])
            {
                if(pt<kseff[1])
                {
                    effks = effksl_1(pt);
                }
                else
                {
                    effks = effksh_1(pt);
                }
            }
            if(eta>etacut[2] && eta<etacut[3])
            {
                if(pt<kseff[2])
                {
                    effks = effksl_2(pt);
                }
                else
                {
                    effks = effksh_2(pt);
                }
            }
            if(eta>etacut[3] && eta<etacut[4])
            {
                if(pt<kseff[3])
                {
                    effks = effksl_3(pt);
                }
                else
                {
                    effks = effksh_3(pt);
                }
            }
            if(eta>etacut[4] && eta<etacut[5])
            {
                if(pt<kseff[4])
                {
                    effks = effksl_4(pt);
                }
                else
                {
                    effks = effksh_4(pt);
                }
            }
            if(eta>etacut[5] && eta<etacut[6])
            {
                if(pt<kseff[5])
                {
                    effks = effksl_5(pt);
                }
                else
                {
                    effks = effksh_5(pt);
                }
            }*/
            
            const reco::Candidate * dau1 = trk.daughter(0);
            const reco::Candidate * dau2 = trk.daughter(1);
            
            double eta_dau1 = dau1->eta();
            double phi_dau1 = dau1->phi();
            double pt_dau1 = dau1->pt();
            
            double eta_dau2 = dau2->eta();
            double phi_dau2 = dau2->phi();
            double pt_dau2 = dau2->pt();
            
            TVector3 pvector;
            pvector.SetPtEtaPhi(pt,eta,phi);
            
            TVector3 pvector_dau1;
            pvector_dau1.SetPtEtaPhi(pt_dau1,eta_dau1,phi_dau1);
            
            TVector3 pvector_dau2;
            pvector_dau2.SetPtEtaPhi(pt_dau2,eta_dau2,phi_dau2);
            
            for(int i=0;i<13;i++)
            {
                if(trk.eta()<=etaMax_trg_ && trk.eta()>=etaMin_trg_ && trk.pt()<=ptcut[i+1] && trk.pt()>=ptcut[i]){
                    hMass_ks[i]->Fill(mass);
                    if(trk.mass()<=mass_ks_max[i] && trk.mass()>=mass_ks_min[i]){
                        pVect_trg_ks[i]->push_back(pvector);
                        pVect_dau_ks[i]->push_back(pvector_dau1);
                        pVect_dau_ks[i]->push_back(pvector_dau2);
                        hPt_ks[i]->Fill(pt,1.0/effks);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_ks[i]->Fill(KET,1.0/effks);
                    }
                    if((trk.mass()<=0.46 && trk.mass()>=0.43) || (trk.mass()<=0.565 && trk.mass()>=0.54)){
                        pVect_trg_ks_bkg[i]->push_back(pvector);
                        pVect_dau_ks_bkg[i]->push_back(pvector_dau1);
                        pVect_dau_ks_bkg[i]->push_back(pvector_dau2);
                        hPt_ks_bkg[i]->Fill(pt,1.0/effks);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_ks_bkg[i]->Fill(KET,1.0/effks);
                    }
                }
            }
        }
        
        for(unsigned it=0; it<v0candidates_la->size(); ++it){
            
            const reco::VertexCompositeCandidate & trk = (*v0candidates_la)[it];
            
            double secvz=-999.9, secvx=-999.9, secvy=-999.9;
            //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
            
            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
            //secvzError = trk.vertexCovariance(2,2); secvxError = trk.vertexCovariance(0,0); secvyError = trk.vertexCovariance(1,1);
            
            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();
            double mass = trk.mass();
            double px = trk.px();
            double py = trk.py();
            double pz = trk.pz();
            
            TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            TVector3 secvec(px,py,pz);
            
            double agl = cos(secvec.Angle(ptosvec));
            
            if(agl<=0.999) continue;
            
            typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
            typedef ROOT::Math::SVector<double, 3> SVector3;
            
            SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
            SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
            
            double dl = ROOT::Math::Mag(distanceVector);
            double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
            
            double dlos = dl/dlerror;
            
            if(dlos<=5) continue;
            
            //efficiency
            double effla = 1.0;
            
            if(pt<1.6)
            {
                effla = efflas(pt);
            }
            else
            {
                effla = efflab(pt);
            }
            
            /*if(eta>etacut[0] && eta<etacut[1])
            {
                if(pt<laeff[0])
                {
                    effla = efflal_0(pt);
                }
                else
                {
                    effla = efflah_0(pt);
                }
            }
            if(eta>etacut[1] && eta<etacut[2])
            {
                if(pt<laeff[1])
                {
                    effla = efflal_1(pt);
                }
                else
                {
                    effla = efflah_1(pt);
                }
            }
            if(eta>etacut[2] && eta<etacut[3])
            {
                if(pt<laeff[2])
                {
                    effla = efflal_2(pt);
                }
                else
                {
                    effla = efflah_2(pt);
                }
            }
            if(eta>etacut[3] && eta<etacut[4])
            {
                if(pt<laeff[3])
                {
                    effla = efflal_3(pt);
                }
                else
                {
                    effla = efflah_3(pt);
                }
            }
            if(eta>etacut[4] && eta<etacut[5])
            {
                if(pt<laeff[4])
                {
                    effla = efflal_4(pt);
                }
                else
                {
                    effla = efflah_4(pt);
                }
            }
            if(eta>etacut[5] && eta<etacut[6])
            {
                if(pt<laeff[5])
                {
                    effla = efflal_5(pt);
                }
                else
                {
                    effla = efflah_5(pt);
                }
            }*/

            
            const reco::Candidate * dau1 = trk.daughter(0);
            const reco::Candidate * dau2 = trk.daughter(1);
            
            double eta_dau1 = dau1->eta();
            double phi_dau1 = dau1->phi();
            double pt_dau1 = dau1->pt();
            
            double eta_dau2 = dau2->eta();
            double phi_dau2 = dau2->phi();
            double pt_dau2 = dau2->pt();
            
            
            TVector3 pvector;
            pvector.SetPtEtaPhi(pt,eta,phi);
            
            TVector3 pvector_dau1;
            pvector_dau1.SetPtEtaPhi(pt_dau1,eta_dau1,phi_dau1);
            
            TVector3 pvector_dau2;
            pvector_dau2.SetPtEtaPhi(pt_dau2,eta_dau2,phi_dau2);
            
            for(int i=0;i<13;i++)
            {
                if(trk.eta()<=etaMax_trg_ && trk.eta()>=etaMin_trg_ && trk.pt()<=ptcut[i+1] && trk.pt()>=ptcut[i]){
                    hMass_la[i]->Fill(mass);
                    if(trk.mass()<=mass_la_max[i] && trk.mass()>=mass_la_min[i]){
                        pVect_trg_la[i]->push_back(pvector);
                        pVect_dau_la[i]->push_back(pvector_dau1);
                        pVect_dau_la[i]->push_back(pvector_dau2);
                        hPt_la[i]->Fill(pt,1.0/effla);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_la[i]->Fill(KET,1.0/effla);
                    }
                    if(trk.mass()<=1.16 && trk.mass()>=1.13){
                        pVect_trg_la_bkg[i]->push_back(pvector);
                        pVect_dau_la_bkg[i]->push_back(pvector_dau1);
                        pVect_dau_la_bkg[i]->push_back(pvector_dau2);
                        hPt_la_bkg[i]->Fill(pt,1.0/effla);
                        double KET = sqrt(mass*mass + pt*pt) - mass;
                        hKET_la_bkg[i]->Fill(KET,1.0/effla);
                    }
                }
            }
        }
        
        
        for(unsigned it=0; it<tracks->size(); ++it){
            
            const reco::Track & trk = (*tracks)[it];
            
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double dzvtx = trk.dz(bestvtx);
            double dxyvtx = trk.dxy(bestvtx);
            double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
            double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
            
            if(!trk.quality(reco::TrackBase::highPurity)) continue;
            if(fabs(trk.ptError())/trk.pt()>0.10) continue;
            if(fabs(dzvtx/dzerror) > 3) continue;
            if(fabs(dxyvtx/dxyerror) > 3) continue;
            
            
            double eta = trk.eta();
            double phi = trk.phi();
            double pt  = trk.pt();
            
            TVector3 pvector;
            pvector.SetPtEtaPhi(pt,eta,phi);
            if(trk.eta()<=etaMax_ass_ && trk.eta()>=etaMin_ass_ && trk.pt()<=ptMax_ass_ && trk.pt()>=ptMin_ass_) pVect_ass->push_back(pvector);
        }
        
        //Calculating signal
        int nMult_ass = (int)pVect_ass->size();
        hMult_ass->Fill(nMult_ass);
        
        for(int i=0; i<13; i++)
        {
            int nMult_trg_ks = (int)pVect_trg_ks[i]->size();
            int nMult_trg_la = (int)pVect_trg_la[i]->size();
            hMult_ks[i]->Fill(nMult_trg_ks);
            hMult_la[i]->Fill(nMult_trg_la);
            for(int ntrg=0;ntrg<nMult_trg_ks;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg_ks[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double effks = 1.0;
                
                if(pt_trg<1.2)
                {
                    effks = effkss(pt_trg);
                }
                else
                {
                    effks = effksb(pt_trg);
                }

                
                /*if(eta_trg>etacut[0] && eta_trg<etacut[1])
                {
                    if(pt_trg<kseff[0])
                    {
                        effks = effksl_0(pt_trg);
                    }
                    else
                    {
                        effks = effksh_0(pt_trg);
                    }
                }
                if(eta_trg>etacut[1] && eta_trg<etacut[2])
                {
                    if(pt_trg<kseff[1])
                    {
                        effks = effksl_1(pt_trg);
                    }
                    else
                    {
                        effks = effksh_1(pt_trg);
                    }
                }
                if(eta_trg>etacut[2] && eta_trg<etacut[3])
                {
                    if(pt_trg<kseff[2])
                    {
                        effks = effksl_2(pt_trg);
                    }
                    else
                    {
                        effks = effksh_2(pt_trg);
                    }
                }
                if(eta_trg>etacut[3] && eta_trg<etacut[4])
                {
                    if(pt_trg<kseff[3])
                    {
                        effks = effksl_3(pt_trg);
                    }
                    else
                    {
                        effks = effksh_3(pt_trg);
                    }
                }
                if(eta_trg>etacut[4] && eta_trg<etacut[5])
                {
                    if(pt_trg<kseff[4])
                    {
                        effks = effksl_4(pt_trg);
                    }
                    else
                    {
                        effks = effksh_4(pt_trg);
                    }
                }
                if(eta_trg>etacut[5] && eta_trg<etacut[6])
                {
                    if(pt_trg<kseff[5])
                    {
                        effks = effksl_5(pt_trg);
                    }
                    else
                    {
                        effks = effksh_5(pt_trg);
                    }
                }*/
                
                
                TVector3 pvector_trg_dau1 = (*pVect_dau_ks[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();
                
                TVector3 pvector_trg_dau2 = (*pVect_dau_ks[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();
                
                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    
                    if(eta_ass==eta_trg_dau1 && phi_ass==phi_trg_dau1) continue;
                    if(eta_ass==eta_trg_dau2 && phi_ass==phi_trg_dau2) continue;
                    
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;
                    
                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_ks[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_ks);
                }
            }
            
            for(int ntrg=0;ntrg<nMult_trg_la;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg_la[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double effla = 1.0;
                
                if(pt_trg<1.6)
                {
                    effla = efflas(pt_trg);
                }
                else
                {
                    effla = efflab(pt_trg);
                }
                
                /*if(eta_trg>etacut[0] && eta_trg<etacut[1])
                {
                    if(pt_trg<laeff[0])
                    {
                        effla = efflal_0(pt_trg);
                    }
                    else
                    {
                        effla = efflah_0(pt_trg);
                    }
                }
                if(eta_trg>etacut[1] && eta_trg<etacut[2])
                {
                    if(pt_trg<laeff[1])
                    {
                        effla = efflal_1(pt_trg);
                    }
                    else
                    {
                        effla = efflah_1(pt_trg);
                    }
                }
                if(eta_trg>etacut[2] && eta_trg<etacut[3])
                {
                    if(pt_trg<laeff[2])
                    {
                        effla = efflal_2(pt_trg);
                    }
                    else
                    {
                        effla = efflah_2(pt_trg);
                    }
                }
                if(eta_trg>etacut[3] && eta_trg<etacut[4])
                {
                    if(pt_trg<laeff[3])
                    {
                        effla = efflal_3(pt_trg);
                    }
                    else
                    {
                        effla = efflah_3(pt_trg);
                    }
                }
                if(eta_trg>etacut[4] && eta_trg<etacut[5])
                {
                    if(pt_trg<laeff[4])
                    {
                        effla = efflal_4(pt_trg);
                    }
                    else
                    {
                        effla = efflah_4(pt_trg);
                    }
                }
                if(eta_trg>etacut[5] && eta_trg<etacut[6])
                {
                    if(pt_trg<laeff[5])
                    {
                        effla = efflal_5(pt_trg);
                    }
                    else
                    {
                        effla = efflah_5(pt_trg);
                    }
                }*/
                
                
                TVector3 pvector_trg_dau1 = (*pVect_dau_la[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();
                
                TVector3 pvector_trg_dau2 = (*pVect_dau_la[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();
                
                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    
                    if(eta_ass==eta_trg_dau1 && phi_ass==phi_trg_dau1) continue;
                    if(eta_ass==eta_trg_dau2 && phi_ass==phi_trg_dau2) continue;
                    
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;
                    
                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_la[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_la);
                }
            }

        }

        for(int i=0; i<13; i++)
        {
            int nMult_trg_ks = (int)pVect_trg_ks_bkg[i]->size();
            int nMult_trg_la = (int)pVect_trg_la_bkg[i]->size();
            hMult_ks_bkg[i]->Fill(nMult_trg_ks);
            hMult_la_bkg[i]->Fill(nMult_trg_la);
            for(int ntrg=0;ntrg<nMult_trg_ks;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg_ks_bkg[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double effks = 1.0;
                if(pt_trg<1.2)
                {
                    effks = effkss(pt_trg);
                }
                else
                {
                    effks = effksb(pt_trg);
                }
                
                
                /*if(eta_trg>etacut[0] && eta_trg<etacut[1])
                 {
                 if(pt_trg<kseff[0])
                 {
                 effks = effksl_0(pt_trg);
                 }
                 else
                 {
                 effks = effksh_0(pt_trg);
                 }
                 }
                 if(eta_trg>etacut[1] && eta_trg<etacut[2])
                 {
                 if(pt_trg<kseff[1])
                 {
                 effks = effksl_1(pt_trg);
                 }
                 else
                 {
                 effks = effksh_1(pt_trg);
                 }
                 }
                 if(eta_trg>etacut[2] && eta_trg<etacut[3])
                 {
                 if(pt_trg<kseff[2])
                 {
                 effks = effksl_2(pt_trg);
                 }
                 else
                 {
                 effks = effksh_2(pt_trg);
                 }
                 }
                 if(eta_trg>etacut[3] && eta_trg<etacut[4])
                 {
                 if(pt_trg<kseff[3])
                 {
                 effks = effksl_3(pt_trg);
                 }
                 else
                 {
                 effks = effksh_3(pt_trg);
                 }
                 }
                 if(eta_trg>etacut[4] && eta_trg<etacut[5])
                 {
                 if(pt_trg<kseff[4])
                 {
                 effks = effksl_4(pt_trg);
                 }
                 else
                 {
                 effks = effksh_4(pt_trg);
                 }
                 }
                 if(eta_trg>etacut[5] && eta_trg<etacut[6])
                 {
                 if(pt_trg<kseff[5])
                 {
                 effks = effksl_5(pt_trg);
                 }
                 else
                 {
                 effks = effksh_5(pt_trg);
                 }
                 }*/
                
                
                TVector3 pvector_trg_dau1 = (*pVect_dau_ks_bkg[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();
                
                TVector3 pvector_trg_dau2 = (*pVect_dau_ks_bkg[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();
                
                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    
                    if(eta_ass==eta_trg_dau1 && phi_ass==phi_trg_dau1) continue;
                    if(eta_ass==eta_trg_dau2 && phi_ass==phi_trg_dau2) continue;
                    
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;
                    
                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_ks_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_ks);
                }
            }
            
            for(int ntrg=0;ntrg<nMult_trg_la;ntrg++)
            {
                TVector3 pvector_trg = (*pVect_trg_la_bkg[i])[ntrg];
                double eta_trg = pvector_trg.Eta();
                double phi_trg = pvector_trg.Phi();
                double pt_trg = pvector_trg.Pt();
                double effla = 1.0;
                if(pt_trg<1.6)
                {
                    effla = efflas(pt_trg);
                }
                else
                {
                    effla = efflab(pt_trg);
                }
                
                /*if(eta_trg>etacut[0] && eta_trg<etacut[1])
                 {
                 if(pt_trg<laeff[0])
                 {
                 effla = efflal_0(pt_trg);
                 }
                 else
                 {
                 effla = efflah_0(pt_trg);
                 }
                 }
                 if(eta_trg>etacut[1] && eta_trg<etacut[2])
                 {
                 if(pt_trg<laeff[1])
                 {
                 effla = efflal_1(pt_trg);
                 }
                 else
                 {
                 effla = efflah_1(pt_trg);
                 }
                 }
                 if(eta_trg>etacut[2] && eta_trg<etacut[3])
                 {
                 if(pt_trg<laeff[2])
                 {
                 effla = efflal_2(pt_trg);
                 }
                 else
                 {
                 effla = efflah_2(pt_trg);
                 }
                 }
                 if(eta_trg>etacut[3] && eta_trg<etacut[4])
                 {
                 if(pt_trg<laeff[3])
                 {
                 effla = efflal_3(pt_trg);
                 }
                 else
                 {
                 effla = efflah_3(pt_trg);
                 }
                 }
                 if(eta_trg>etacut[4] && eta_trg<etacut[5])
                 {
                 if(pt_trg<laeff[4])
                 {
                 effla = efflal_4(pt_trg);
                 }
                 else
                 {
                 effla = efflah_4(pt_trg);
                 }
                 }
                 if(eta_trg>etacut[5] && eta_trg<etacut[6])
                 {
                 if(pt_trg<laeff[5])
                 {
                 effla = efflal_5(pt_trg);
                 }
                 else
                 {
                 effla = efflah_5(pt_trg);
                 }
                 }*/
                
                
                TVector3 pvector_trg_dau1 = (*pVect_dau_la_bkg[i])[2*ntrg];
                double eta_trg_dau1 = pvector_trg_dau1.Eta();
                double phi_trg_dau1 = pvector_trg_dau1.Phi();
                
                TVector3 pvector_trg_dau2 = (*pVect_dau_la_bkg[i])[2*ntrg+1];
                double eta_trg_dau2 = pvector_trg_dau2.Eta();
                double phi_trg_dau2 = pvector_trg_dau2.Phi();
                
                for(int nass=0;nass<nMult_ass;nass++)
                {
                    TVector3 pvector_ass = (*pVect_ass)[nass];
                    double eta_ass = pvector_ass.Eta();
                    double phi_ass = pvector_ass.Phi();
                    
                    if(eta_ass==eta_trg_dau1 && phi_ass==phi_trg_dau1) continue;
                    if(eta_ass==eta_trg_dau2 && phi_ass==phi_trg_dau2) continue;
                    
                    double deltaEta=eta_ass-eta_trg;
                    double deltaPhi=phi_ass-phi_trg;
                    if(deltaPhi>PI)
                        deltaPhi=deltaPhi-2*PI;
                    if(deltaPhi<-PI)
                        deltaPhi=deltaPhi+2*PI;
                    if(deltaPhi>-PI && deltaPhi<-PI/2.)
                        deltaPhi=deltaPhi+2*PI;
                    
                    //if(deltaEta==0 && deltaPhi==0) continue;
                    hSignal_la_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg_la);
                }
            }
            
        }
        
        
        for(int i=0; i<13; i++)
        {
            pVectVect_trg_ks[i]->push_back(*pVect_trg_ks[i]);
            pVectVect_trg_la[i]->push_back(*pVect_trg_la[i]);
            delete pVect_trg_ks[i];
            delete pVect_trg_la[i];
            delete pVect_dau_ks[i];
            delete pVect_dau_la[i];
            pVectVect_trg_ks_bkg[i]->push_back(*pVect_trg_ks_bkg[i]);
            pVectVect_trg_la_bkg[i]->push_back(*pVect_trg_la_bkg[i]);
            delete pVect_trg_ks_bkg[i];
            delete pVect_trg_la_bkg[i];
            delete pVect_dau_ks_bkg[i];
            delete pVect_dau_la_bkg[i];
        }
        
        pVectVect_ass->push_back(*pVect_ass);
        zvtxVect->push_back(bestvz);
        delete pVect_ass;
    }
}


// ------------ method called once each job just before starting event
//loop  ------------
void
DemoAnalyzerAll220::beginJob()
{
    edm::Service<TFileService> fs;
    
    TH1D::SetDefaultSumw2();
    
    
    hMult = fs->make<TH1D>("mult",";N",300,0,300);
    hMult_ass = fs->make<TH1D>("mult_ass",";N",600,0,600);
    
    for(int i=0; i<13; i++)
    {
        hKET_ks[i] = fs->make<TH1D>(Form("KETkshort_pt%d",i),";GeV",25000,0,12.5);
        hKET_la[i] = fs->make<TH1D>(Form("KETlambda_pt%d",i),";GeV",25000,0,12.5);
        hPt_ks[i] = fs->make<TH1D>(Form("Ptkshort_pt%d",i),";GeV",25000,0,12.5);
        hPt_la[i] = fs->make<TH1D>(Form("Ptlambda_pt%d",i),";GeV",25000,0,12.5);
        hMass_ks[i] = fs->make<TH1D>(Form("masskshort_pt%d",i),";GeV",2000,0,1.0);
        hMass_la[i] = fs->make<TH1D>(Form("masslambda_pt%d",i),";GeV",2000,0.5,1.5);
        hMult_ks[i] = fs->make<TH1D>(Form("mult_ks_pt%d",i),";N",250,0,250);
        hSignal_ks[i] = fs->make<TH2D>(Form("signalkshort_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_ks[i] = fs->make<TH2D>(Form("backgroundkshort_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hMult_la[i] = fs->make<TH1D>(Form("mult_la_pt%d",i),";N",250,0,250);
        hSignal_la[i] = fs->make<TH2D>(Form("signallambda_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_la[i] = fs->make<TH2D>(Form("backgroundlambda_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        pVectVect_trg_ks[i] = new vector< vector<TVector3> >;
        pVectVect_trg_la[i] = new vector< vector<TVector3> >;
        hKET_ks_bkg[i] = fs->make<TH1D>(Form("KETkshort_bkg_pt%d",i),";GeV",25000,0,12.5);
        hKET_la_bkg[i] = fs->make<TH1D>(Form("KETlambda_bkg_pt%d",i),";GeV",25000,0,12.5);
        hPt_ks_bkg[i] = fs->make<TH1D>(Form("Ptkshort_bkg_pt%d",i),";GeV",25000,0,12.5);
        hPt_la_bkg[i] = fs->make<TH1D>(Form("Ptlambda_bkg_pt%d",i),";GeV",25000,0,12.5);
        hMult_ks_bkg[i] = fs->make<TH1D>(Form("mult_ks_bkg_pt%d",i),";N",250,0,250);
        hSignal_ks_bkg[i] = fs->make<TH2D>(Form("signalkshort_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_ks_bkg[i] = fs->make<TH2D>(Form("backgroundkshort_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hMult_la_bkg[i] = fs->make<TH1D>(Form("mult_la_bkg_pt%d",i),";N",250,0,250);
        hSignal_la_bkg[i] = fs->make<TH2D>(Form("signallambda_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        hBackground_la_bkg[i] = fs->make<TH2D>(Form("backgroundlambda_bkg_pt%d",i),";#Delta#eta;#Delta#phi",33,-4.95,4.95,31,-(0.5-1.0/32)*PI,(1.5-1.0/32)*PI);
        pVectVect_trg_ks_bkg[i] = new vector< vector<TVector3> >;
        pVectVect_trg_la_bkg[i] = new vector< vector<TVector3> >;
    }
    
    pVectVect_ass = new vector< vector<TVector3> >;
    zvtxVect = new vector<double>;
    
}

// ------------ method called once each job just after ending the event
//loop  ------------
void
DemoAnalyzerAll220::endJob() {
    //Calculating background
    int nevttotal_ass = (int)pVectVect_ass->size();
    
    /*double etacut[7] = {-2.4,-1.6,-0.8,0,0.8,1.6,2.4};
    double kseff[6] = {1.2,1.2,1.2,1.2,1.2,1.2};
    double laeff[6] = {2,2,2,2,2,1.8};*/

    for(int i=0;i<13;i++)
    {
        int nevttotal_trg_ks = (int)pVectVect_trg_ks[i]->size();
        int nevttotal_trg_la = (int)pVectVect_trg_la[i]->size();
        
        for(int nround=0;nround<bkgFactor_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_ks; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }
                
                vector<TVector3> pVectTmp_trg = (*pVectVect_trg_ks[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double effks = 1.0;
                    
                    if(pt_trg<1.2)
                    {
                        effks = effkss(pt_trg);
                    }
                    else
                    {
                        effks = effksb(pt_trg);
                    }

                    
                    /*if(eta_trg>etacut[0] && eta_trg<etacut[1])
                    {
                        if(pt_trg<kseff[0])
                        {
                            effks = effksl_0(pt_trg);
                        }
                        else
                        {
                            effks = effksh_0(pt_trg);
                        }
                    }
                    if(eta_trg>etacut[1] && eta_trg<etacut[2])
                    {
                        if(pt_trg<kseff[1])
                        {
                            effks = effksl_1(pt_trg);
                        }
                        else
                        {
                            effks = effksh_1(pt_trg);
                        }
                    }
                    if(eta_trg>etacut[2] && eta_trg<etacut[3])
                    {
                        if(pt_trg<kseff[2])
                        {
                            effks = effksl_2(pt_trg);
                        }
                        else
                        {
                            effks = effksh_2(pt_trg);
                        }
                    }
                    if(eta_trg>etacut[3] && eta_trg<etacut[4])
                    {
                        if(pt_trg<kseff[3])
                        {
                            effks = effksl_3(pt_trg);
                        }
                        else
                        {
                            effks = effksh_3(pt_trg);
                        }
                    }
                    if(eta_trg>etacut[4] && eta_trg<etacut[5])
                    {
                        if(pt_trg<kseff[4])
                        {
                            effks = effksl_4(pt_trg);
                        }
                        else
                        {
                            effks = effksh_4(pt_trg);
                        }
                    }
                    if(eta_trg>etacut[5] && eta_trg<etacut[6])
                    {
                        if(pt_trg<kseff[5])
                        {
                            effks = effksl_5(pt_trg);
                        }
                        else
                        {
                            effks = effksh_5(pt_trg);
                        }
                    }*/
                    
                    
                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        
                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;
                        
                        if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_ks[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg);
                    }
                }
            }
        }
        
        for(int nround=0;nround<bkgFactor_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_la; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }
                
                vector<TVector3> pVectTmp_trg = (*pVectVect_trg_la[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double effla = 1.0;
                    
                    if(pt_trg<1.6)
                    {
                        effla = efflas(pt_trg);
                    }
                    else
                    {
                        effla = efflab(pt_trg);
                    }

                    
                    /*if(eta_trg>etacut[0] && eta_trg<etacut[1])
                    {
                        if(pt_trg<laeff[0])
                        {
                            effla = efflal_0(pt_trg);
                        }
                        else
                        {
                            effla = efflah_0(pt_trg);
                        }
                    }
                    if(eta_trg>etacut[1] && eta_trg<etacut[2])
                    {
                        if(pt_trg<laeff[1])
                        {
                            effla = efflal_1(pt_trg);
                        }
                        else
                        {
                            effla = efflah_1(pt_trg);
                        }
                    }
                    if(eta_trg>etacut[2] && eta_trg<etacut[3])
                    {
                        if(pt_trg<laeff[2])
                        {
                            effla = efflal_2(pt_trg);
                        }
                        else
                        {
                            effla = efflah_2(pt_trg);
                        }
                    }
                    if(eta_trg>etacut[3] && eta_trg<etacut[4])
                    {
                        if(pt_trg<laeff[3])
                        {
                            effla = efflal_3(pt_trg);
                        }
                        else
                        {
                            effla = efflah_3(pt_trg);
                        }
                    }
                    if(eta_trg>etacut[4] && eta_trg<etacut[5])
                    {
                        if(pt_trg<laeff[4])
                        {
                            effla = efflal_4(pt_trg);
                        }
                        else
                        {
                            effla = efflah_4(pt_trg);
                        }
                    }
                    if(eta_trg>etacut[5] && eta_trg<etacut[6])
                    {
                        if(pt_trg<laeff[5])
                        {
                            effla = efflal_5(pt_trg);
                        }
                        else
                        {
                            effla = efflah_5(pt_trg);
                        }
                    }*/
                    
                    
                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        
                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;
                        
                        if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_la[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg);
                    }
                }
            }
        }

    }
    
    for(int i=0;i<13;i++)
    {
        int nevttotal_trg_ks = (int)pVectVect_trg_ks_bkg[i]->size();
        int nevttotal_trg_la = (int)pVectVect_trg_la_bkg[i]->size();
        
        for(int nround=0;nround<bkgFactor_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_ks; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }
                
                vector<TVector3> pVectTmp_trg = (*pVectVect_trg_ks_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double effks = 1.0;
                    if(pt_trg<1.2)
                    {
                        effks = effkss(pt_trg);
                    }
                    else
                    {
                        effks = effksb(pt_trg);
                    }
                    
                    
                    /*if(eta_trg>etacut[0] && eta_trg<etacut[1])
                     {
                     if(pt_trg<kseff[0])
                     {
                     effks = effksl_0(pt_trg);
                     }
                     else
                     {
                     effks = effksh_0(pt_trg);
                     }
                     }
                     if(eta_trg>etacut[1] && eta_trg<etacut[2])
                     {
                     if(pt_trg<kseff[1])
                     {
                     effks = effksl_1(pt_trg);
                     }
                     else
                     {
                     effks = effksh_1(pt_trg);
                     }
                     }
                     if(eta_trg>etacut[2] && eta_trg<etacut[3])
                     {
                     if(pt_trg<kseff[2])
                     {
                     effks = effksl_2(pt_trg);
                     }
                     else
                     {
                     effks = effksh_2(pt_trg);
                     }
                     }
                     if(eta_trg>etacut[3] && eta_trg<etacut[4])
                     {
                     if(pt_trg<kseff[3])
                     {
                     effks = effksl_3(pt_trg);
                     }
                     else
                     {
                     effks = effksh_3(pt_trg);
                     }
                     }
                     if(eta_trg>etacut[4] && eta_trg<etacut[5])
                     {
                     if(pt_trg<kseff[4])
                     {
                     effks = effksl_4(pt_trg);
                     }
                     else
                     {
                     effks = effksh_4(pt_trg);
                     }
                     }
                     if(eta_trg>etacut[5] && eta_trg<etacut[6])
                     {
                     if(pt_trg<kseff[5])
                     {
                     effks = effksl_5(pt_trg);
                     }
                     else
                     {
                     effks = effksh_5(pt_trg);
                     }
                     }*/
                    
                    
                    
                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        
                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;
                        
                        if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_ks_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg);
                    }
                }
            }
        }
        
        for(int nround=0;nround<bkgFactor_;nround++)
        {
            int ncount = 0;
            for(int nevt_trg=0; nevt_trg<nevttotal_trg_la; nevt_trg++)
            {
                int nevt_ass = gRandom->Integer(nevttotal_ass);
                if(nevt_trg == nevt_ass) { nevt_trg--; continue; }
                if(fabs((*zvtxVect)[nevt_trg]-(*zvtxVect)[nevt_ass])>0.5) {
                    nevt_trg--;
                    ncount++;
                    if(ncount>5000) {nevt_trg++; ncount = 0;}
                    continue; }
                
                vector<TVector3> pVectTmp_trg = (*pVectVect_trg_la_bkg[i])[nevt_trg];
                vector<TVector3> pVectTmp_ass = (*pVectVect_ass)[nevt_ass];
                int nMult_trg = pVectTmp_trg.size();
                int nMult_ass = pVectTmp_ass.size();
                
                for(int ntrg=0;ntrg<nMult_trg;ntrg++)
                {
                    TVector3 pvectorTmp_trg = pVectTmp_trg[ntrg];
                    double eta_trg = pvectorTmp_trg.Eta();
                    double phi_trg = pvectorTmp_trg.Phi();
                    double pt_trg = pvectorTmp_trg.Pt();
                    double effla = 1.0;
                    if(pt_trg<1.6)
                    {
                        effla = efflas(pt_trg);
                    }
                    else
                    {
                        effla = efflab(pt_trg);
                    }
                    
                    
                    /*if(eta_trg>etacut[0] && eta_trg<etacut[1])
                     {
                     if(pt_trg<laeff[0])
                     {
                     effla = efflal_0(pt_trg);
                     }
                     else
                     {
                     effla = efflah_0(pt_trg);
                     }
                     }
                     if(eta_trg>etacut[1] && eta_trg<etacut[2])
                     {
                     if(pt_trg<laeff[1])
                     {
                     effla = efflal_1(pt_trg);
                     }
                     else
                     {
                     effla = efflah_1(pt_trg);
                     }
                     }
                     if(eta_trg>etacut[2] && eta_trg<etacut[3])
                     {
                     if(pt_trg<laeff[2])
                     {
                     effla = efflal_2(pt_trg);
                     }
                     else
                     {
                     effla = efflah_2(pt_trg);
                     }
                     }
                     if(eta_trg>etacut[3] && eta_trg<etacut[4])
                     {
                     if(pt_trg<laeff[3])
                     {
                     effla = efflal_3(pt_trg);
                     }
                     else
                     {
                     effla = efflah_3(pt_trg);
                     }
                     }
                     if(eta_trg>etacut[4] && eta_trg<etacut[5])
                     {
                     if(pt_trg<laeff[4])
                     {
                     effla = efflal_4(pt_trg);
                     }
                     else
                     {
                     effla = efflah_4(pt_trg);
                     }
                     }
                     if(eta_trg>etacut[5] && eta_trg<etacut[6])
                     {
                     if(pt_trg<laeff[5])
                     {
                     effla = efflal_5(pt_trg);
                     }
                     else
                     {
                     effla = efflah_5(pt_trg);
                     }
                     }*/
                    
                    
                    for(int nass=0;nass<nMult_ass;nass++)
                    {
                        TVector3 pvectorTmp_ass = pVectTmp_ass[nass];
                        double eta_ass = pvectorTmp_ass.Eta();
                        double phi_ass = pvectorTmp_ass.Phi();
                        
                        double deltaEta=eta_ass-eta_trg;
                        double deltaPhi=phi_ass-phi_trg;
                        if(deltaPhi>PI)
                            deltaPhi=deltaPhi-2*PI;
                        if(deltaPhi<-PI)
                            deltaPhi=deltaPhi+2*PI;
                        if(deltaPhi>-PI && deltaPhi<-PI/2.)
                            deltaPhi=deltaPhi+2*PI;
                        
                        if(deltaEta==0 && deltaPhi==0) continue;
                        hBackground_la_bkg[i]->Fill(deltaEta,deltaPhi,1.0/nMult_trg);
                    }
                }
            }
        }
        
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzerAll220);


