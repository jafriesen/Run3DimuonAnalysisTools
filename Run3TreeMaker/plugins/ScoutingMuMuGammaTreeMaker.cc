// -*- C++ -*-
//
// Package:    Run3DimuonAnalysisTools/ScoutingMuMuGammaTreeMaker
// Class:      ScoutingMuMuGammaTreeMaker
//
/**\class ScoutingMuMuGammaTreeMaker ScoutingMuMuGammaTreeMaker.cc Run3DimuonAnalysisTools/ScoutingMuMuGammaTreeMaker/plugins/ScoutingMuMuGammaTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/

// system include files
#include <memory>
#include <TTree.h>
#include <TLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class ScoutingMuMuGammaTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingMuMuGammaTreeMaker(const edm::ParameterSet&);
  ~ScoutingMuMuGammaTreeMaker() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
    //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;   

  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>             triggerResultsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon> >      muonsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> >  electronsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >    primaryVerticesToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >    verticesToken;
  const edm::EDGetTokenT<double>                          rhoToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> >  photonsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle> >  pfcandsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPFJet> >  pfjetsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingTrack> >  tracksToken;

  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;

  bool doL1;
  triggerExpression::Data triggerCache_;

  edm::InputTag                algInputTag_;
  edm::InputTag                extInputTag_;
  edm::EDGetToken              algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;
  std::vector<bool>            l1Result_;

  TTree* tree;

  float trackIso1;
  float trackIso2;
  int nValidPixelHits1;
  int nValidPixelHits2;
  int nTrackerLayersWithMeasurement1;
  int nTrackerLayersWithMeasurement2;
  float trk_chi21;
  float trk_chi22;

  float mass;
  float pt;
  float dr;
  float pt1;
  float pt2;

  float eta1;
  float eta2;
  float phi1;
  float phi2;

  float rho;

  int npvtx;

  float vtxChi2;

  float avgPrimaryX;
  float avgPrimaryY;
  float avgPrimaryZ;

  float vtxX;
  float vtxY;
  float vtxZ; 

  float vtxXError;
  float vtxYError;
  float vtxZError;

  float Lxy;
  float LxyErr;
  float LxySig;

  int nPhotons;

  std::vector<float> photonPt;
  std::vector<float> photonEta;
  std::vector<float> photonPhi;
  std::vector<float> photonM;
  std::vector<float> photonSigmaIetaIeta;
  std::vector<float> photonHOverE;
  std::vector<float> photonEcalIso;
  std::vector<float> photonHcalIso;
  std::vector<float> photonTrkIso;
  std::vector<float> photonR9;
  std::vector<float> photonSMin;
  std::vector<float> photonSMaj;
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
ScoutingMuMuGammaTreeMaker::ScoutingMuMuGammaTreeMaker(const edm::ParameterSet& iConfig):
    triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
    triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),
    muonsToken               (consumes<std::vector<Run3ScoutingMuon> >             (iConfig.getParameter<edm::InputTag>("muons"))),
    electronsToken           (consumes<std::vector<Run3ScoutingElectron> >         (iConfig.getParameter<edm::InputTag>("electrons"))),
    primaryVerticesToken     (consumes<std::vector<Run3ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("primaryVertices"))),
    verticesToken            (consumes<std::vector<Run3ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("displacedVertices"))),
    rhoToken                 (consumes<double>                                 (iConfig.getParameter<edm::InputTag>("rho"))), 
    photonsToken             (consumes<std::vector<Run3ScoutingPhoton> >         (iConfig.getParameter<edm::InputTag>("photons"))),
    pfcandsToken             (consumes<std::vector<Run3ScoutingParticle> >         (iConfig.getParameter<edm::InputTag>("pfcands"))),
    pfjetsToken              (consumes<std::vector<Run3ScoutingPFJet> >            (iConfig.getParameter<edm::InputTag>("pfjets"))),
    tracksToken              (consumes<std::vector<Run3ScoutingTrack> >            (iConfig.getParameter<edm::InputTag>("tracks"))),
    doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1")            : false)
{
    usesResource("TFileService");
    if (doL1) {
        algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
        extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
        algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
        l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
        l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
    }
    else {
        l1Seeds_ = std::vector<std::string>();
        l1GtUtils_ = 0;
    }
}

ScoutingMuMuGammaTreeMaker::~ScoutingMuMuGammaTreeMaker() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void ScoutingMuMuGammaTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;

  Handle<vector<Run3ScoutingVertex> > verticesH;
  iEvent.getByToken(verticesToken, verticesH);
  int ndvtx = verticesH->size();

  if ( ndvtx == 0 ) return;

  Handle<vector<Run3ScoutingMuon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);

  if ( muonsH->size() < 2 ) return;

  bool fillTree = false;
  int vtxIndx = 0;
  int idx[2];

  //std::cout << std::endl;

  // iterate through displaced vertex indices
  for ( unsigned int iVtx = 0; iVtx < verticesH->size(); iVtx++ ) {

    //std::cout<<std::endl<<"vertex "<<iVtx<< " vtxChi2: " << (verticesH->at(iVtx)).chi2()<<std::endl;

    int charge = 1;
    std::vector<int> muIdx;

    // iterate through muon collection
    for ( unsigned int iMuon = 0; iMuon < muonsH->size(); iMuon++ ) {
      //std::cout<<"    muon "<<iMuon<<" pt: "<<(muonsH->at(iMuon)).pt()<<std::endl;
      if ( (muonsH->at(iMuon)).pt() > 3 ) {
        // iterate through muon vertex indices
        for ( auto& muVtxIndx : (muonsH->at(iMuon)).vtxIndx() ) {
          //std::cout<<"        muVtxIndx "<<muVtxIndx<<std::endl;

          // if muon pt > 3 and muon has the current vertex index, save the muon index and charge for later
          if ( (unsigned)muVtxIndx == iVtx ) {
            muIdx.push_back(iMuon);
            charge *= (muonsH->at(iMuon)).charge();

            //std::cout<<"            iMuon "<<iMuon<<" charge "<<(muonsH->at(iMuon)).charge()<<std::endl;
          }
        }
      }
    }

    //std::cout<<"muIdx.size() "<<muIdx.size()<<" charge "<<charge<<std::endl;
    if (muIdx.size() == 2 && charge < 0) {
      if ( fillTree == false || vtxChi2 > (verticesH->at(iVtx)).chi2() ) {
        vtxChi2 = (verticesH->at(iVtx)).chi2();
        vtxIndx = iVtx;
        idx[0] = muIdx[0];
        idx[1] = muIdx[1];
        fillTree = true;

        //std::cout << "saving: " << idx[0] << idx[1] << std::endl;
      }
    }
  }

  if (fillTree) {

    trackIso1 = muonsH->at(idx[0]).trackIso();
    trackIso2 = muonsH->at(idx[1]).trackIso();   
    nValidPixelHits1 = muonsH->at(idx[0]).nValidPixelHits();
    nValidPixelHits2 = muonsH->at(idx[1]).nValidPixelHits();
    nTrackerLayersWithMeasurement1 = muonsH->at(idx[0]).nTrackerLayersWithMeasurement();
    nTrackerLayersWithMeasurement2 = muonsH->at(idx[1]).nTrackerLayersWithMeasurement();
    trk_chi21 = muonsH->at(idx[0]).trk_chi2();
    trk_chi22 = muonsH->at(idx[1]).trk_chi2();

    pt1=muonsH->at(idx[0]).pt();
    pt2=muonsH->at(idx[1]).pt();

    eta1=muonsH->at(idx[0]).eta();
    eta2=muonsH->at(idx[1]).eta();
    phi1=muonsH->at(idx[0]).phi();
    phi2=muonsH->at(idx[1]).phi();
    
    TLorentzVector mu1;
    mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.105658);

    TLorentzVector mu2;
    mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.105658);

    TLorentzVector dimu = mu1+mu2;
    mass=dimu.M();
    pt=dimu.Pt();
    dr=mu1.DeltaR(mu2);

    //std::cout<<"muon1 pt: "<<pt1<<" eta: "<<eta1<<" phi: "<<phi1<<std::endl;
    //std::cout<<"muon2 pt: "<<pt2<<" eta: "<<eta2<<" phi: "<<phi2<<std::endl;

    Handle<double> rhoH;
    iEvent.getByToken(rhoToken, rhoH);
    rho=*rhoH;

    Handle<vector<Run3ScoutingVertex> > primaryVerticesH;
    iEvent.getByToken(primaryVerticesToken, primaryVerticesH);

    std::vector<float> pVtxX;
    std::vector<float> pVtxY;
    std::vector<float> pVtxZ;

    npvtx = 0;
    for (auto vtx_iter = primaryVerticesH->begin(); vtx_iter != primaryVerticesH->end(); ++vtx_iter) {
      //std::cout<<"primary x: "<<vtx_iter->x() <<" y: "<<  vtx_iter->y()<<" z: "<<  vtx_iter->z()<<" ex: "<<vtx_iter->xError()<<" ey: "<<vtx_iter->yError()<<" ez: "<<vtx_iter->zError()<<std::endl;
      npvtx++;
      pVtxX.push_back(vtx_iter->x());
      pVtxY.push_back(vtx_iter->y());
      pVtxZ.push_back(vtx_iter->z());
    }

    avgPrimaryX = ( pVtxX.empty() ) ? 0 : ( std::reduce(pVtxX.begin(), pVtxX.end(), 0.0) / pVtxX.size() );
    avgPrimaryY = ( pVtxY.empty() ) ? 0 : ( std::reduce(pVtxY.begin(), pVtxY.end(), 0.0) / pVtxY.size() );
    avgPrimaryZ = ( pVtxZ.empty() ) ? 0 : ( std::reduce(pVtxZ.begin(), pVtxZ.end(), 0.0) / pVtxZ.size() );

    //std::cout << "npvtx: " << npvtx << " avgPrimaryX: " << avgPrimaryX << " avgPrimaryY: " << avgPrimaryY << " avgPrimaryZ: " << avgPrimaryZ << std::endl;

    vtxX = ((verticesH->at(vtxIndx)).x());
    vtxY = ((verticesH->at(vtxIndx)).y());
    vtxZ = ((verticesH->at(vtxIndx)).z());

    double dx = ((verticesH->at(vtxIndx)).x()) - avgPrimaryX;
    double dy = ((verticesH->at(vtxIndx)).y()) - avgPrimaryY;

    vtxXError =(verticesH->at(vtxIndx)).xError();
    vtxYError =(verticesH->at(vtxIndx)).yError();
    vtxZError =(verticesH->at(vtxIndx)).zError();

    Lxy = sqrt(dx*dx + dy*dy);
    LxyErr = sqrt(dx*dx*vtxXError*vtxXError + dy*dy*vtxYError*vtxYError) / Lxy;
    LxySig = Lxy/LxyErr;

    //std::cout<<"Lxy: "<<Lxy<<" LxyErr: "<<LxyErr<<" LxySig: "<<LxySig<<std::endl;

    Handle<vector<Run3ScoutingPhoton> > photonsH;
    iEvent.getByToken(photonsToken, photonsH);

    photonPt.clear();
    photonEta.clear();
    photonPhi.clear();
    photonM.clear();
    photonSigmaIetaIeta.clear();
    photonHOverE.clear();
    photonEcalIso.clear();
    photonHcalIso.clear();
    photonTrkIso.clear();
    photonR9.clear();
    photonSMin.clear();
    photonSMaj.clear();

    nPhotons = 0;
    for (auto photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
      //std::cout << "photonPt " <<  photons_iter->pt() << std::endl;
      if (photons_iter->pt() > 3) {
        photonPt.push_back(photons_iter->pt());
        photonEta.push_back(photons_iter->eta());
        photonPhi.push_back(photons_iter->phi());
        photonM.push_back(photons_iter->m());
        photonSigmaIetaIeta.push_back(photons_iter->sigmaIetaIeta());
        photonHOverE.push_back(photons_iter->hOverE());
        photonEcalIso.push_back(photons_iter->ecalIso());
        photonHcalIso.push_back(photons_iter->hcalIso());
        photonTrkIso.push_back(photons_iter->trkIso());
        photonR9.push_back(photons_iter->r9());
        photonSMin.push_back(photons_iter->sMin());
        photonSMaj.push_back(photons_iter->sMaj());
        nPhotons++;
      }
    }

    Handle<vector<Run3ScoutingElectron> > electronsH;
    iEvent.getByToken(electronsToken, electronsH);

    for (auto electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) {
      //std::cout << "photonPt " <<  photons_iter->pt() << std::endl;
      if (electrons_iter->pt() > 3) {
        photonPt.push_back(electrons_iter->pt());
        photonEta.push_back(electrons_iter->eta());
        photonPhi.push_back(electrons_iter->phi());
        photonM.push_back(electrons_iter->m());
        photonSigmaIetaIeta.push_back(electrons_iter->sigmaIetaIeta());
        photonHOverE.push_back(electrons_iter->hOverE());
        photonEcalIso.push_back(electrons_iter->ecalIso());
        photonHcalIso.push_back(electrons_iter->hcalIso());
        photonTrkIso.push_back(electrons_iter->trackIso());
        photonR9.push_back(electrons_iter->r9());
        photonSMin.push_back(electrons_iter->sMin());
        photonSMaj.push_back(electrons_iter->sMaj());
      }
    }


    //std::cout << photonPt.size() << " photons" << std::endl;

    l1Result_.clear();
    if (doL1) {
        l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
        for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
            bool l1htbit = 0;
            l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
            l1Result_.push_back( l1htbit );
        }
    }


    int eventNum = iEvent.id().event();
    int lumiSec = iEvent.luminosityBlock();
    int runNum = iEvent.id().run();
    std::cout<<runNum<<":"<<lumiSec<<":"<<eventNum<<std::endl;
    std::cout<<"  Lxy: "<<Lxy<<" LxyErr: "<<LxyErr<<" LxySig: "<<LxySig<<std::endl;
    std::cout<<"  muon1 pt: "<<pt1<<" eta: "<<eta1<<" phi: "<<phi1<<" nValidPixelHits: "<<nValidPixelHits1<<std::endl;
    std::cout<<"  muon2 pt: "<<pt2<<" eta: "<<eta2<<" phi: "<<phi2<<" nValidPixelHits: "<<nValidPixelHits2<<std::endl;
    std::cout<<"  dimu mass: "<<mass<<" dR: "<<dr<<" dPhi: "<<mu1.DeltaPhi(mu2)<<std::endl;
    TLorentzVector vtx4;
    vtx4.SetXYZT(vtxX-avgPrimaryX,vtxY-avgPrimaryY,vtxZ-avgPrimaryZ,0);
    std::cout<<"  vtx x: "<<vtxX<<" y: "<<vtxY<<" vtxZ: "<<vtxZ<<" phi: "<<vtx4.Phi()<<" eta: "<<vtx4.Eta()<<std::endl;
    std::cout<<"  dimu vtx dPhi: "<<vtx4.DeltaPhi(dimu)<<" dR: "<<vtx4.DeltaR(dimu)<<std::endl;
    std::cout<<"  ndvtx: "<<ndvtx<<" muon 1 ndvtx: "<<(muonsH->at(idx[0]).vtxIndx()).size()<<" muon 2 ndvtx: "<<(muonsH->at(idx[1]).vtxIndx()).size()<<std::endl;


    //std::cout<<"tree filling with mass: "<<mass<<", pt: "<<pt<<std::endl;
    tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void ScoutingMuMuGammaTreeMaker::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"      , "tree");

    tree->Branch("trackIso1", &trackIso1, "trackIso1/F");
    tree->Branch("trackIso2", &trackIso2, "trackIso2/F");
    tree->Branch("nValidPixelHits1", &nValidPixelHits1, "nValidPixelHits1/I");
    tree->Branch("nValidPixelHits2", &nValidPixelHits2, "nValidPixelHits2/I");
    tree->Branch("nTrackerLayersWithMeasurement1", &nTrackerLayersWithMeasurement1, "nTrackerLayersWithMeasurement1/I");
    tree->Branch("nTrackerLayersWithMeasurement2", &nTrackerLayersWithMeasurement2, "nTrackerLayersWithMeasurement2/I");
    tree->Branch("trk_chi21", &trk_chi21, "trk_chi21/F");
    tree->Branch("trk_chi22", &trk_chi22, "trk_chi22/F");

    tree->Branch("mass"                , &mass                        , "mass/F"    );
    tree->Branch("pt"                  , &pt                          , "pt/F"      );
    tree->Branch("dr"                  , &dr                          , "dr/F"      );
    tree->Branch("pt1"                 , &pt1                         , "pt1/F"     );
    tree->Branch("pt2"                 , &pt2                         , "pt2/F"     );
    tree->Branch("eta1"                , &eta1                        , "eta1/F"    );
    tree->Branch("eta2"                , &eta2                        , "eta2/F"    );
    tree->Branch("phi1"                , &phi1                        , "phi1/F"    );
    tree->Branch("phi2"                , &phi2                        , "phi2/F"    );
    tree->Branch("rho"                 , &rho                         , "rho/F"     );

    tree->Branch("npvtx"               , &npvtx                       , "npvtx/I"   );

    tree->Branch("vtxChi2"             , &vtxChi2                     , "vtxChi2/F" );
    tree->Branch("Lxy"                 , &Lxy                         , "Lxy/F"     );
    tree->Branch("LxyErr"              , &LxyErr                      , "LxyErr/F"  );
    tree->Branch("LxySig"              , &LxySig                      , "LxySig/F"  );

    tree->Branch("avgPrimaryX"           , &avgPrimaryX                   , "avgPrimaryX/F");
    tree->Branch("avgPrimaryY"           , &avgPrimaryY                   , "avgPrimaryY/F");
    tree->Branch("avgPrimaryZ"           , &avgPrimaryZ                   , "avgPrimaryZ/F");

    tree->Branch("vtxX"           , &vtxX                   , "vtxX/F");
    tree->Branch("vtxY"           , &vtxY                   , "vtxY/F");
    tree->Branch("vtxZ"           , &vtxZ                   , "vtxZ/F");

    tree->Branch("vtxXError"           , &vtxXError                   , "vtxXError/F");
    tree->Branch("vtxYError"           , &vtxYError                   , "vtxYError/F");
    tree->Branch("vtxZError"           , &vtxZError                   , "vtxZError/F");

    tree->Branch("l1Result", "std::vector<bool>"             ,&l1Result_, 32000, 0  );

    tree->Branch("nPhotons"               , &nPhotons                       , "nPhotons/I"   );

    tree->Branch("photonPt", "std::vector<float>", &photonPt, 32000, 0);
    tree->Branch("photonEta", "std::vector<float>", &photonEta, 32000, 0);
    tree->Branch("photonPhi", "std::vector<float>", &photonPhi, 32000, 0);
    tree->Branch("photonM", "std::vector<float>", &photonM, 32000, 0);
    tree->Branch("photonSigmaIetaIeta", "std::vector<float>", &photonSigmaIetaIeta, 32000, 0);
    tree->Branch("photonHOverE", "std::vector<float>", &photonHOverE, 32000, 0);
    tree->Branch("photonEcalIso", "std::vector<float>", &photonEcalIso, 32000, 0);
    tree->Branch("photonHcalIso", "std::vector<float>", &photonHcalIso, 32000, 0);
    tree->Branch("photonTrkIso", "std::vector<float>", &photonTrkIso, 32000, 0);
    tree->Branch("photonR9", "std::vector<float>", &photonR9, 32000, 0);
    tree->Branch("photonSMin", "std::vector<float>", &photonSMin, 32000, 0);
    tree->Branch("photonSMaj", "std::vector<float>", &photonSMaj, 32000, 0);
}

// ------------ method called once each job just after ending the event loop  ------------
void ScoutingMuMuGammaTreeMaker::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingMuMuGammaTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingMuMuGammaTreeMaker);
