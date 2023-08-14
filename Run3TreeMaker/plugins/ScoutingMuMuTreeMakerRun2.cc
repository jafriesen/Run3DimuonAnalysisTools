// -*- C++ -*-
//
// Package:    Run3DimuonAnalysisTools/ScoutingMuMuTreeMakerRun2
// Class:      ScoutingMuMuTreeMakerRun2
//
/**\class ScoutingMuMuTreeMakerRun2 ScoutingMuMuTreeMakerRun2.cc Run3DimuonAnalysisTools/ScoutingMuMuTreeMakerRun2/plugins/ScoutingMuMuTreeMakerRun2.cc

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

#include "DataFormats/Scouting/interface/ScoutingElectron.h"
#include "DataFormats/Scouting/interface/ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/ScoutingVertex.h"
#include "DataFormats/Scouting/interface/ScoutingTrack.h"
#include "DataFormats/Scouting/interface/ScoutingMuon.h"
#include "DataFormats/Scouting/interface/ScoutingParticle.h"

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

class ScoutingMuMuTreeMakerRun2 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingMuMuTreeMakerRun2(const edm::ParameterSet&);
  ~ScoutingMuMuTreeMakerRun2() override;

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
  const edm::EDGetTokenT<std::vector<ScoutingMuon> >      muonsToken;
  const edm::EDGetTokenT<std::vector<ScoutingVertex> >    primaryVerticesToken;
  const edm::EDGetTokenT<std::vector<ScoutingVertex> >    verticesToken;
  const edm::EDGetTokenT<double>                          rhoToken;
  const edm::EDGetTokenT<std::vector<ScoutingPFJet> >  pfjetsToken;
  const edm::EDGetTokenT<std::vector<ScoutingTrack> >  tracksToken;

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

  float vtxXError;
  float vtxYError;
  float vtxZError;

  float Lxy;
  float LxyErr;
  float LxySig;

  // beam spot corresponding to 2018D
  const float BSx = 0.0965;
  const float BSy = -0.062;
  const float phi0 = -0.571;
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
ScoutingMuMuTreeMakerRun2::ScoutingMuMuTreeMakerRun2(const edm::ParameterSet& iConfig):
    triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
    triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),
    muonsToken               (consumes<std::vector<ScoutingMuon> >             (iConfig.getParameter<edm::InputTag>("muons"))),
    primaryVerticesToken     (consumes<std::vector<ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("primaryVertices"))),
    verticesToken            (consumes<std::vector<ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("displacedVertices"))),
    rhoToken                 (consumes<double>                                 (iConfig.getParameter<edm::InputTag>("rho"))), 
    pfjetsToken              (consumes<std::vector<ScoutingPFJet> >            (iConfig.getParameter<edm::InputTag>("pfjets"))),
    tracksToken              (consumes<std::vector<ScoutingTrack> >            (iConfig.getParameter<edm::InputTag>("tracks"))),

    doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1")            : false)
{
    usesResource("TFileService");
    if (doL1) {
        algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
        extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
        algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
        l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
        //l1GtUtils_ = new l1t::L1TGlobalUtil(iConfig,consumesCollector());
        l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);

    }
    else {
          l1Seeds_ = std::vector<std::string>();
          l1GtUtils_ = 0;
    }
}

ScoutingMuMuTreeMakerRun2::~ScoutingMuMuTreeMakerRun2() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void ScoutingMuMuTreeMakerRun2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;

  Handle<vector<ScoutingVertex> > verticesH;
  iEvent.getByToken(verticesToken, verticesH);

  if ( verticesH->size() == 0 ) return;

  Handle<vector<ScoutingMuon> > muonsH;
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
    trk_chi21 = muonsH->at(idx[0]).chi2();
    trk_chi22 = muonsH->at(idx[1]).chi2();

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

    //std::cout<<"pt: "<<pt1<<", "<<pt2<<std::endl;

    Handle<double> rhoH;
    iEvent.getByToken(rhoToken, rhoH);
    rho=*rhoH;

    Handle<vector<ScoutingVertex> > primaryVerticesH;
    iEvent.getByToken(primaryVerticesToken, primaryVerticesH);

    npvtx = 0;
    for (auto vtx_iter = primaryVerticesH->begin(); vtx_iter != primaryVerticesH->end(); ++vtx_iter) {
      //std::cout<<"primary x: "<<vtx_iter->x() <<" y: "<<  vtx_iter->y()<<" ex: "<<vtx_iter->xError()<<" ey: "<<vtx_iter->yError()<<std::endl;
      npvtx++;
    }

    //std::cout << "npvtx: " << npvtx << " avgPrimaryX: " << avgPrimary[0] << " avgPrimaryY: " << avgPrimary[1] << std::endl;

    double dx = ((verticesH->at(vtxIndx)).x()) - BSx;
    double dy = ((verticesH->at(vtxIndx)).y()) - BSy;

    vtxXError =(verticesH->at(vtxIndx)).xError();
    vtxYError =(verticesH->at(vtxIndx)).yError();
    vtxZError =(verticesH->at(vtxIndx)).zError();

    Lxy = sqrt(dx*dx + dy*dy);
    LxyErr = sqrt(dx*dx*vtxXError*vtxXError + dy*dy*vtxYError*vtxYError) / Lxy;
    LxySig = Lxy/LxyErr;

    //std::cout<<"Lxy: "<<Lxy<<" LxyErr: "<<LxyErr<<" LxySig: "<<LxySig<<std::endl;

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

    //std::cout<<"tree filling with mass: "<<mass<<", pt: "<<pt<<std::endl;
    tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void ScoutingMuMuTreeMakerRun2::beginJob() {
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

    tree->Branch("vtxXError"           , &vtxXError                   , "vtxXError/F");
    tree->Branch("vtxYError"           , &vtxYError                   , "vtxYError/F");
    tree->Branch("vtxZError"           , &vtxZError                   , "vtxZError/F");

    tree->Branch("l1Result", "std::vector<bool>"             ,&l1Result_, 32000, 0  );

}

// ------------ method called once each job just after ending the event loop  ------------
void ScoutingMuMuTreeMakerRun2::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingMuMuTreeMakerRun2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(ScoutingMuMuTreeMakerRun2);
