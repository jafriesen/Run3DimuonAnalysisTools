// -*- C++ -*-
//
// Package:    Run3ScoutingAnalysisTools/MuMuGammaTreeMaker
// Class:      MuMuGammaTreeMaker
//
/**\class MuMuGammaTreeMaker MuMuGammaTreeMaker.cc Run3ScoutingAnalysisTools/MuMuGammaTreeMaker/plugins/MuMuGammaTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/

// system include files
#include <memory>
#include <TTree.h>
#include <TLorentzVector.h>
#include "TMath.h"
#include <TPRegexp.h>

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

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//#include "TMtuple.hh"
#include "MMGUtils.hh"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class MuMuGammaTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit MuMuGammaTreeMaker(const edm::ParameterSet&);
  ~MuMuGammaTreeMaker() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;   

  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>             triggerResultsToken;
  const edm::EDGetTokenT<std::vector<pat::Muon> >      muonsToken;
  const edm::EDGetTokenT<std::vector<pat::Electron> >  electronsToken;
  const edm::EDGetTokenT<std::vector<reco::Vertex> >    primaryVerticesToken;
  const edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate> >    verticesToken;
  const edm::EDGetTokenT<double>                          rhoToken;
  const edm::EDGetTokenT<std::vector<pat::Photon> >         photonsToken;
  const edm::EDGetTokenT<std::vector<pat::PackedCandidate> >         pfCandsToken;
  const edm::EDGetTokenT<std::vector<reco::GenParticle> >         prunedGenToken;
  const edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >    packedGenToken;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> esToken;

  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;

  bool doL1;
  bool doGEN;
  triggerExpression::Data triggerCache_;

  edm::InputTag                algInputTag_;
  edm::InputTag                extInputTag_;
  edm::EDGetToken              algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;
  std::vector<bool>            l1Result_;
  std::vector<bool>            hltResult_;

  TTree* tree;

  int eventNum;
  int lumiSec;
  int runNum;

  float pfIso1;
  float pfIso2;

  float mass;
  float pt;
  float dr;
  float pt1;
  float pt2;

  float eta1;
  float eta2;
  float phi1;
  float phi2;

  float dxy1;
  float dxy2;
  float dz1;
  float dz2;
  float trkChi21;
  float trkChi22;
  float trkNdof1;
  float trkNdof2;

  // MVA Variables
  float segmentCompatibility[2];
  float chi2LocalMomentum[2];
  float chi2LocalPosition[2];
  float glbTrackProbability[2];
  float iValidFraction[2];
  float layersWithMeasurement[2];
  float trkKink[2];
  float log2PlusGlbKink[2];
  float timeAtIpInOutErr[2];
  float outerChi2[2];
  float innerChi2[2];
  float trkRelChi2[2];
  float vMuonHitComb[2];
  float qProd[2];
  float mva[2];

  float probVtx;
  float vtxX;
  float vtxY;
  float vtxZ;
  float vtxXError;
  float vtxYError;
  float vtxZError;
  float vtx_chi2;

  int npv;
  float pvX;
  float pvY;
  float pvZ;

  int nPhotons;

  std::vector<bool> muonID1;
  std::vector<bool> muonID2;

  std::vector<float> slimmedPhotonPt;
  std::vector<float> slimmedPhotonEta;
  std::vector<float> slimmedPhotonPhi;
  std::vector<float> slimmedPhotonM;
  std::vector<float> slimmedPhotonSigmaIetaIeta;
  std::vector<float> slimmedPhotonHOverE;
  std::vector<float> slimmedPhotonEcalIso;
  std::vector<float> slimmedPhotonHcalIso;
  std::vector<float> slimmedPhotonTrkIso;
  std::vector<float> slimmedPhotonR9;

  std::vector<float> pfCandPhotonDr;
  std::vector<float> pfCandPhotonIso;
  std::vector<float> pfCandPhotonPt;
  std::vector<float> pfCandPhotonEta;
  std::vector<float> pfCandPhotonPhi;
  std::vector<float> pfCandPhotonEnergy;
  std::vector<float> pfCandPhotonEt;
  std::vector<float> pfCandPhotonEt2;

  int motherID1;
  int motherID2;

  int motherGenID;
  int mupGenID;
  int mumGenID;

  std::vector<int> matchedDaughtersIDs;
  float matchedPhotonPt;
  float matchedPhotonEta;
  float matchedPhotonPhi;

  float mu1_id;
  float mu2_id;

  // flags for GEN matched decays
  bool isEta2MuMu;
  bool isEta2MuMuGamma;
  bool isEtaPrime2MuMu;
  bool isEtaPrime2MuMuGamma;
  bool isOmega2MuMu;
  bool isOmega2Pi0MuMu;
  bool isRho2MuMu;
  // bool isRho2Pi0MuMu; //not observed?
  bool isPhi2MuMu;
  bool isPhi2KK;
    
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
MuMuGammaTreeMaker::MuMuGammaTreeMaker(const edm::ParameterSet& iConfig):
    triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
    triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),
    muonsToken               (consumes<std::vector<pat::Muon> >             (iConfig.getParameter<edm::InputTag>("muons"))),
    electronsToken           (consumes<std::vector<pat::Electron> >         (iConfig.getParameter<edm::InputTag>("electrons"))),
    primaryVerticesToken     (consumes<std::vector<reco::Vertex> >           (iConfig.getParameter<edm::InputTag>("primaryVertices"))),
    verticesToken            (consumes<std::vector<reco::VertexCompositePtrCandidate> >           (iConfig.getParameter<edm::InputTag>("displacedVertices"))),
    //rhoToken                 (consumes<double>                                 (iConfig.getParameter<edm::InputTag>("rho"))), 
    photonsToken             (consumes<std::vector<pat::Photon> >         (iConfig.getParameter<edm::InputTag>("photons"))),
    pfCandsToken             (consumes<std::vector<pat::PackedCandidate> >         (iConfig.getParameter<edm::InputTag>("pfcands"))),
    prunedGenToken  (consumes<std::vector<reco::GenParticle> >      (iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
    packedGenToken  (consumes<std::vector<pat::PackedGenParticle> > (iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
    esToken(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
    doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1")            : false),
    doGEN                    (iConfig.existsAs<bool>("doGEN")              ?    iConfig.getParameter<bool>  ("doGEN")           : false),
    doFullGEN                (iConfig.existsAs<bool>("doFullGEN")          ?    iConfig.getParameter<bool>  ("doFullGEN")       : false)
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

MuMuGammaTreeMaker::~MuMuGammaTreeMaker() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void MuMuGammaTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;

  if (doFullGEN) {

    Handle<vector<reco::GenParticle> > prunedGenParticles;
    iEvent.getByToken(prunedGenToken, prunedGenParticles);

    Handle<vector<pat::PackedGenParticle> > packedGenParticles;
    iEvent.getByToken(packedGenToken, packedGenParticles);

    float muGenDr[2] = {999., 999.};
    TLorentzVector mupGen, mumGen, daught3Gen;


    for (auto genp = prunedGenParticles->begin(); genp != prunedGenParticles->end(); ++genp) {            
      if (abs(genp->pdgId())==221 or abs(genp->pdgId())==113 or abs(genp->pdgId())==223 or abs(genp->pdgId())==331 or abs(genp->pdgId()==333)) {
        // std::cout<<"genp id: "<<genp->pdgId()<<" pt: "<<genp->pt()<<" eta: "<<genp->eta()<<" status: "<<genp->status();
        int nDaughterMuons = 0;
        int daughterMuIdx[2] = {-1, -1}; // mup, mum

        for (int i=0; i<(int)genp->numberOfDaughters(); ++i) {
          // check if muons and save idx of daughters
          if ( abs(genp->daughter(i)->pdgId()) == 13 and genp->daughter(i)->pt() > 2 and abs(genp->daughter(i)->eta()) < 2.5 ){
            nDaughterMuons+=1;
            daughterMuIdx[genp->daughter(i)->pdgId() > 0] = i;
          }
        }
        
        if (nDaughterMuons==2) {
          auto daught_mup =  &(*genp->daughter(daughterMuIdx[0]));
          auto daught_mum =  &(*genp->daughter(daughterMuIdx[1]));

          daught_mup_pt = daught_mup->pt();
          daught_mup_eta = daught_mup->eta();
          daught_mup_phi = daught_mup->phi();
          daught_mum_pt = daught_mum->pt();
          daught_mum_eta = daught_mum->eta();
          daught_mum_phi = daught_mum->phi();

          mup_gen.SetPtEtaPhiM( daught_mup->pt(), daught_mup->eta(), daught_mup->phi(), mu_mass);
          mum_gen.SetPtEtaPhiM( daught_mum->pt(), daught_mum->eta(), daught_mum->phi(), mu_mass);

          mass_gen = (mup_gen+mum_gen).M();

          matchedDaughtersIDs.clear();
          std::cout<<"Mother with 2 muons ID == " << genp->pdgId() <<std::endl;
          for (int i=0; i<(int)genp->numberOfDaughters(); ++i) {
            matchedDaughtersIDs.push_back(genp->daughter(i)->pdgId());
            std::cout<<"Daughter "<< i << "  ID == " << genp->daughter(i)->pdgId()<<std::endl;
          }

          //check for photons
          bool aPhoton = false;
          for (auto pgp = packedGenParticles->begin(); pgp != packedGenParticles->end(); ++pgp) {

            if (pgp->pdgId()==22 and pgp->motherRef()->pdgId()==genp->pdgId() and
                abs((pgp->motherRef()->pt() - genp->pt()) / genp->pt()) < 0.01 and
                abs((pgp->motherRef()->eta() - genp->eta()) / genp->eta()) < 0.01 and
                abs((pgp->motherRef()->phi() - genp->phi()) / genp->phi()) < 0.01
            ) {
              std::cout<<"packed photon "<<pgp->pt()<<" mother pt: "<<pgp->motherRef()->pt()<<std::endl;

              matchedPhotonPt  = pgp->pt();
              matchedPhotonEta = pgp->eta();
              matchedPhotonPhi = pgp->phi();
              aPhoton = true;
            }
          }
          // isEta2MuMu
          // isEta2MuMuGamma
          // isEtaPrime2MuMu
          // isEtaPrime2MuMuGamma
          // isOmega2MuMu
          // isOmega2Pi0MuMu

          // isRho2MuMu
          // isPhi2MuMu
          // isPhi2KK
          if      ((abs(genp->pdgId()) == 221) and !(aPhoton)) isEta2MuMu = true;
          else if ((abs(genp->pdgId()) == 221) and aPhoton)    isEta2MuMuGamma = true;
          else if ((abs(genp->pdgId()) == 331) and !(aPhoton)) isEtaPrime2MuMu = true;
          else if ((abs(genp->pdgId()) == 331) and aPhoton)    isEtaPrime2MuMuGamma = true;
          else if ((abs(genp->pdgId()) == 223) and !(aPhoton)) isOmega2MuMu = true;
          // else if ((abs(genp->pdgId()) == 223) and aPhoton)    isOmega2Pi0MuMu = true;    // this one is not so easy     
          else if ((abs(genp->pdgId()) == 113) and !(aPhoton)) isRho2MuMu = true;
          else if ((abs(genp->pdgId()) == 333) and !(aPhoton)) isPhi2MuMu = true;
        }
      }
    }    
  }

  
  Handle<vector<reco::Vertex> > primaryVerticesH;
  iEvent.getByToken(primaryVerticesToken, primaryVerticesH);
  npv = primaryVerticesH->size();
  reco::Vertex pv = *primaryVerticesH->begin();

  if ( npv == 0 ) return;

  Handle<vector<pat::Muon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);

  if ( muonsH->size() < 2 ) return;

  bool fillTree = false;
  int idx[2];

  //std::cout << std::endl;
  //std::cout << "primary vtx 1: " << pv.x() << " " << pv.y() << " " << pv.z() << std::endl;

  /*
    int charge1 = (muonsH->at(0)).charge();
    idx[0] = 0;

    for ( unsigned int iMuon = 0; iMuon < muonsH->size(); iMuon++ ) {
      if ( ((muonsH->at(iMuon)).charge() * charge1) < 0 ) {
        //std::cout << "idx[1] " << iMuon << std::endl;
        idx[1] = iMuon;
        fillTree = true;
        break;
      }
    }
  */
  edm::ESHandle<TransientTrackBuilder> theB = iSetup.getHandle(esToken);
  KalmanVertexFitter kvf(true);
  TransientVertex tv;
  float bestProbVtx = 0.005; // loose minimum probability as in dimuon HLT path

  reco::Vertex vertex;
  for ( unsigned int iMuon = 0; iMuon < muonsH->size(); iMuon++ ) {
    if ( (muonsH->at(iMuon)).pt() > 3 ) {
      for ( unsigned int jMuon = iMuon+1; jMuon < muonsH->size(); jMuon++ ) {
        if ( (muonsH->at(jMuon)).pt() > 3 && ( (muonsH->at(iMuon)).charge() * (muonsH->at(jMuon)).charge() < 0 ) ) {
          reco::Track part_1, part_2;
          part_1 = *((muonsH->at(iMuon)).bestTrack());
          part_2 = *((muonsH->at(jMuon)).bestTrack());
          vector<reco::TransientTrack> transient_tracks{};
          transient_tracks.push_back(theB->build(fix_track(&part_1)));
          transient_tracks.push_back(theB->build(fix_track(&part_2)));
          tv = kvf.vertex(transient_tracks);

          if (!tv.isValid()) {
            //std::cout << "ij " << iMuon << jMuon << "Vertex not valid." << std::endl;
          } else {
            vertex = reco::Vertex(tv);
            float vertexProb = TMath::Prob( vertex.chi2() , vertex.ndof() );

            if (vertexProb > bestProbVtx) {
              fillTree = true;
              idx[0] = iMuon;
              idx[1] = jMuon;
              bestProbVtx = vertexProb;
            }

            //vtxX = vertex.x();
            //vtxY = vertex.y();
            //vtx_chi2 = vertex.normalizedChi2();
            //vtxZ = vertex.z();

            //std::cout << "ij " << iMuon << jMuon << " vtx: " << vtxX << " " << vtxY << " " << vtxZ << " chi2 " << vertex.chi2() << " ndof " << vertex.ndof() << " probVtx " << probVtx << " normalizedChi2 " << vtx_chi2 <<std::endl;
          }
        }
      }
    }
  }

  //if (fillTree) {
  //    for (auto vtx_iter = primaryVerticesH->begin(); vtx_iter != primaryVerticesH->end(); ++vtx_iter) {
  //      std::cout<<"primary vtx: "<<vtx_iter->x() <<" "<<  vtx_iter->y()<<" "<<  vtx_iter->z() << " probVtx " << TMath::Prob( vtx_iter->chi2() , vtx_iter->ndof() ) << " normalizedChi2 " << vtx_iter->normalizedChi2() << std::endl;
  //    }
  //}

  if (fillTree) {
    eventNum = iEvent.id().event();
    lumiSec = iEvent.luminosityBlock();
    runNum = iEvent.id().run();

    //std::cout <<idx[0]<< " muon1 vtx: " << (muonsH->at(idx[0])).vx() << " " << (muonsH->at(idx[0])).vy() << " " << (muonsH->at(idx[0])).vz() << std::endl;
    //std::cout <<idx[1]<< " muon2 vtx: " << (muonsH->at(idx[1])).vx() << " " << (muonsH->at(idx[1])).vy() << " " << (muonsH->at(idx[1])).vz() << std::endl;
    //std::cout << "charge " << (muonsH->at(idx[0])).charge() << (muonsH->at(idx[1])).charge() << std::endl;

    // Compute the pfIso for the muon. Note: PUCorr = 0.5*muons_iter->puChargedHadronIso()                                                                              
    // -----------------------------------------------------------------------------------                                                                              
    pfIso1 = (muonsH->at(idx[0]).chargedHadronIso() + std::max(muonsH->at(idx[0]).photonIso() + muonsH->at(idx[0]).neutralHadronIso() - 0.5*muonsH->at(idx[0]).puChargedHadronIso(), 0.0))/muonsH->at(idx[0]).pt();
    pfIso2 = (muonsH->at(idx[1]).chargedHadronIso() + std::max(muonsH->at(idx[1]).photonIso() + muonsH->at(idx[1]).neutralHadronIso() - 0.5*muonsH->at(idx[1]).puChargedHadronIso(), 0.0))/muonsH->at(idx[1]).pt();

    pt1=muonsH->at(idx[0]).pt();
    pt2=muonsH->at(idx[1]).pt();

    eta1=muonsH->at(idx[0]).eta();
    eta2=muonsH->at(idx[1]).eta();
    phi1=muonsH->at(idx[0]).phi();
    phi2=muonsH->at(idx[1]).phi();

    dxy1 = muonsH->at(idx[0]).muonBestTrack()->dxy();
    dxy2 = muonsH->at(idx[1]).muonBestTrack()->dxy();
    dz1 = muonsH->at(idx[0]).muonBestTrack()->dz();
    dz2 = muonsH->at(idx[1]).muonBestTrack()->dz();
    trkChi21 = muonsH->at(idx[0]).muonBestTrack()->chi2();
    trkChi22 = muonsH->at(idx[1]).muonBestTrack()->chi2();
    trkNdof1 = muonsH->at(idx[0]).muonBestTrack()->ndof();
    trkNdof2 = muonsH->at(idx[1]).muonBestTrack()->ndof();

    for (int i : idx) {
      chi2LocalMomentum[i] = muonsH->at(idx[i]).combinedQuality().chi2LocalMomentum;
      chi2LocalPosition[i] = muonsH->at(idx[i]).combinedQuality().chi2LocalPosition;
      glbTrackProbability[i] = muonsH->at(idx[i]).combinedQuality().glbTrackProbability;

      trkKink[i] = muonsH->at(idx[i]).combinedQuality().trkKink;
      log2PlusGlbKink[i] = TMath::Log(2 + muonsH->at(idx[i]).combinedQuality().glbKink);
      trkRelChi2[i] = muonsH->at(idx[i]).combinedQuality().trkRelChi2;
      segmentCompatibility[i] = muonsH->at(idx[i]).segmentCompatibility();

      timeAtIpInOutErr[i] = muonsH->at(idx[i]).time().timeAtIpInOutErr;

      reco::TrackRef gTrack = muonsH->at(idx[i]).globalTrack();
      reco::TrackRef iTrack = muonsH->at(idx[i]).innerTrack();
      reco::TrackRef oTrack = muonsH->at(idx[i]).outerTrack();

      if (!(iTrack.isNonnull() and oTrack.isNonnull() and gTrack.isNonnull())) {
        std::cout << "event " << eventNum << ": null track" << std::endl;
        iValidFraction[i] = -1;
        innerChi2[i] = -1;
        layersWithMeasurement[i] = -1;
        outerChi2[i] = -1;
        qProd[i] = -1;
        vMuonHitComb[i] = -1;
        mva[i] = -1;
      }
      else {
        iValidFraction[i] = iTrack->validFraction();
        innerChi2[i] = iTrack->normalizedChi2();
        layersWithMeasurement[i] = iTrack->hitPattern().trackerLayersWithMeasurement();
        outerChi2[i] = oTrack->normalizedChi2();

        qProd[i] = iTrack->charge() * oTrack->charge();

        vMuonHitComb[i] = validMuonHitComb( muonsH->at(idx[i]) );

        mva[i] = muonsH->at(idx[i]).softMvaValue();
      }
    }

    //std::cout<<dxy1<<" "<<dz1<<" "<<trkChi21<<" "<<trkNdof1<<std::endl;
    
    TLorentzVector mu1;
    mu1.SetPtEtaPhiM(pt1,eta1,phi1,mu_mass);

    TLorentzVector mu2;
    mu2.SetPtEtaPhiM(pt2,eta2,phi2,mu_mass);

    TLorentzVector dimu = mu1+mu2;
    mass=dimu.M();
    pt=dimu.Pt();
    dr=mu1.DeltaR(mu2);

    motherID1=0; motherID2=0; 
    simType1=999; simType2=999;
    simExtType1=999; simExtType2=999;
    
    mu1_id        = (muonsH->at(idx[0])).pdgId();
    mu2_id        = (muonsH->at(idx[1])).pdgId();

    if (doGEN) {
        motherID1=muonsH->at(idx[0]).simMotherPdgId(); motherID2=muonsH->at(idx[1]).simMotherPdgId();
        simType1=muonsH->at(idx[0]).simType(); simType2=muonsH->at(idx[1]).simType();
        simExtType1=muonsH->at(idx[0]).simExtType(); simExtType2=muonsH->at(idx[1]).simExtType();
    }
    
    /*
    std::cout<<"pt: "<<pt1<<", "<<pt2<<" eta: "<<eta1<<", "<<eta2<<" phi: "<<std::endl;
    std::cout<<"motherID1: "<<muonsH->at(idx[0]).simMotherPdgId()<<" motherID2: "<<muonsH->at(idx[1]).simMotherPdgId()<<std::endl;
    std::cout<<"simType1: "<<muonsH->at(idx[0]).simType()<<" simType2: "<<muonsH->at(idx[1]).simType()<<std::endl;
    std::cout<<"simExtType1: "<<muonsH->at(idx[0]).simExtType()<<" simExtType2: "<<muonsH->at(idx[1]).simExtType()<<std::endl;
    */
    //std::cout<<"pfIso: "<<pfIso1<<", "<<pfIso2<<std::endl;
    //std::cout << "id1 " << (muonsH->at(idx[0])).isHighPtMuon(pv) << (muonsH->at(idx[0])).isLooseMuon() << (muonsH->at(idx[0])).isMediumMuon() << (muonsH->at(idx[0])).isSoftMuon(pv) << (muonsH->at(idx[0])).isTightMuon(pv) << std::endl;
    //std::cout << "id2 " << (muonsH->at(idx[1])).isHighPtMuon(pv) << (muonsH->at(idx[1])).isLooseMuon() << (muonsH->at(idx[1])).isMediumMuon() << (muonsH->at(idx[1])).isSoftMuon(pv) << (muonsH->at(idx[1])).isTightMuon(pv) << std::endl;

    muonID1.clear();
    muonID2.clear();

    muonID1.push_back((muonsH->at(idx[0])).isHighPtMuon(pv));
    muonID1.push_back((muonsH->at(idx[0])).isLooseMuon());
    muonID1.push_back((muonsH->at(idx[0])).isMediumMuon());
    muonID1.push_back((muonsH->at(idx[0])).isSoftMuon(pv));
    muonID1.push_back((muonsH->at(idx[0])).isTightMuon(pv));

    muonID2.push_back((muonsH->at(idx[1])).isHighPtMuon(pv));
    muonID2.push_back((muonsH->at(idx[1])).isLooseMuon());
    muonID2.push_back((muonsH->at(idx[1])).isMediumMuon());
    muonID2.push_back((muonsH->at(idx[1])).isSoftMuon(pv));
    muonID2.push_back((muonsH->at(idx[1])).isTightMuon(pv));

    pvX = pv.x();
    pvY = pv.y();
    pvZ = pv.z();

    vtxX = vertex.x();
    vtxY = vertex.y();
    vtxZ = vertex.z();
    vtx_chi2 = vertex.normalizedChi2();
    vtxXError = vertex.xError();
    vtxYError = vertex.yError();
    vtxZError = vertex.zError();
    probVtx = bestProbVtx;

    //Handle<double> rhoH;
    //iEvent.getByToken(rhoToken, rhoH);
    //rho=*rhoH;

    Handle<vector<pat::Photon> > photonsH;
    iEvent.getByToken(photonsToken, photonsH);

    slimmedPhotonPt.clear();
    slimmedPhotonEta.clear();
    slimmedPhotonPhi.clear();
    slimmedPhotonM.clear();
    slimmedPhotonSigmaIetaIeta.clear();
    slimmedPhotonHOverE.clear();
    slimmedPhotonEcalIso.clear();
    slimmedPhotonHcalIso.clear();
    slimmedPhotonTrkIso.clear();
    slimmedPhotonR9.clear();
    //slimmedPhotonIso.clear();

    nPhotons = 0;
    for (auto photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
      //std::cout << "  slimmedPhotonPt " <<  photons_iter->pt() << " slimmedPhotonEta " << photons_iter->eta() << " slimmedPhotonPhi " << photons_iter->phi() << " isPFlowPhoton " << photons_iter->isPFlowPhoton() << std::endl;
      slimmedPhotonPt.push_back(photons_iter->pt());
      slimmedPhotonEta.push_back(photons_iter->eta());
      slimmedPhotonPhi.push_back(photons_iter->phi());
      slimmedPhotonM.push_back(photons_iter->mass());
      slimmedPhotonSigmaIetaIeta.push_back(photons_iter->sigmaIetaIeta());
      slimmedPhotonHOverE.push_back(photons_iter->hadronicOverEm());
      slimmedPhotonEcalIso.push_back(photons_iter->ecalIso());
      slimmedPhotonHcalIso.push_back(photons_iter->hcalIso());
      slimmedPhotonTrkIso.push_back(photons_iter->trackIso());
      slimmedPhotonR9.push_back(photons_iter->r9());

      nPhotons++;
    }

    Handle<vector<pat::PackedCandidate> > pfCandH;
    iEvent.getByToken(pfCandsToken, pfCandH);

    pfCandPhotonDr.clear();
    pfCandPhotonIso.clear();
    pfCandPhotonPt.clear();
    pfCandPhotonEta.clear();
    pfCandPhotonPhi.clear();
    pfCandPhotonEnergy.clear();
    pfCandPhotonEt.clear();
    pfCandPhotonEt2.clear();

    for (auto pfCand_iter = pfCandH->begin(); pfCand_iter != pfCandH->end(); ++pfCand_iter) {
      if (pfCand_iter->pdgId() != 22) continue;
      if (pfCand_iter->pt() < 1.0) continue;

      float pfCandDimuDr = deltaR(dimu.Eta(), dimu.Phi(), pfCand_iter->eta(), pfCand_iter->phi());
      if (pfCandDimuDr < 0.5) {
        pfCandPhotonDr.push_back(pfCandDimuDr);
        pfCandPhotonIso.push_back(photonPfIso03(*pfCand_iter,pfCandH)/pfCand_iter->pt());
        pfCandPhotonPt.push_back(pfCand_iter->pt());
        pfCandPhotonEta.push_back(pfCand_iter->eta());
        pfCandPhotonPhi.push_back(pfCand_iter->phi());
        pfCandPhotonEnergy.push_back(pfCand_iter->energy());
        pfCandPhotonEt.push_back(pfCand_iter->et());
        pfCandPhotonEt2.push_back(pfCand_iter->et2());
        //std::cout << " photon candidate: photonDr " << pfCandPhotonDr.back() << " photonIso " << pfCandPhotonIso.back() << " photonPt " << pfCandPhotonPt.back() << " photonEta " << pfCandPhotonEta.back() << " photonPhi " << pfCandPhotonPhi.back() << " photonEnergy " << pfCandPhotonEnergy.back() << " photonEt " << pfCandPhotonEt.back() << " photonEt2 " << pfCandPhotonEt2.back() << std::endl;
      }
    }

    //std::cout << slimmedPhotonPt.size() << " photons" << std::endl;

    l1Result_.clear();
    if (doL1) {
        l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
        for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
            bool l1htbit = 0;
            l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
            l1Result_.push_back( l1htbit );
        }
    }

    Handle<TriggerResults> triggerResultsH;
    iEvent.getByToken(triggerResultsToken, triggerResultsH);
    hltResult_.clear();
    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        hltResult_.push_back(triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]]));
    }

    motherGenID = 0; 
    mupGenID = 0; mumGenID = 0;
    isEta2MuMu            = false;
    isEta2MuMuGamma       = false;
    isEtaPrime2MuMu       = false;
    isEtaPrime2MuMuGamma  = false;
    isOmega2MuMu         = false;
    isOmega2Pi0MuMu       = false;
    isRho2MuMu            = false;
    isPhi2MuMu            = false;
    isPhi2KK              = false;

    if (doGEN and (mass < 2.0)) {

        Handle<vector<reco::GenParticle> > prunedGenParticles;
        iEvent.getByToken(prunedGenToken, prunedGenParticles);

        Handle<vector<pat::PackedGenParticle> > packedGenParticles;
        iEvent.getByToken(packedGenToken, packedGenParticles);

        float dr_mup_gen = 999.;
        float dr_mum_gen = 999.;
        TLorentzVector mup_gen, mum_gen, daught3_gen; 
	      TLorentzVector *mup_reco;
	      TLorentzVector *mum_reco;
        if( mu1_id == 13) {
          mum_reco = &mu1;
          mup_reco = &mu2;
        }
        else {
          mum_reco = &mu2;
          mup_reco = &mu1;
        }

        for (auto genp = prunedGenParticles->begin(); genp != prunedGenParticles->end(); ++genp) {            
            if (abs(genp->pdgId())==221 or abs(genp->pdgId())==113 or abs(genp->pdgId())==223 or abs(genp->pdgId())==331 or abs(genp->pdgId()==333)) {
                // std::cout<<"genp id: "<<genp->pdgId()<<" pt: "<<genp->pt()<<" eta: "<<genp->eta()<<" status: "<<genp->status();
                    int nDaughterMuons = 0;
                    int mup_idx = 99;
                    int mum_idx = 99;
                    // int daught3_idx = 99;

                    for (int i=0; i<(int)genp->numberOfDaughters(); ++i) {

                      // check if muons and save idx of daughters
                      if (genp->daughter(i)->pdgId()==13){
                        mum_idx = i;
                        nDaughterMuons+=1;
                      } else if  (genp->daughter(i)->pdgId()== -13){
                        mup_idx = i;
                        nDaughterMuons+=1;
                      }
            
                    }
                    if (nDaughterMuons==2) {
                      // try to match pair of reco muons with a pair of daughters
                      auto daught_mup =  &(*genp->daughter(mup_idx));
                      auto daught_mum =  &(*genp->daughter(mum_idx));

                      mup_gen.SetPtEtaPhiM( daught_mup->pt(), daught_mup->eta(), daught_mup->phi(),mu_mass);
                      mum_gen.SetPtEtaPhiM( daught_mum->pt(), daught_mum->eta(), daught_mum->phi(),mu_mass);
                      dr_mup_gen = mup_gen.DeltaR(*mup_reco);
                      dr_mum_gen = mum_gen.DeltaR(*mum_reco);

                      matchedDaughtersIDs.clear();
                      std::cout<<"Mother with 2 muons ID == " << genp->pdgId() <<std::endl;
                      for (int i=0; i<(int)genp->numberOfDaughters(); ++i) {
                        matchedDaughtersIDs.push_back(genp->daughter(i)->pdgId());
                        std::cout<<"Daughter "<< i << "  ID == " << genp->daughter(i)->pdgId()<<std::endl;
                      }

                      if ((dr_mup_gen < MIN_DR_TRUTH) and (dr_mum_gen < MIN_DR_TRUTH)){
                        motherGenID = genp->pdgId();
                        mupGenID = daught_mup->pdgId();
                        mumGenID = daught_mum->pdgId();
                        std::cout<<"Matched "<< daught_mup->pt() << "  " << mup_reco->Pt() << "  " << daught_mum->pt() << "  " << mum_reco->Pt() << "for mother ID == " << motherGenID <<std::endl;
                        
                      }
                      
                      //check for photons
                      bool aPhoton = false;
                      for (auto pgp = packedGenParticles->begin(); pgp != packedGenParticles->end(); ++pgp) {
                            if (pgp->pdgId()==22 and pgp->motherRef()->pdgId()==genp->pdgId()) {
                                std::cout<<"packed photon "<<pgp->pt()<<" mother pt: "<<pgp->motherRef()->pt()<<std::endl;
                                matchedPhotonPt  = pgp->pt();
                                matchedPhotonEta = pgp->eta();
                                matchedPhotonPhi = pgp->phi();
                                aPhoton = true;
                                
                            }
                        }
                      // isEta2MuMu
                      // isEta2MuMuGamma
                      // isEtaPrime2MuMu
                      // isEtaPrime2MuMuGamma
                      // isOmega2MuMu
                      // isOmega2Pi0MuMu

                      // isRho2MuMu
                      // isPhi2MuMu
                      // isPhi2KK
                      if      ((abs(genp->pdgId()) == 221) and !(aPhoton)) isEta2MuMu = true;
                      else if ((abs(genp->pdgId()) == 221) and aPhoton)    isEta2MuMuGamma = true;
                      else if ((abs(genp->pdgId()) == 331) and !(aPhoton)) isEtaPrime2MuMu = true;
                      else if ((abs(genp->pdgId()) == 331) and aPhoton)    isEtaPrime2MuMuGamma = true;
                      else if ((abs(genp->pdgId()) == 223) and !(aPhoton)) isOmega2MuMu = true;
                      // else if ((abs(genp->pdgId()) == 223) and aPhoton)    isOmega2Pi0MuMu = true;    // this one is not so easy     
                      else if ((abs(genp->pdgId()) == 113) and !(aPhoton)) isRho2MuMu = true;
                      else if ((abs(genp->pdgId()) == 333) and !(aPhoton)) isPhi2MuMu = true;
                    }
            }
        }    
    }
    tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void MuMuGammaTreeMaker::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"      , "tree");

    tree->Branch("eventNum"            , &eventNum                    , "eventNum/I");
    tree->Branch("lumiSec"             , &lumiSec                     , "lumiSec/I");
    tree->Branch("runNum"              , &runNum                      , "runNum/I");

    tree->Branch("mass"                , &mass                        , "mass/F"    );
    tree->Branch("pt"                  , &pt                          , "pt/F"      );
    tree->Branch("dr"                  , &dr                          , "dr/F"      );
    tree->Branch("pt1"                 , &pt1                         , "pt1/F"     );
    tree->Branch("pt2"                 , &pt2                         , "pt2/F"     );
    tree->Branch("eta1"                , &eta1                        , "eta1/F"    );
    tree->Branch("eta2"                , &eta2                        , "eta2/F"    );
    tree->Branch("phi1"                , &phi1                        , "phi1/F"    );
    tree->Branch("phi2"                , &phi2                        , "phi2/F"    );
    tree->Branch("pfIso1"              , &pfIso1                      , "pfIso1/F"  );
    tree->Branch("pfIso2"              , &pfIso2                      , "pfIso2/F"  );

    tree->Branch("segmentCompatibility1",     &segmentCompatibility[0],     "segmentCompatibility1/F"   );
    tree->Branch("chi2LocalMomentum1",        &chi2LocalMomentum[0],        "chi2LocalMomentum1/F"      );
    tree->Branch("chi2LocalPosition1",        &chi2LocalPosition[0],        "chi2LocalPosition1/F"      );
    tree->Branch("glbTrackProbability1",      &glbTrackProbability[0],      "glbTrackProbability1/F"    );
    tree->Branch("iValidFraction1",           &iValidFraction[0],           "iValidFraction1/F"         );
    tree->Branch("layersWithMeasurement1",    &layersWithMeasurement[0],    "layersWithMeasurement1/F"  );
    tree->Branch("trkKink1",                  &trkKink[0],                  "trkKink1/F"                );
    tree->Branch("log2PlusGlbKink1",          &log2PlusGlbKink[0],          "log2PlusGlbKink1/F"        );
    tree->Branch("timeAtIpInOutErr1",         &timeAtIpInOutErr[0],         "timeAtIpInOutErr1/F"       );
    tree->Branch("outerChi21",                &outerChi2[0],                "outerChi21/F"              );
    tree->Branch("innerChi21",                &innerChi2[0],                "innerChi21/F"              );
    tree->Branch("trkRelChi21",               &trkRelChi2[0],               "trkRelChi21/F"             );
    tree->Branch("vMuonHitComb1",             &vMuonHitComb[0],             "vMuonHitComb1/F"           );
    tree->Branch("qProd1",                    &qProd[0],                    "qProd1/F"                  );
    tree->Branch("mva1",                      &mva[0],                      "mva1/F"                    );

    tree->Branch("segmentCompatibility2",     &segmentCompatibility[1],     "segmentCompatibility2/F"   );
    tree->Branch("chi2LocalMomentum2",        &chi2LocalMomentum[1],        "chi2LocalMomentum2/F"      );
    tree->Branch("chi2LocalPosition2",        &chi2LocalPosition[1],        "chi2LocalPosition2/F"      );
    tree->Branch("glbTrackProbability2",      &glbTrackProbability[1],      "glbTrackProbability2/F"    );
    tree->Branch("iValidFraction2",           &iValidFraction[1],           "iValidFraction2/F"         );
    tree->Branch("layersWithMeasurement2",    &layersWithMeasurement[1],    "layersWithMeasurement2/F"  );
    tree->Branch("trkKink2",                  &trkKink[1],                  "trkKink2/F"                );
    tree->Branch("log2PlusGlbKink2",          &log2PlusGlbKink[1],          "log2PlusGlbKink2/F"        );
    tree->Branch("timeAtIpInOutErr2",         &timeAtIpInOutErr[1],         "timeAtIpInOutErr2/F"       );
    tree->Branch("outerChi22",                &outerChi2[1],                "outerChi22/F"              );
    tree->Branch("innerChi22",                &innerChi2[1],                "innerChi22/F"              );
    tree->Branch("trkRelChi22",               &trkRelChi2[1],               "trkRelChi22/F"             );
    tree->Branch("vMuonHitComb2",             &vMuonHitComb[1],             "vMuonHitComb2/F"           );
    tree->Branch("qProd2",                    &qProd[1],                    "qProd2/F"                  );
    tree->Branch("mva2",                      &mva[1],                      "mva2/F"                    );

    tree->Branch("dxy1"              , &dxy1                      , "dxy1/F"  );
    tree->Branch("dxy2"              , &dxy2                      , "dxy2/F"  );
    tree->Branch("dz1"              , &dz1                      , "dz1/F"  );
    tree->Branch("dz2"              , &dz2                      , "dz2/F"  );
    tree->Branch("trkChi21"              , &trkChi21                      , "trkChi21/F"  );
    tree->Branch("trkChi22"              , &trkChi22                      , "trkChi22/F"  );
    tree->Branch("trkNdof1"              , &trkNdof1                      , "trkNdof1/F"  );
    tree->Branch("trkNdof2"              , &trkNdof2                      , "trkNdof2/F"  );

    tree->Branch("motherID1"              , &motherID1                      , "motherID1/I"  );
    tree->Branch("motherID2"              , &motherID2                      , "motherID2/I"  );
    tree->Branch("simType1"              , &simType1                      , "simType1/I"  );
    tree->Branch("simType2"              , &simType2                      , "simType2/I"  );
    tree->Branch("simExtType1"              , &simExtType1                      , "simExtType1/I"  );
    tree->Branch("simExtType2"              , &simExtType2                      , "simExtType2/I"  );

    tree->Branch("matchedDaughtersIDs", "std::vector<int>", &matchedDaughtersIDs, 32000, 0);
    tree->Branch("matchedPhotonPt"              , &matchedPhotonPt                      , "matchedPhotonPt/F"  );
    tree->Branch("matchedPhotonEta"              , &matchedPhotonEta                      , "matchedPhotonEta/F"  );
    tree->Branch("matchedPhotonPhi"              , &matchedPhotonPhi                      , "matchedPhotonPhi/F"  );

    tree->Branch("isEta2MuMu",               &isEta2MuMu,             "isEta2MuMu/b");    
    tree->Branch("isEta2MuMuGamma",          &isEta2MuMuGamma,            "isEta2MuMuGamma/b");
    tree->Branch("isEtaPrime2MuMu",          &isEtaPrime2MuMu,            "isEtaPrime2MuMu/b");
    tree->Branch("isEtaPrime2MuMuGamma",     &isEtaPrime2MuMuGamma,             "isEtaPrime2MuMuGamma/b");
    tree->Branch("isOmega2MuMu",             &isOmega2MuMu,             "isOmega2MuMu/b");
    tree->Branch("isOmega2Pi0MuMu",          &isOmega2Pi0MuMu,            "isOmega2Pi0MuMu/b");
    tree->Branch("isRho2MuMu",               &isRho2MuMu,             "isRho2MuMu/b");
    tree->Branch("isPhi2MuMu",               &isPhi2MuMu,             "isPhi2MuMu/b");
    tree->Branch("isPhi2KK",                 &isPhi2KK,             "isPhi2KK/b");

    //tree->Branch("rho"                 , &rho                         , "rho/F"     );

    tree->Branch("probVtx"            , &probVtx                      , "probVtx/F"  );
    tree->Branch("vtxX"               , &vtxX                         , "vtxX/F"  );
    tree->Branch("vtxY"               , &vtxY                         , "vtxY/F"  );
    tree->Branch("vtxZ"               , &vtxZ                         , "vtxZ/F"  );
    tree->Branch("vtxXError"          , &vtxXError                    , "vtxXError/F"  );
    tree->Branch("vtxYError"          , &vtxYError                    , "vtxYError/F"  );
    tree->Branch("vtxZError"          , &vtxZError                    , "vtxZError/F"  );
    tree->Branch("vtx_chi2"           , &vtx_chi2                     , "vtx_chi2/F"  );

    tree->Branch("npv"                , &npv                          , "npv/I"  );
    tree->Branch("pvX"                , &pvX                          , "pvX/F"  );
    tree->Branch("pvY"                , &pvY                          , "pvY/F"  );
    tree->Branch("pvZ"                , &pvZ                          , "pvZ/F"  );

    tree->Branch("muonID1", "std::vector<bool>", &muonID1, 32000, 0);
    tree->Branch("muonID2", "std::vector<bool>", &muonID2, 32000, 0);

    tree->Branch("l1Result", "std::vector<bool>"             ,&l1Result_, 32000, 0  );
    tree->Branch("hltResult", "std::vector<bool>"             ,&hltResult_, 32000, 0  );

    tree->Branch("nPhotons"               , &nPhotons                       , "nPhotons/I"   );

    tree->Branch("slimmedPhotonPt", "std::vector<float>", &slimmedPhotonPt, 32000, 0);
    tree->Branch("slimmedPhotonEta", "std::vector<float>", &slimmedPhotonEta, 32000, 0);
    tree->Branch("slimmedPhotonPhi", "std::vector<float>", &slimmedPhotonPhi, 32000, 0);
    tree->Branch("slimmedPhotonM", "std::vector<float>", &slimmedPhotonM, 32000, 0);
    tree->Branch("slimmedPhotonSigmaIetaIeta", "std::vector<float>", &slimmedPhotonSigmaIetaIeta, 32000, 0);
    tree->Branch("slimmedPhotonHOverE", "std::vector<float>", &slimmedPhotonHOverE, 32000, 0);
    tree->Branch("slimmedPhotonEcalIso", "std::vector<float>", &slimmedPhotonEcalIso, 32000, 0);
    tree->Branch("slimmedPhotonHcalIso", "std::vector<float>", &slimmedPhotonHcalIso, 32000, 0);
    tree->Branch("slimmedPhotonTrkIso", "std::vector<float>", &slimmedPhotonTrkIso, 32000, 0);
    tree->Branch("slimmedPhotonR9", "std::vector<float>", &slimmedPhotonR9, 32000, 0);

    tree->Branch("pfCandPhotonDr", "std::vector<float>", &pfCandPhotonDr, 32000, 0);
    tree->Branch("pfCandPhotonIso", "std::vector<float>", &pfCandPhotonIso, 32000, 0);
    tree->Branch("pfCandPhotonPt", "std::vector<float>", &pfCandPhotonPt, 32000, 0);
    tree->Branch("pfCandPhotonEta", "std::vector<float>", &pfCandPhotonEta, 32000, 0);
    tree->Branch("pfCandPhotonPhi", "std::vector<float>", &pfCandPhotonPhi, 32000, 0);
    tree->Branch("pfCandPhotonEnergy", "std::vector<float>", &pfCandPhotonEnergy, 32000, 0);
    tree->Branch("pfCandPhotonEt", "std::vector<float>", &pfCandPhotonEt, 32000, 0);
    tree->Branch("pfCandPhotonEt2", "std::vector<float>", &pfCandPhotonEt2, 32000, 0);
}

// ------------ method called once each job just after ending the event loop  ------------
void MuMuGammaTreeMaker::endJob() {
  // please remove this method if not needed
}

void MuMuGammaTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    // HLT paths
    triggerPathsVector.push_back("HLT_DoubleMu4_3_LowMass_v*");
    triggerPathsVector.push_back("HLT_DoubleMu4_LowMass_Displaced_v*");

    HLTConfigProvider hltConfig;
    bool changedConfig = false;
    hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);

    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        triggerPathsMap[triggerPathsVector[i]] = -1;
    }

    for(size_t i = 0; i < triggerPathsVector.size(); i++){
        TPRegexp pattern(triggerPathsVector[i]);
        for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
            std::string pathName = hltConfig.triggerNames()[j];
            if(TString(pathName).Contains(pattern)){
                triggerPathsMap[triggerPathsVector[i]] = j;
            }
        }
    }
}

void MuMuGammaTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void MuMuGammaTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void MuMuGammaTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuMuGammaTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(MuMuGammaTreeMaker);
