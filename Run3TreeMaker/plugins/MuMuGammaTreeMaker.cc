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
  bool isAncestor(const reco::GenParticle* ancestor, const reco::Candidate* particle);
  bool isMatched(const reco::Candidate* gen_particle, const TLorentzVector* reco_vector, float cand_mass);
  bool isPi0(const std::vector<float>& photonsPt, const std::vector<float>& photonsEta, const std::vector<float>& photonsPhi);

private:
  static constexpr float kMinMuonPt = 3.0f;
  static constexpr float kMinVtxProb = 0.005f; // loose minimum probability as in dimuon HLT path
  static constexpr float kMinPfCandPhotonPt = 1.0f;
  static constexpr float kMaxPfCandDimuDeltaR = 0.5f;
  enum class PDGID : int {
    ID_PHOTON = 22,
    ID_ETA = 221,
    ID_ETA_PRIME = 331,
    ID_OMEGA = 223,
    ID_RHO = 113,
    ID_PHI = 333
  };

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>             triggerResultsToken;
  const edm::EDGetTokenT<std::vector<pat::Muon> >      muonsToken;
  const edm::EDGetTokenT<std::vector<reco::Vertex> >    primaryVerticesToken;
  const edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate> >    verticesToken;
  const edm::EDGetTokenT<std::vector<pat::Photon> >         photonsToken;
  const edm::EDGetTokenT<std::vector<pat::PackedCandidate> >         pfCandsToken;
  const edm::EDGetTokenT<std::vector<reco::GenParticle> >         prunedGenToken;
  const edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >    packedGenToken;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> esToken;

  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;

  bool doL1;
  bool doGEN;
  bool doFullGEN;
  triggerExpression::Data triggerCache_;

  edm::InputTag                algInputTag_;
  edm::InputTag                extInputTag_;
  edm::EDGetToken              algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;
  std::vector<bool>            l1Result_;
  std::vector<bool>            hltResult_;

  TTree* tree;

  // event info
  int eventNum;
  int lumiSec;
  int runNum;

  // dimuon variables
  float mumu_mass;
  float mumu_pt;
  float mumu_deltaR;

  // primary vertex variables
  reco::Vertex pv;
  int npv;
  float pv_pos[3];

  // fitted vertex variables
  reco::Vertex vertex;
  float vtx_prob;
  float vtx_pos[3];
  float vtx_posError[3];
  float vtx_chi2;

  // muon variables
  TLorentzVector mu_v[2];
  std::vector<bool> mu_ID[2];
  int   mu_idx[2]; // let standard be 0 = mu-, 1 = mu+
  float mu_pfIso[2];
  float mu_pt[2];
  float mu_eta[2];
  float mu_phi[2];
  float mu_dxy[2];
  float mu_dz[2];
  float mu_trkChi2[2];
  float mu_trkNdof[2];
  float mu_pdgId[2];
  // muon MVA variables
  float mu_segmentCompatibility[2];
  float mu_chi2LocalMomentum[2];
  float mu_chi2LocalPosition[2];
  float mu_glbTrackProbability[2];
  float mu_iValidFraction[2];
  float mu_layersWithMeasurement[2];
  float mu_trkKink[2];
  float mu_log2PlusGlbKink[2];
  float mu_timeAtIpInOutErr[2];
  float mu_outerChi2[2];
  float mu_innerChi2[2];
  float mu_trkRelChi2[2];
  float mu_vMuonHitComb[2];
  float mu_qProd[2];
  float mu_mva[2];

  // photon variables
  // slimmed photons
  std::vector<float> slimmedPhoton_pt;
  std::vector<float> slimmedPhoton_eta;
  std::vector<float> slimmedPhoton_phi;
  std::vector<float> slimmedPhoton_mass;
  std::vector<float> slimmedPhoton_sigmaIetaIeta;
  std::vector<float> slimmedPhoton_hOverE;
  std::vector<float> slimmedPhoton_ecalIso;
  std::vector<float> slimmedPhoton_hcalIso;
  std::vector<float> slimmedPhoton_trkIso;
  std::vector<float> slimmedPhoton_r9;
  // pfCand photons
  std::vector<float> pfCandPhoton_deltaR;
  std::vector<float> pfCandPhoton_iso;
  std::vector<float> pfCandPhoton_pt;
  std::vector<float> pfCandPhoton_eta;
  std::vector<float> pfCandPhoton_phi;
  std::vector<float> pfCandPhoton_energy;
  std::vector<float> pfCandPhoton_et;
  std::vector<float> pfCandPhoton_et2;

  /*
  // doGEN variables
  int motherID[2];
  int motherGenID;
  int genID[2];
  int simType[2];
  int simExtType[2];
  std::vector<int> matchedDaughtersIDs;
  std::vector<float> matchedPhotonPt;
  std::vector<float> matchedPhotonEta;
  std::vector<float> matchedPhotonPhi;
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
  */

  int motherID[2];
  int simType[2];
  int simExtType[2];
  int matchedGenIdx;

  std::vector< std::vector<int> > gen_matchedDaughtersIDs;
  std::vector< std::vector<float> > gen_matchedPhotonPt;
  std::vector< std::vector<float> > gen_matchedPhotonEta;
  std::vector< std::vector<float> > gen_matchedPhotonPhi;

  std::vector<int> gen_motherID;
  std::vector<float> gen_mumu_mass;
  std::vector<float> gen_mu_pt[2];
  std::vector<float> gen_mu_eta[2];
  std::vector<float> gen_mu_phi[2];
  std::vector<bool> gen_mu_recoMatch[2];

  std::vector<bool> gen_isEta2MuMu;
  std::vector<bool> gen_isEta2MuMuGamma;
  std::vector<bool> gen_isEtaPrime2MuMu;
  std::vector<bool> gen_isEtaPrime2MuMuGamma;
  std::vector<bool> gen_isOmega2MuMu;
  std::vector<bool> gen_isOmega2Pi0MuMu;
  std::vector<bool> gen_isRho2MuMu;
  // std::vector<bool> isRho2Pi0MuMu;
  std::vector<bool> gen_isPhi2MuMu;
  std::vector<bool> gen_isPhi2KK;
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
    primaryVerticesToken     (consumes<std::vector<reco::Vertex> >           (iConfig.getParameter<edm::InputTag>("primaryVertices"))),
    verticesToken            (consumes<std::vector<reco::VertexCompositePtrCandidate> >           (iConfig.getParameter<edm::InputTag>("displacedVertices"))),
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

//Check recursively if any ancestor of particle is the given one
bool MuMuGammaTreeMaker::isAncestor(const reco::GenParticle* ancestor, const reco::Candidate* particle) {
  // cast to the base class to make direct comparison
  const reco::Candidate* ancestorPtr = ancestor;
  //particle is already the ancestor
          if(ancestorPtr == particle ) return true;

  //otherwise loop on mothers, if any and return true if the ancestor is found
          for(size_t i=0;i< particle->numberOfMothers();i++)
          {
            const reco::Candidate* motherPtr = particle->mother(i);
            if(isAncestor(ancestor, motherPtr)) return true;
          }
  //if we did not return yet, then particle and ancestor are not relatives
          return false;
}

//check if invariant mass of 2 photons is close to pi0
bool MuMuGammaTreeMaker::isPi0(const std::vector<float>& photonsPt, const std::vector<float>& photonsEta, const std::vector<float>& photonsPhi) {
  TLorentzVector photon1;
  TLorentzVector photon2;
  TLorentzVector diPhoton;
  float diPhotonMass;

  photon1.SetPtEtaPhiM( photonsPt[0], photonsEta[0], photonsPhi[0], 0.);
  photon2.SetPtEtaPhiM( photonsPt[1], photonsEta[1], photonsPhi[1], 0.);
  diPhoton = photon1 + photon2;
  diPhotonMass = diPhoton.M();

  return ((diPhotonMass >= (PI0_MASS - PI0_MASS_SHIFT)) and (diPhotonMass <= (PI0_MASS + PI0_MASS_SHIFT)));

}

// check if a vector of gen particle and reco vector match (by dR)
bool MuMuGammaTreeMaker::isMatched(const reco::Candidate* gen_particle, const TLorentzVector* reco_vector, float cand_mass) {
  bool is_matched = false;
  TLorentzVector gen_vec;
  gen_vec.SetPtEtaPhiM( gen_particle->pt(), gen_particle->eta(), gen_particle->phi(), cand_mass);
  float dr_gen_reco = gen_vec.DeltaR(*reco_vector);
  is_matched = dr_gen_reco <= MIN_DR_TRUTH;

  return is_matched;
}


// ------------ method called for each event  ------------
void MuMuGammaTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;

  // RECO
  bool fillTree = false;
  Handle<vector<pat::Muon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);
  if ( muonsH->size() > 1 ) { // check for 2+ muons
    Handle<vector<reco::Vertex> > primaryVerticesH;
    iEvent.getByToken(primaryVerticesToken, primaryVerticesH);
    npv = primaryVerticesH->size();
    if ( npv > 0 ) { // check for primary vertex
      pv = *primaryVerticesH->begin();

      float bestProbVtx = 0;
      edm::ESHandle<TransientTrackBuilder> theB = iSetup.getHandle(esToken);
      KalmanVertexFitter kvf(true);
      TransientVertex tv;
      for ( unsigned int iMuon = 0; iMuon < muonsH->size(); iMuon++ ) {
        if ( (muonsH->at(iMuon)).pt() > kMinMuonPt ) {
          for ( unsigned int jMuon = iMuon+1; jMuon < muonsH->size(); jMuon++ ) {
            if ( (muonsH->at(jMuon)).pt() > kMinMuonPt && ( (muonsH->at(iMuon)).charge() * (muonsH->at(jMuon)).charge() < 0 ) ) {
              reco::Track part_1 = *((muonsH->at(iMuon)).bestTrack());
              reco::Track part_2 = *((muonsH->at(jMuon)).bestTrack());
              vector<reco::TransientTrack> transient_tracks{};
              transient_tracks.push_back(theB->build(fix_track(&part_1)));
              transient_tracks.push_back(theB->build(fix_track(&part_2)));
              tv = kvf.vertex(transient_tracks);

              if (!tv.isValid()) {
                //std::cout << "ij " << iMuon << jMuon << "Vertex not valid." << std::endl;
              } else {
                reco::Vertex currentVtx = reco::Vertex(tv);
                float currentVtxProb = TMath::Prob( vertex.chi2() , vertex.ndof() );

                if (currentVtxProb > kMinVtxProb && currentVtxProb > bestProbVtx) {
                  vertex = currentVtx;
                  bestProbVtx = currentVtxProb;
                  bool iMuonIsPositive = (muonsH->at(iMuon)).charge() > 0;
                  mu_idx[iMuonIsPositive] = iMuon; // if iMuon is negative/positive, save to mu_idx[0]/mu_idx[1]
                  mu_idx[!iMuonIsPositive] = jMuon; // if iMuon is negative/positive, jMuon is positive/negative, save jMuon to mu_idx[1]/mu_idx[0]
                }
              }
            }
          }
        }
      }

      if (bestProbVtx > kMinVtxProb) { // will fill tree if true
        fillTree = true;
        vtx_prob = bestProbVtx;
        vtx_chi2 = vertex.normalizedChi2();
        vtx_pos[0] = vertex.x();
        vtx_pos[1] = vertex.y();
        vtx_pos[2] = vertex.z();
        vtx_posError[0] = vertex.xError();
        vtx_posError[1] = vertex.yError();
        vtx_posError[2] = vertex.zError();

        pv_pos[0] = pv.x();
        pv_pos[1] = pv.y();
        pv_pos[2] = pv.z();

        eventNum = iEvent.id().event();
        lumiSec = iEvent.luminosityBlock();
        runNum = iEvent.id().run();

        for ( int i : {0, 1} ) {
          auto iMu = muonsH->at(mu_idx[i]);

          mu_pt[i] = iMu.pt();
          mu_eta[i] = iMu.eta();
          mu_phi[i] = iMu.phi();
          mu_v[i].SetPtEtaPhiM(mu_pt[i],mu_eta[i],mu_phi[i],mu_mass);

          // Compute the pfIso for the muon. Note: PUCorr = 0.5*muons_iter->puChargedHadronIso()                                                                              
          // -----------------------------------------------------------------------------------                                                                              
          mu_pfIso[i] = (iMu.chargedHadronIso() + std::max(iMu.photonIso() + iMu.neutralHadronIso() - 0.5 * iMu.puChargedHadronIso() , 0.0)) / iMu.pt();

          mu_dxy[i] = iMu.muonBestTrack()->dxy();
          mu_dz[i] = iMu.muonBestTrack()->dz();
          mu_trkChi2[i] = iMu.muonBestTrack()->chi2();
          mu_trkNdof[i] = iMu.muonBestTrack()->ndof();

          mu_pdgId[i] = iMu.pdgId();

          mu_ID[i].clear();
          mu_ID[i].push_back(iMu.isHighPtMuon(pv));
          mu_ID[i].push_back(iMu.isLooseMuon());
          mu_ID[i].push_back(iMu.isMediumMuon());
          mu_ID[i].push_back(iMu.isSoftMuon(pv));
          mu_ID[i].push_back(iMu.isTightMuon(pv));

          // MVA variables
          mu_chi2LocalMomentum[i] = iMu.combinedQuality().chi2LocalMomentum;
          mu_chi2LocalPosition[i] = iMu.combinedQuality().chi2LocalPosition;
          mu_glbTrackProbability[i] = iMu.combinedQuality().glbTrackProbability;
          mu_trkKink[i] = iMu.combinedQuality().trkKink;
          mu_log2PlusGlbKink[i] = TMath::Log(2 + iMu.combinedQuality().glbKink);
          mu_trkRelChi2[i] = iMu.combinedQuality().trkRelChi2;
          mu_segmentCompatibility[i] = iMu.segmentCompatibility();
          mu_timeAtIpInOutErr[i] = iMu.time().timeAtIpInOutErr;
          reco::TrackRef gTrack = iMu.globalTrack();
          reco::TrackRef iTrack = iMu.innerTrack();
          reco::TrackRef oTrack = iMu.outerTrack();
          if (!(iTrack.isNonnull() and oTrack.isNonnull() and gTrack.isNonnull())) {
            std::cout << "event " << eventNum << ": null track" << std::endl;
            mu_iValidFraction[i] = -1;
            mu_innerChi2[i] = -1;
            mu_layersWithMeasurement[i] = -1;
            mu_outerChi2[i] = -1;
            mu_qProd[i] = -1;
            mu_vMuonHitComb[i] = -1;
            mu_mva[i] = -1;
          } else {
            mu_iValidFraction[i] = iTrack->validFraction();
            mu_innerChi2[i] = iTrack->normalizedChi2();
            mu_layersWithMeasurement[i] = iTrack->hitPattern().trackerLayersWithMeasurement();
            mu_outerChi2[i] = oTrack->normalizedChi2();
            mu_qProd[i] = iTrack->charge() * oTrack->charge();
            mu_vMuonHitComb[i] = validMuonHitComb( iMu );
            mu_mva[i] = iMu.softMvaValue();
          }
        }

        TLorentzVector dimu = mu_v[0] + mu_v[1];
        mumu_mass = dimu.M();
        mumu_pt = dimu.Pt();
        mumu_deltaR = mu_v[0].DeltaR(mu_v[1]);

        Handle<vector<pat::Photon> > photonsH;
        iEvent.getByToken(photonsToken, photonsH);
        slimmedPhoton_pt.clear();
        slimmedPhoton_eta.clear();
        slimmedPhoton_phi.clear();
        slimmedPhoton_mass.clear();
        slimmedPhoton_sigmaIetaIeta.clear();
        slimmedPhoton_hOverE.clear();
        slimmedPhoton_ecalIso.clear();
        slimmedPhoton_hcalIso.clear();
        slimmedPhoton_trkIso.clear();
        slimmedPhoton_r9.clear();
        for (auto photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
          slimmedPhoton_pt.push_back(photons_iter->pt());
          slimmedPhoton_eta.push_back(photons_iter->eta());
          slimmedPhoton_phi.push_back(photons_iter->phi());
          slimmedPhoton_mass.push_back(photons_iter->mass());
          slimmedPhoton_sigmaIetaIeta.push_back(photons_iter->sigmaIetaIeta());
          slimmedPhoton_hOverE.push_back(photons_iter->hadronicOverEm());
          slimmedPhoton_ecalIso.push_back(photons_iter->ecalIso());
          slimmedPhoton_hcalIso.push_back(photons_iter->hcalIso());
          slimmedPhoton_trkIso.push_back(photons_iter->trackIso());
          slimmedPhoton_r9.push_back(photons_iter->r9());
        }

        Handle<vector<pat::PackedCandidate> > pfCandH;
        iEvent.getByToken(pfCandsToken, pfCandH);
        pfCandPhoton_deltaR.clear();
        pfCandPhoton_iso.clear();
        pfCandPhoton_pt.clear();
        pfCandPhoton_eta.clear();
        pfCandPhoton_phi.clear();
        pfCandPhoton_energy.clear();
        pfCandPhoton_et.clear();
        pfCandPhoton_et2.clear();
        for (auto pfCand_iter = pfCandH->begin(); pfCand_iter != pfCandH->end(); ++pfCand_iter) {
          if (abs(pfCand_iter->pdgId()) != 22 or pfCand_iter->pt() < kMinPfCandPhotonPt) continue;
          float pfCandDimuDr = deltaR(dimu.Eta(), dimu.Phi(), pfCand_iter->eta(), pfCand_iter->phi());
          if (pfCandDimuDr < kMaxPfCandDimuDeltaR) {
            pfCandPhoton_deltaR.push_back(pfCandDimuDr);
            pfCandPhoton_iso.push_back(photonPfIso03(*pfCand_iter,pfCandH)/pfCand_iter->pt());
            pfCandPhoton_pt.push_back(pfCand_iter->pt());
            pfCandPhoton_eta.push_back(pfCand_iter->eta());
            pfCandPhoton_phi.push_back(pfCand_iter->phi());
            pfCandPhoton_energy.push_back(pfCand_iter->energy());
            pfCandPhoton_et.push_back(pfCand_iter->et());
            pfCandPhoton_et2.push_back(pfCand_iter->et2());
          }
        }

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
      }
    }
  }

  if ( doFullGEN or (doGEN and fillTree) ) {
    if ( fillTree ) {
      motherID[0]=muonsH->at(mu_idx[0]).simMotherPdgId();
      motherID[1]=muonsH->at(mu_idx[1]).simMotherPdgId();
      simType[0]=muonsH->at(mu_idx[0]).simType();
      simType[1]=muonsH->at(mu_idx[1]).simType();
      simExtType[0]=muonsH->at(mu_idx[0]).simExtType();
      simExtType[1]=muonsH->at(mu_idx[1]).simExtType();
    }

    gen_motherID.clear();
    matchedGenIdx = -1;
    gen_reco_match[0].clear();
    gen_reco_match[1].clear();

    gen_mumu_mass.clear();
    gen_mu_pt[0].clear();
    gen_mu_pt[1].clear();
    gen_mu_eta[0].clear();
    gen_mu_eta[1].clear();
    gen_mu_phi[0].clear();
    gen_mu_phi[1].clear();

    gen_isEta2MuMu.clear();
    gen_isEta2MuMuGamma.clear();
    gen_isEtaPrime2MuMu.clear();
    gen_isEtaPrime2MuMuGamma.clear();
    gen_isOmega2MuMu.clear();
    gen_isOmega2Pi0MuMu.clear();
    gen_isRho2MuMu.clear();
    gen_isPhi2MuMu.clear();
    gen_isPhi2KK.clear();

    bool isEta2MuMu            = false;
    bool isEta2MuMuGamma       = false;
    bool isEtaPrime2MuMu       = false;
    bool isEtaPrime2MuMuGamma  = false;
    bool isOmega2MuMu          = false;
    bool isOmega2Pi0MuMu       = false;
    bool isRho2MuMu            = false;
    bool isPhi2MuMu            = false;
    bool isPhi2KK              = false;

    gen_matchedDaughtersIDs.clear();
    gen_matchedPhotonPt.clear();
    gen_matchedPhotonEta.clear();
    gen_matchedPhotonPhi.clear();

    std::vector<int> matchedDaughterIDs;
    std::vector<float> matchedPhotonPt;
    std::vector<float> matchedPhotonEta;
    std::vector<float> matchedPhotonPhi;

    Handle<vector<reco::GenParticle> > prunedGenParticles;
    iEvent.getByToken(prunedGenToken, prunedGenParticles);

    Handle<vector<pat::PackedGenParticle> > packedGenParticles;
    iEvent.getByToken(packedGenToken, packedGenParticles);

    int test = 0;
    for (auto genp = prunedGenParticles->begin(); genp != prunedGenParticles->end(); ++genp) {            
      if (abs(genp->pdgId())==221 or abs(genp->pdgId())==113 or abs(genp->pdgId())==223 or abs(genp->pdgId())==331 or abs(genp->pdgId()==333)) {
      // std::cout<<"genp id: "<<genp->pdgId()<<" pt: "<<genp->pt()<<" eta: "<<genp->eta()<<" status: "<<genp->status();
        int nDaughterMuons = 0;
        int mup_idx = 99;
        int mum_idx = 99;

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
          matchedDaughtersIDs.clear();
          test++;
          std::cout<<"test"<<test<<std::endl;
          // try to match pair of reco muons with a pair of daughters
          auto daught_mup =  &(*genp->daughter(mup_idx));
          auto daught_mum =  &(*genp->daughter(mum_idx));
          
          // std::cout<<"Mother with 2 muons ID == " << genp->pdgId() <<std::endl;
          for (int i=0; i<(int)genp->numberOfDaughters(); ++i) {
            matchedDaughtersIDs.push_back(genp->daughter(i)->pdgId());
            // std::cout<<"Daughter "<< i << "  ID == " << genp->daughter(i)->pdgId()<<std::endl;
          }
          gen_matchedDaughtersIDs.push_back(matchedDaughtersIDs);

          bool isMatched1 = false;
          bool isMatched2 = false;
          if (fillTree) {
            isMatched1 = isMatched(daught_mup, &mu_v[1], mu_mass);
            isMatched2 = isMatched(daught_mum, &mu_v[0], mu_mass);
          }

          if ( doFullGEN or ( isMatched1 and isMatched2 ) ){

            if (isMatched1 and isMatched2) matchedGenIdx = gen_motherID.size();

            gen_reco_match[0].push_back(isMatched1);
            gen_reco_match[1].push_back(isMatched2);

            gen_motherID.push_back(genp->pdgId());

            std::cout << (daught_mup + daught_mum).M() << std::endl;
            gen_mumu_mass.push_back((daught_mup + daught_mum).M());
            gen_mu_pt[0].push_back(daught_mum->pt());
            gen_mu_pt[1].push_back(daught_mup->pt());
            gen_mu_eta[0].push_back(daught_mum->eta());
            gen_mu_eta[1].push_back(daught_mup->eta());
            gen_mu_phi[0].push_back(daught_mum->phi());
            gen_mu_phi[1].push_back(daught_mup->phi());

            // std::cout<<"Matched "<< daught_mup->pt() << "  " << mup_reco->Pt() << "  " << daught_mum->pt() << "  " << mum_reco->Pt() << "for mother ID == " << motherGenID <<std::endl;
            //check for photons
            int nPhotons = 0;
            matchedPhotonPt.clear();
            matchedPhotonEta.clear();
            matchedPhotonPhi.clear();
            // look for photons
            for (auto pgp = packedGenParticles->begin(); pgp != packedGenParticles->end(); ++pgp) {
              // derefference and cast to the base class to make direct comparison
              const reco::Candidate* pgpPtr = &(*pgp);
              if (pgp->pdgId()==22 and isAncestor( &(*genp) , pgpPtr)) {
                  std::cout<<"packed photon "<<pgp->pt()<<" mother pt: "<<pgp->motherRef()->pt()<<std::endl;
                  matchedPhotonPt.push_back( pgp->pt() );
                  matchedPhotonEta.push_back( pgp->eta() );
                  matchedPhotonPhi.push_back( pgp->phi() );
                  nPhotons++;   
              }
            }
            gen_matchedPhotonPt.push_back(matchedPhotonPt);
            gen_matchedPhotonEta.push_back(matchedPhotonEta);
            gen_matchedPhotonPhi.push_back(matchedPhotonPhi);
            bool aPi0 = false;
            // check if 2 photons make a pi0
            if (nPhotons == 2){
              aPi0 = isPi0( matchedPhotonPt, matchedPhotonEta, matchedPhotonPhi);
              // if (aPi0) std::cout<<"Matched pi0 to two photons!" <<std::endl;
            }
            // isPhi2KK
            if      ((abs(genp->pdgId()) == 221) and (nPhotons == 0)) isEta2MuMu = true;
            else if ((abs(genp->pdgId()) == 221) and (nPhotons == 1))    isEta2MuMuGamma = true;
            else if ((abs(genp->pdgId()) == 331) and (nPhotons == 0)) isEtaPrime2MuMu = true;
            else if ((abs(genp->pdgId()) == 331) and (nPhotons == 1))    isEtaPrime2MuMuGamma = true;
            else if ((abs(genp->pdgId()) == 223) and (nPhotons == 0)) isOmega2MuMu = true;
            else if ((abs(genp->pdgId()) == 223) and aPi0)    isOmega2Pi0MuMu = true;      
            else if ((abs(genp->pdgId()) == 113) and (nPhotons == 0)) isRho2MuMu = true;
            else if ((abs(genp->pdgId()) == 333) and (nPhotons == 0)) isPhi2MuMu = true;

            gen_isEta2MuMu.push_back(isEta2MuMu);
            gen_isEta2MuMuGamma.push_back(isEta2MuMuGamma);
            gen_isEtaPrime2MuMu.push_back(isEtaPrime2MuMu);
            gen_isEtaPrime2MuMuGamma.push_back(isEtaPrime2MuMuGamma);
            gen_isOmega2MuMu.push_back(isOmega2MuMu);
            gen_isOmega2Pi0MuMu.push_back(isOmega2Pi0MuMu);
            gen_isRho2MuMu.push_back(isRho2MuMu);
            gen_isPhi2MuMu.push_back(isPhi2MuMu);
            gen_isPhi2KK.push_back(isPhi2KK);
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

    tree->Branch("mumu_mass"                , &mumu_mass                        , "mumu_mass/F"    );
    tree->Branch("mumu_pt"                  , &mumu_pt                          , "mumu_pt/F"      );
    tree->Branch("mumu_deltaR"              , &mumu_deltaR                      , "mumu_deltaR/F"      );

    tree->Branch("mu1_pt"                 , &mu_pt[0]                         , "mu1_pt/F"     );
    tree->Branch("mu1_eta"                , &mu_eta[0]                        , "mu1_eta/F"    );
    tree->Branch("mu1_phi"                , &mu_phi[0]                        , "mu1_phi/F"    );
    tree->Branch("mu1_pfIso"              , &mu_pfIso[0]                      , "mu1_pfIso/F"  );

    tree->Branch("mu2_pt"                 , &mu_pt[1]                         , "mu2_pt/F"     );
    tree->Branch("mu2_eta"                , &mu_eta[1]                        , "mu2_eta/F"    );
    tree->Branch("mu2_phi"                , &mu_phi[1]                        , "mu2_phi/F"    );
    tree->Branch("mu2_pfIso"              , &mu_pfIso[1]                      , "mu2_pfIso/F"  );

    tree->Branch("mu1_segmentCompatibility",     &mu_segmentCompatibility[0],     "mu1_segmentCompatibility/F"   );
    tree->Branch("mu1_chi2LocalMomentum",        &mu_chi2LocalMomentum[0],        "mu1_chi2LocalMomentum/F"      );
    tree->Branch("mu1_chi2LocalPosition",        &mu_chi2LocalPosition[0],        "mu1_chi2LocalPosition/F"      );
    tree->Branch("mu1_glbTrackProbability",      &mu_glbTrackProbability[0],      "mu1_glbTrackProbability/F"    );
    tree->Branch("mu1_iValidFraction",           &mu_iValidFraction[0],           "mu1_iValidFraction/F"         );
    tree->Branch("mu1_layersWithMeasurement",    &mu_layersWithMeasurement[0],    "mu1_layersWithMeasurement/F"  );
    tree->Branch("mu1_trkKink",                  &mu_trkKink[0],                  "mu1_trkKink/F"                );
    tree->Branch("mu1_log2PlusGlbKink",          &mu_log2PlusGlbKink[0],          "mu1_log2PlusGlbKink/F"        );
    tree->Branch("mu1_timeAtIpInOutErr",         &mu_timeAtIpInOutErr[0],         "mu1_timeAtIpInOutErr/F"       );
    tree->Branch("mu1_outerChi2",                &mu_outerChi2[0],                "mu1_outerChi2/F"              );
    tree->Branch("mu1_innerChi2",                &mu_innerChi2[0],                "mu1_innerChi2/F"              );
    tree->Branch("mu1_trkRelChi2",               &mu_trkRelChi2[0],               "mu1_trkRelChi2/F"             );
    tree->Branch("mu1_vMuonHitComb",             &mu_vMuonHitComb[0],             "mu1_vMuonHitComb/F"           );
    tree->Branch("mu1_qProd",                    &mu_qProd[0],                    "mu1_qProd/F"                  );
    tree->Branch("mu1_mva",                      &mu_mva[0],                      "mu1_mva/F"                    );

    tree->Branch("mu2_segmentCompatibility",     &mu_segmentCompatibility[1],     "mu2_segmentCompatibility/F"   );
    tree->Branch("mu2_chi2LocalMomentum",        &mu_chi2LocalMomentum[1],        "mu2_chi2LocalMomentum/F"      );
    tree->Branch("mu2_chi2LocalPosition",        &mu_chi2LocalPosition[1],        "mu2_chi2LocalPosition/F"      );
    tree->Branch("mu2_glbTrackProbability",      &mu_glbTrackProbability[1],      "mu2_glbTrackProbability/F"    );
    tree->Branch("mu2_iValidFraction",           &mu_iValidFraction[1],           "mu2_iValidFraction/F"         );
    tree->Branch("mu2_layersWithMeasurement",    &mu_layersWithMeasurement[1],    "mu2_layersWithMeasurement/F"  );
    tree->Branch("mu2_trkKink",                  &mu_trkKink[1],                  "mu2_trkKink/F"                );
    tree->Branch("mu2_log2PlusGlbKink",          &mu_log2PlusGlbKink[1],          "mu2_log2PlusGlbKink/F"        );
    tree->Branch("mu2_timeAtIpInOutErr",         &mu_timeAtIpInOutErr[1],         "mu2_timeAtIpInOutErr/F"       );
    tree->Branch("mu2_outerChi2",                &mu_outerChi2[1],                "mu2_outerChi2/F"              );
    tree->Branch("mu2_innerChi2",                &mu_innerChi2[1],                "mu2_innerChi2/F"              );
    tree->Branch("mu2_trkRelChi2",               &mu_trkRelChi2[1],               "mu2_trkRelChi2/F"             );
    tree->Branch("mu2_vMuonHitComb",             &mu_vMuonHitComb[1],             "mu2_vMuonHitComb/F"           );
    tree->Branch("mu2_qProd",                    &mu_qProd[1],                    "mu2_qProd/F"                  );
    tree->Branch("mu2_mva",                      &mu_mva[1],                      "mu2_mva/F"                    );

    tree->Branch("mu1_dxy"              , &mu_dxy[0]                      , "mu1_dxy/F"  );
    tree->Branch("mu1_dz"               , &mu_dz[0]                       , "mu1_dz/F"  );
    tree->Branch("mu1_trkChi2"          , &mu_trkChi2[0]                  , "mu1_trkChi2/F"  );
    tree->Branch("mu1_trkNdof"          , &mu_trkNdof[0]                  , "mu1_trkNdof/F"  );
    
    tree->Branch("mu2_dxy"              , &mu_dxy[1]                      , "mu2_dxy/F"  );
    tree->Branch("mu2_dz"               , &mu_dz[1]                       , "mu2_dz/F"  );
    tree->Branch("mu2_trkChi2"          , &mu_trkChi2[1]                  , "mu2_trkChi2/F"  );
    tree->Branch("mu2_trkNdof"          , &mu_trkNdof[1]                  , "mu2_trkNdof/F"  );

    tree->Branch("motherID1"              , &motherID[0]                     , "motherID1/I"  );
    tree->Branch("motherID2"              , &motherID[1]                      , "motherID2/I" );
    tree->Branch("simType1"              , &simType[0]                      , "simType1/I"  );
    tree->Branch("simType2"              , &simType[1]                      , "simType2/I"  );
    tree->Branch("simExtType1"              , &simExtType[0]                      , "simExtType1/I"  );
    tree->Branch("simExtType2"              , &simExtType[1]                      , "simExtType2/I"  );

    tree->Branch("matchedGenIdx"              , &matchedGenIdx                    , "matchedGenIdx/I"  );

    tree->Branch("gen_motherID", "std::vector<int>", &gen_motherID, 32000, 0);
    tree->Branch("gen_mumu_mass", "std::vector<float>", &gen_mumu_mass, 32000, 0);

    tree->Branch("gen_mu1_pt", "std::vector<float>", &gen_mu_pt[0], 32000, 0);
    tree->Branch("gen_mu1_eta", "std::vector<float>", &gen_mu_eta[0], 32000, 0);
    tree->Branch("gen_mu1_phi", "std::vector<float>", &gen_mu_phi[0], 32000, 0);
    tree->Branch("gen_mu1_recoMatch", "std::vector<bool>", &gen_mu_recoMatch[0], 32000, 0);

    tree->Branch("gen_mu2_pt", "std::vector<float>", &gen_mu_pt[1], 32000, 0);
    tree->Branch("gen_mu2_eta", "std::vector<float>", &gen_mu_eta[1], 32000, 0);
    tree->Branch("gen_mu2_phi", "std::vector<float>", &gen_mu_phi[1], 32000, 0);
    tree->Branch("gen_mu2_recoMatch", "std::vector<bool>", &gen_mu_recoMatch[1], 32000, 0);

    tree->Branch("gen_matchedDaughtersIDs", "std::vector<std::vector<int> >", &gen_matchedDaughtersIDs, 32000, 0);
    tree->Branch("gen_matchedPhotonPt", "std::vector<std::vector<float> >", &gen_matchedPhotonPt, 32000, 0);
    tree->Branch("gen_matchedPhotonEta", "std::vector<std::vector<float> >", &gen_matchedPhotonEta, 32000, 0);
    tree->Branch("gen_matchedPhotonPhi", "std::vector<std::vector<float> >", &gen_matchedPhotonPhi, 32000, 0);

    tree->Branch("gen_isEta2MuMu", "std::vector<bool>", &gen_isEta2MuMu, 32000, 0);
    tree->Branch("gen_isEta2MuMuGamma", "std::vector<bool>", &gen_isEta2MuMuGamma, 32000, 0);
    tree->Branch("gen_isEtaPrime2MuMu", "std::vector<bool>", &gen_isEtaPrime2MuMu, 32000, 0);
    tree->Branch("gen_isEtaPrime2MuMuGamma", "std::vector<bool>", &gen_isEtaPrime2MuMuGamma, 32000, 0);
    tree->Branch("gen_isOmega2MuMu", "std::vector<bool>", &gen_isOmega2MuMu, 32000, 0);
    tree->Branch("gen_isOmega2Pi0MuMu", "std::vector<bool>", &gen_isOmega2Pi0MuMu, 32000, 0);
    tree->Branch("gen_isRho2MuMu", "std::vector<bool>", &gen_isRho2MuMu, 32000, 0);
    tree->Branch("gen_isPhi2MuMu", "std::vector<bool>", &gen_isPhi2MuMu, 32000, 0);
    tree->Branch("gen_isPhi2KK", "std::vector<bool>", &gen_isPhi2KK, 32000, 0);

    //tree->Branch("rho"                 , &rho                         , "rho/F"     );

    tree->Branch("probVtx"            , &vtx_prob                      , "probVtx/F"  );
    tree->Branch("vtxX"               , &vtx_pos[0]                         , "vtxX/F"  );
    tree->Branch("vtxY"               , &vtx_pos[1]                          , "vtxY/F"  );
    tree->Branch("vtxZ"               , &vtx_pos[2]                          , "vtxZ/F"  );
    tree->Branch("vtxXError"          , &vtx_posError[0]                     , "vtxXError/F"  );
    tree->Branch("vtxYError"          , &vtx_posError[1]                    , "vtxYError/F"  );
    tree->Branch("vtxZError"          , &vtx_posError[2]                    , "vtxZError/F"  );
    tree->Branch("vtx_chi2"           , &vtx_chi2                     , "vtx_chi2/F"  );

    tree->Branch("npv"                , &npv                          , "npv/I"  );
    tree->Branch("pvX"                , &pv_pos[0]                          , "pvX/F"  );
    tree->Branch("pvY"                , &pv_pos[1]                          , "pvY/F"  );
    tree->Branch("pvZ"                , &pv_pos[2]                          , "pvZ/F"  );

    tree->Branch("mu1_ID", "std::vector<bool>", &mu_ID[0], 32000, 0);
    tree->Branch("mu2_ID", "std::vector<bool>", &mu_ID[1], 32000, 0);

    tree->Branch("l1Result", "std::vector<bool>"             ,&l1Result_, 32000, 0  );
    tree->Branch("hltResult", "std::vector<bool>"             ,&hltResult_, 32000, 0  );

    tree->Branch("slimmedPhoton_pt", "std::vector<float>", &slimmedPhoton_pt, 32000, 0);
    tree->Branch("slimmedPhoton_eta", "std::vector<float>", &slimmedPhoton_eta, 32000, 0);
    tree->Branch("slimmedPhoton_phi", "std::vector<float>", &slimmedPhoton_phi, 32000, 0);
    tree->Branch("slimmedPhoton_mass", "std::vector<float>", &slimmedPhoton_mass, 32000, 0);
    tree->Branch("slimmedPhoton_sigmaIetaIeta", "std::vector<float>", &slimmedPhoton_sigmaIetaIeta, 32000, 0);
    tree->Branch("slimmedPhoton_hOverE", "std::vector<float>", &slimmedPhoton_hOverE, 32000, 0);
    tree->Branch("slimmedPhoton_ecalIso", "std::vector<float>", &slimmedPhoton_ecalIso, 32000, 0);
    tree->Branch("slimmedPhoton_hcalIso", "std::vector<float>", &slimmedPhoton_hcalIso, 32000, 0);
    tree->Branch("slimmedPhoton_trkIso", "std::vector<float>", &slimmedPhoton_trkIso, 32000, 0);
    tree->Branch("slimmedPhoton_r9", "std::vector<float>", &slimmedPhoton_r9, 32000, 0);

    tree->Branch("pfCandPhoton_deltaR", "std::vector<float>", &pfCandPhoton_deltaR, 32000, 0);
    tree->Branch("pfCandPhoton_iso", "std::vector<float>", &pfCandPhoton_iso, 32000, 0);
    tree->Branch("pfCandPhoton_pt", "std::vector<float>", &pfCandPhoton_pt, 32000, 0);
    tree->Branch("pfCandPhoton_eta", "std::vector<float>", &pfCandPhoton_eta, 32000, 0);
    tree->Branch("pfCandPhoton_phi", "std::vector<float>", &pfCandPhoton_phi, 32000, 0);
    tree->Branch("pfCandPhoton_energy", "std::vector<float>", &pfCandPhoton_energy, 32000, 0);
    tree->Branch("pfCandPhoton_et", "std::vector<float>", &pfCandPhoton_et, 32000, 0);
    tree->Branch("pfCandPhoton_et2", "std::vector<float>", &pfCandPhoton_et2, 32000, 0);
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
