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

    motherGenID = 0; 
    genID[0] = 0; genID[1] = 0;
    isEta2MuMu            = false;
    isEta2MuMuGamma       = false;
    isEtaPrime2MuMu       = false;
    isEtaPrime2MuMuGamma  = false;
    isOmega2MuMu         = false;
    isOmega2Pi0MuMu       = false;
    isRho2MuMu            = false;
    isPhi2MuMu            = false;
    isPhi2KK              = false;

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
          test++;
          std::cout<<"test"<<test<<std::endl;
          // try to match pair of reco muons with a pair of daughters
          auto daught_mup =  &(*genp->daughter(mup_idx));
          auto daught_mum =  &(*genp->daughter(mum_idx));

          matchedDaughtersIDs.clear();
          // std::cout<<"Mother with 2 muons ID == " << genp->pdgId() <<std::endl;
          for (int i=0; i<(int)genp->numberOfDaughters(); ++i) {
            matchedDaughtersIDs.push_back(genp->daughter(i)->pdgId());
            // std::cout<<"Daughter "<< i << "  ID == " << genp->daughter(i)->pdgId()<<std::endl;
          }

          if ( doFullGEN or (fillTree and (isMatched(daught_mup, &mu_v[1], mu_mass) and isMatched(daught_mum, &mu_v[0], mu_mass))) ){
            motherGenID = genp->pdgId();
            genID[1] = daught_mup->pdgId();
            genID[0] = daught_mum->pdgId();
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
                  matchedPhotonPt.push_back( pgp->pt());
                  matchedPhotonEta.push_back( pgp->eta());
                  matchedPhotonPhi.push_back( pgp->phi());
                  nPhotons++;
                  
              }
            }
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

    tree->Branch("mass"                , &mumu_mass                        , "mass/F"    );
    tree->Branch("pt"                  , &mumu_pt                          , "pt/F"      );
    tree->Branch("dr"                  , &mumu_deltaR                          , "dr/F"      );
    tree->Branch("pt1"                 , &mu_pt[0]                         , "pt1/F"     );
    tree->Branch("pt2"                 , &mu_pt[1]                         , "pt2/F"     );
    tree->Branch("eta1"                , &mu_eta[0]                        , "eta1/F"    );
    tree->Branch("eta2"                , &mu_eta[1]                        , "eta2/F"    );
    tree->Branch("phi1"                , &mu_phi[0]                        , "phi1/F"    );
    tree->Branch("phi2"                , &mu_phi[1]                        , "phi2/F"    );
    tree->Branch("pfIso1"              , &mu_pfIso[0]                      , "pfIso1/F"  );
    tree->Branch("pfIso2"              , &mu_pfIso[1]                      , "pfIso2/F"  );

    tree->Branch("segmentCompatibility1",     &mu_segmentCompatibility[0],     "segmentCompatibility1/F"   );
    tree->Branch("chi2LocalMomentum1",        &mu_chi2LocalMomentum[0],        "chi2LocalMomentum1/F"      );
    tree->Branch("chi2LocalPosition1",        &mu_chi2LocalPosition[0],        "chi2LocalPosition1/F"      );
    tree->Branch("glbTrackProbability1",      &mu_glbTrackProbability[0],      "glbTrackProbability1/F"    );
    tree->Branch("iValidFraction1",           &mu_iValidFraction[0],           "iValidFraction1/F"         );
    tree->Branch("layersWithMeasurement1",    &mu_layersWithMeasurement[0],    "layersWithMeasurement1/F"  );
    tree->Branch("trkKink1",                  &mu_trkKink[0],                  "trkKink1/F"                );
    tree->Branch("log2PlusGlbKink1",          &mu_log2PlusGlbKink[0],          "log2PlusGlbKink1/F"        );
    tree->Branch("timeAtIpInOutErr1",         &mu_timeAtIpInOutErr[0],         "timeAtIpInOutErr1/F"       );
    tree->Branch("outerChi21",                &mu_outerChi2[0],                "outerChi21/F"              );
    tree->Branch("innerChi21",                &mu_innerChi2[0],                "innerChi21/F"              );
    tree->Branch("trkRelChi21",               &mu_trkRelChi2[0],               "trkRelChi21/F"             );
    tree->Branch("vMuonHitComb1",             &mu_vMuonHitComb[0],             "vMuonHitComb1/F"           );
    tree->Branch("qProd1",                    &mu_qProd[0],                    "qProd1/F"                  );
    tree->Branch("mva1",                      &mu_mva[0],                      "mva1/F"                    );

    tree->Branch("segmentCompatibility2",     &mu_segmentCompatibility[1],     "segmentCompatibility2/F"   );
    tree->Branch("chi2LocalMomentum2",        &mu_chi2LocalMomentum[1],        "chi2LocalMomentum2/F"      );
    tree->Branch("chi2LocalPosition2",        &mu_chi2LocalPosition[1],        "chi2LocalPosition2/F"      );
    tree->Branch("glbTrackProbability2",      &mu_glbTrackProbability[1],      "glbTrackProbability2/F"    );
    tree->Branch("iValidFraction2",           &mu_iValidFraction[1],           "iValidFraction2/F"         );
    tree->Branch("layersWithMeasurement2",    &mu_layersWithMeasurement[1],    "layersWithMeasurement2/F"  );
    tree->Branch("trkKink2",                  &mu_trkKink[1],                  "trkKink2/F"                );
    tree->Branch("log2PlusGlbKink2",          &mu_log2PlusGlbKink[1],          "log2PlusGlbKink2/F"        );
    tree->Branch("timeAtIpInOutErr2",         &mu_timeAtIpInOutErr[1],         "timeAtIpInOutErr2/F"       );
    tree->Branch("outerChi22",                &mu_outerChi2[1],                "outerChi22/F"              );
    tree->Branch("innerChi22",                &mu_innerChi2[1],                "innerChi22/F"              );
    tree->Branch("trkRelChi22",               &mu_trkRelChi2[1],               "trkRelChi22/F"             );
    tree->Branch("vMuonHitComb2",             &mu_vMuonHitComb[1],             "vMuonHitComb2/F"           );
    tree->Branch("qProd2",                    &mu_qProd[1],                    "qProd2/F"                  );
    tree->Branch("mva2",                      &mu_mva[1],                      "mva2/F"                    );

    tree->Branch("dxy1"              , &mu_dxy[0]                      , "dxy1/F"  );
    tree->Branch("dxy2"              , &mu_dxy[1]                      , "dxy2/F"  );
    tree->Branch("dz1"              , &mu_dz[0]                      , "dz1/F"  );
    tree->Branch("dz2"              , &mu_dz[1]                      , "dz2/F"  );
    tree->Branch("trkChi21"              , &mu_trkChi2[0]                      , "trkChi21/F"  );
    tree->Branch("trkChi22"              , &mu_trkChi2[1]                      , "trkChi22/F"  );
    tree->Branch("trkNdof1"              , &mu_trkNdof[0]                      , "trkNdof1/F"  );
    tree->Branch("trkNdof2"              , &mu_trkNdof[1]                      , "trkNdof2/F"  );

    tree->Branch("motherID1"              , &motherID[0]                     , "motherID1/I"  );
    tree->Branch("motherID2"              , &motherID[1]                      , "motherID2/I"  );
    tree->Branch("simType1"              , &simType[0]                      , "simType1/I"  );
    tree->Branch("simType2"              , &simType[1]                      , "simType2/I"  );
    tree->Branch("simExtType1"              , &simExtType[0]                      , "simExtType1/I"  );
    tree->Branch("simExtType2"              , &simExtType[1]                      , "simExtType2/I"  );

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

    tree->Branch("muonID1", "std::vector<bool>", &mu_ID[0], 32000, 0);
    tree->Branch("muonID2", "std::vector<bool>", &mu_ID[1], 32000, 0);

    tree->Branch("l1Result", "std::vector<bool>"             ,&l1Result_, 32000, 0  );
    tree->Branch("hltResult", "std::vector<bool>"             ,&hltResult_, 32000, 0  );

    tree->Branch("slimmedPhotonPt", "std::vector<float>", &slimmedPhoton_pt, 32000, 0);
    tree->Branch("slimmedPhotonEta", "std::vector<float>", &slimmedPhoton_eta, 32000, 0);
    tree->Branch("slimmedPhotonPhi", "std::vector<float>", &slimmedPhoton_phi, 32000, 0);
    tree->Branch("slimmedPhotonM", "std::vector<float>", &slimmedPhoton_mass, 32000, 0);
    tree->Branch("slimmedPhotonSigmaIetaIeta", "std::vector<float>", &slimmedPhoton_sigmaIetaIeta, 32000, 0);
    tree->Branch("slimmedPhotonHOverE", "std::vector<float>", &slimmedPhoton_hOverE, 32000, 0);
    tree->Branch("slimmedPhotonEcalIso", "std::vector<float>", &slimmedPhoton_ecalIso, 32000, 0);
    tree->Branch("slimmedPhotonHcalIso", "std::vector<float>", &slimmedPhoton_hcalIso, 32000, 0);
    tree->Branch("slimmedPhotonTrkIso", "std::vector<float>", &slimmedPhoton_trkIso, 32000, 0);
    tree->Branch("slimmedPhotonR9", "std::vector<float>", &slimmedPhoton_r9, 32000, 0);

    tree->Branch("pfCandPhotonDr", "std::vector<float>", &pfCandPhoton_deltaR, 32000, 0);
    tree->Branch("pfCandPhotonIso", "std::vector<float>", &pfCandPhoton_iso, 32000, 0);
    tree->Branch("pfCandPhotonPt", "std::vector<float>", &pfCandPhoton_pt, 32000, 0);
    tree->Branch("pfCandPhotonEta", "std::vector<float>", &pfCandPhoton_eta, 32000, 0);
    tree->Branch("pfCandPhotonPhi", "std::vector<float>", &pfCandPhoton_phi, 32000, 0);
    tree->Branch("pfCandPhotonEnergy", "std::vector<float>", &pfCandPhoton_energy, 32000, 0);
    tree->Branch("pfCandPhotonEt", "std::vector<float>", &pfCandPhoton_et, 32000, 0);
    tree->Branch("pfCandPhotonEt2", "std::vector<float>", &pfCandPhoton_et2, 32000, 0);
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
