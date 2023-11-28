#ifndef __UTILS__
#define __UTILS__

#include <vector>
using std::vector;
#include <TLorentzVector.h>

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TVectorD.h"    // for fixing tracks
#include "TMatrixDSym.h" // for fixing tracks
#include "DataFormats/Math/interface/deltaR.h"

// PAT
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

const float ele_mass = 0.000511; //GeV
const float mu_mass = 0.105658; //GeV
const float pi_mass = 0.140; //GeV
const float MIN_DR_TRUTH = 0.05;
TransientVertex computeVertex(pat::Muon & coll_1, pat::Muon & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf);

//struct for saving the used tracks and muons
// -- used as return type for computeVertices functions
struct VertexTracks {
    vector<reco::Track> tracksP;
    vector<reco::Track> tracksN;
    vector<pat::Muon> muonsP;
    vector<pat::Muon> muonsN;
};

//fix tracks with non-pos-def covariance matrices -- needed to prevent crashing
reco::Track fix_track(const reco::Track *tk, double delta=1e-8);

//got this code from Sergey Polikarpov
reco::Track fix_track(const reco::TrackRef& tk);

double photonPfIso03(pat::PackedCandidate pho, edm::Handle<pat::PackedCandidateCollection> pfcands);

#endif // __UTILS__
