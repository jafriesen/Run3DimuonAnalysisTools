#include "MMGUtils.hh"

//helper function to calculate the vertices just given the transient tracks list
// returns probability of the vertex fit (-1 if the fit failed).
float calcVertices(vector<reco::TransientTrack> transient_tracks, TransientVertex tv, std::string type) {
    //float vxy = -9999;
    //float sigma_vxy = -9999;
    //float vtx_chi2;
    //float vz = -9999;
    float prob = -999.;

    //only process valid transient vertices
    if (tv.isValid()) {
        reco::Vertex vertex = reco::Vertex(tv);
        //vxy = sqrt(vertex.x()*vertex.x() + vertex.y()*vertex.y());
        //sigma_vxy = (1/vxy)*sqrt(vertex.x()*vertex.x()*vertex.xError()*vertex.xError() +
        //        vertex.y()*vertex.y()*vertex.yError()*vertex.yError());
        //sigma_vxy = (1/vxy)*(vertex.x()*vertex.xError() + vertex.y()*vertex.yError());
        //float vtx_chi2 = vertex.normalizedChi2();
        //vz = vertex.z();

        //get probability for this vertex
        prob = TMath::Prob( vertex.chi2() , vertex.ndof() );
    }
    return prob;  
} 

TransientVertex computeVertex(pat::Muon & coll_1, pat::Muon & coll_2, std::string type, edm::ESHandle<TransientTrackBuilder> theB, KalmanVertexFitter kvf) {
  reco::Track part_1, part_2;
  part_1 = *(coll_1.bestTrack());
  part_2 = *(coll_2.bestTrack());
  //first build the transient vertex and transient tracks.
  //float dr = -9999;
  //float zv; 
  vector<reco::TransientTrack> transient_tracks{};
  transient_tracks.push_back(theB->build(fix_track(&part_1)));
  transient_tracks.push_back(theB->build(fix_track(&part_2)));
  TransientVertex tv = kvf.vertex(transient_tracks);
  //float probVtx = calcVertices(transient_tracks, tv, type);
  //if ( probVtx > 0.1 ) {
  //  dr = reco::deltaR(part_1, part_2);
  //  nt.mumuVtxDr_.push_back(dr);
  //}
  return tv;
}


//got this code from Sergey Polikarpov
reco::Track fix_track(const reco::TrackRef& tk)
{
    reco::Track t = reco::Track(*tk);
    return fix_track(&t);
}

/* Check for a not positive definite covariance matrix. If the covariance matrix is not positive definite, we force it to be positive definite by
 * adding the minimum eigenvalue to the diagonal of the covariance matrix plus `delta`.
 * See https://nhigham.com/2020/12/22/what-is-a-modified-cholesky-factorization/ */
reco::Track fix_track(const reco::Track *tk, double delta)
{
    unsigned int i, j;
    double min_eig = 1;

    /* Get the original covariance matrix. */
    reco::TrackBase::CovarianceMatrix cov = tk->covariance();

    /* Convert it from an SMatrix to a TMatrixD so we can get the eigenvalues. */
    TMatrixDSym new_cov(cov.kRows);
    for (i = 0; i < cov.kRows; i++) {
        for (j = 0; j < cov.kRows; j++) {
            /* Need to check for nan or inf, because for some reason these
             * cause a segfault when calling Eigenvectors().
             *
             * No idea what to do here or why this happens. */
            if (std::isnan(cov(i,j)) || std::isinf(cov(i,j)))
                cov(i,j) = 1e-6;
            new_cov(i,j) = cov(i,j);
        }
    }

    /* Get the eigenvalues. */
    TVectorD eig(cov.kRows);
    new_cov.EigenVectors(eig);
    for (i = 0; i < cov.kRows; i++)
        if (eig(i) < min_eig)
            min_eig = eig(i);

    /* If the minimum eigenvalue is less than zero, then subtract it from the
     * diagonal and add `delta`. */
    if (min_eig < 0) {
        for (i = 0; i < cov.kRows; i++)
            cov(i,i) -= min_eig - delta;
    }

  return reco::Track(tk->chi2(), tk->ndof(), tk->referencePoint(), tk->momentum(), tk->charge(), cov, tk->algo(), (reco::TrackBase::TrackQuality) tk->qualityMask());
}


double photonPfIso03(pat::PackedCandidate pho, edm::Handle<pat::PackedCandidateCollection> pfcands) {

    double ptsum=0.0;

    for (const pat::PackedCandidate &pfc : *pfcands) {

        double dr = deltaR(pho.p4(), pfc.p4());

        if (dr>=0.3) continue;

        if (pfc.charge()!=0 && abs(pfc.pdgId())==211 && pfc.pt()>0.2) {
            if (dr>0.0001) ptsum+=pfc.pt();
        } 
        else if (pfc.charge()==0 && (abs(pfc.pdgId())==22||abs(pfc.pdgId())==130) && pfc.pt()>0.5) {
            if (dr>0.01) ptsum+=pfc.pt();
        }

    }
    return ptsum;
}

float validMuonHitComb(const pat::Muon& muon) {

    const reco::HitPattern& gMpattern = muon.globalTrack()->hitPattern();

    std::vector<int> fvDThits{0, 0, 0, 0};
    std::vector<int> fvRPChits{0, 0, 0, 0};
    std::vector<int> fvCSChits{0, 0, 0, 0};

    float vMuonHitComb = 0;

    for (int i = 0; i < gMpattern.numberOfAllHits(reco::HitPattern::TRACK_HITS); i++) {
    uint32_t hit = gMpattern.getHitPattern(reco::HitPattern::TRACK_HITS, i);
    if (!gMpattern.validHitFilter(hit))
      continue;

    int muStation0 = gMpattern.getMuonStation(hit) - 1;
    if (muStation0 >= 0 && muStation0 < 4) {
      if (gMpattern.muonDTHitFilter(hit))
        fvDThits[muStation0]++;
      if (gMpattern.muonRPCHitFilter(hit))
        fvRPChits[muStation0]++;
      if (gMpattern.muonCSCHitFilter(hit))
        fvCSChits[muStation0]++;
    }
    }

    for (unsigned int station = 0; station < 4; ++station) {
        vMuonHitComb += (fvDThits[station]) / 2.;
        vMuonHitComb += fvRPChits[station];

        if (fvCSChits[station] > 6) {
          vMuonHitComb += 6;
        } else {
          vMuonHitComb += fvCSChits[station];
        }
    }

    return vMuonHitComb;
}
