#ifndef ZHfunctions_H
#define ZHfunctions_H

#include <cmath>
#include <vector>
#include <math.h>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "ReconstructedParticle2MC.h"


namespace FCCAnalyses { namespace ZHfunctions {


// build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
// technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
struct resonanceBuilder_mass_recoil {
    float m_resonance_mass;
    float m_recoil_mass;
    float chi2_recoil_frac;
    float ecm;
    bool m_use_MC_Kinematics;
    resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics);
    Vec_rp operator()(Vec_rp legs, Vec_i recind, Vec_i mcind, Vec_rp reco, Vec_mc mc, Vec_i parents, Vec_i daugthers) ;
};

resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics) {m_resonance_mass = arg_resonance_mass, m_recoil_mass = arg_recoil_mass, chi2_recoil_frac = arg_chi2_recoil_frac, ecm = arg_ecm, m_use_MC_Kinematics = arg_use_MC_Kinematics;}

Vec_rp resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil::operator()(Vec_rp legs, Vec_i recind, Vec_i mcind, Vec_rp reco, Vec_mc mc, Vec_i parents, Vec_i daugthers) {

    Vec_rp result;
    result.reserve(3);
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = legs.size();
  
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            rp reso;
            reso.charge = 0;
            TLorentzVector reso_lv; 
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                    reso.charge += legs[i].charge;
                    TLorentzVector leg_lv;

                    if(m_use_MC_Kinematics) { // MC kinematics
                        int track_index = legs[i].tracks_begin;   // index in the Track array
                        int mc_index = ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
                        if (mc_index >= 0 && mc_index < mc.size()) {
                            leg_lv.SetXYZM(mc.at(mc_index).momentum.x, mc.at(mc_index).momentum.y, mc.at(mc_index).momentum.z, mc.at(mc_index).mass);
                        }
                    }
                    else { // reco kinematics
                         leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
                    }

                    reso_lv += leg_lv;
                }
            }

            if(reso.charge != 0) continue; // neglect non-zero charge pairs
            reso.momentum.x = reso_lv.Px();
            reso.momentum.y = reso_lv.Py();
            reso.momentum.z = reso_lv.Pz();
            reso.mass = reso_lv.M();
            result.emplace_back(reso);
            pairs.push_back(pair);

        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: resonanceBuilder_mass_recoil, at least two leptons required." << std::endl;
        exit(1);
    }
  
    if(result.size() > 1) {
  
        Vec_rp bestReso;
        
        int idx_min = -1;
        float d_min = 9e9;
        for (int i = 0; i < result.size(); ++i) {
            
            // calculate recoil
            auto recoil_p4 = TLorentzVector(0, 0, 0, ecm);
            TLorentzVector tv1;
            tv1.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
            recoil_p4 -= tv1;
      
            auto recoil_fcc = edm4hep::ReconstructedParticleData();
            recoil_fcc.momentum.x = recoil_p4.Px();
            recoil_fcc.momentum.y = recoil_p4.Py();
            recoil_fcc.momentum.z = recoil_p4.Pz();
            recoil_fcc.mass = recoil_p4.M();
            
            TLorentzVector tg;
            tg.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
        
            float boost = tg.P();
            float mass = std::pow(result.at(i).mass - m_resonance_mass, 2); // mass
            float rec = std::pow(recoil_fcc.mass - m_recoil_mass, 2); // recoil
            float d = (1.0-chi2_recoil_frac)*mass + chi2_recoil_frac*rec;
            
            if(d < d_min) {
                d_min = d;
                idx_min = i;
            }

     
        }
        if(idx_min > -1) { 
            bestReso.push_back(result.at(idx_min));
            auto & l1 = legs[pairs[idx_min][0]];
            auto & l2 = legs[pairs[idx_min][1]];
            bestReso.emplace_back(l1);
            bestReso.emplace_back(l2);
        }
        else {
            std::cout << "ERROR: resonanceBuilder_mass_recoil, no mininum found." << std::endl;
            exit(1);
        }
        return bestReso;
    }
    else {
        auto & l1 = legs[0];
        auto & l2 = legs[1];
        result.emplace_back(l1);
        result.emplace_back(l2);
        return result;
    }
}    




struct sel_iso {
    sel_iso(float arg_max_iso);
    float m_max_iso = .25;
    Vec_rp operator() (Vec_rp in, Vec_f iso);
  };

sel_iso::sel_iso(float arg_max_iso) : m_max_iso(arg_max_iso) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sel_iso::operator() (Vec_rp in, Vec_f iso) {
    Vec_rp result;
    result.reserve(in.size());
    for (size_t i = 0; i < in.size(); ++i) {
        auto & p = in[i];
        if (iso[i] < m_max_iso) {
            result.emplace_back(p);
        }
    }
    return result;
}


struct sel_btag {
    sel_btag(float arg_min_btag);
    float m_min_btag = 0.5;
    ROOT::VecOps::RVec<float> operator() (ROOT::VecOps::RVec<float> btag);
};

sel_btag::sel_btag(float arg_min_btag) : m_min_btag(arg_min_btag) {}
ROOT::VecOps::RVec<float> sel_btag::operator() (ROOT::VecOps::RVec<float> btag) {
    ROOT::VecOps::RVec<float> result;
    result.reserve(btag.size());
    for (auto &p : btag) {
        if (p > m_min_btag) {
            result.emplace_back(p);
        }
    }
    return result;
}


// compute the cone isolation for reco particles
struct coneIsolation {

    coneIsolation(float arg_dr_min, float arg_dr_max);
    double deltaR(double eta1, double phi1, double eta2, double phi2) { return TMath::Sqrt(TMath::Power(eta1-eta2, 2) + (TMath::Power(phi1-phi2, 2))); };

    float dr_min = 0;
    float dr_max = 0.4;
    Vec_f operator() (Vec_rp in, Vec_rp rps) ;
};

coneIsolation::coneIsolation(float arg_dr_min, float arg_dr_max) : dr_min(arg_dr_min), dr_max( arg_dr_max ) { };
Vec_f coneIsolation::coneIsolation::operator() (Vec_rp in, Vec_rp rps) {
  
    Vec_f result;
    result.reserve(in.size());

    std::vector<ROOT::Math::PxPyPzEVector> lv_reco;
    std::vector<ROOT::Math::PxPyPzEVector> lv_charged;
    std::vector<ROOT::Math::PxPyPzEVector> lv_neutral;

    for(size_t i = 0; i < rps.size(); ++i) {

        ROOT::Math::PxPyPzEVector tlv;
        tlv.SetPxPyPzE(rps.at(i).momentum.x, rps.at(i).momentum.y, rps.at(i).momentum.z, rps.at(i).energy);
        
        if(rps.at(i).charge == 0) lv_neutral.push_back(tlv);
        else lv_charged.push_back(tlv);
    }
    
    for(size_t i = 0; i < in.size(); ++i) {

        ROOT::Math::PxPyPzEVector tlv;
        tlv.SetPxPyPzE(in.at(i).momentum.x, in.at(i).momentum.y, in.at(i).momentum.z, in.at(i).energy);
        lv_reco.push_back(tlv);
    }

    
    // compute the isolation (see https://github.com/delphes/delphes/blob/master/modules/Isolation.cc#L154) 
    for (auto & lv_reco_ : lv_reco) {
    
        double sumNeutral = 0.0;
        double sumCharged = 0.0;
    
        // charged
        for (auto & lv_charged_ : lv_charged) {
    
            double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_charged_.Eta(), lv_charged_.Phi());
            if(dr > dr_min && dr < dr_max) sumCharged += lv_charged_.P();
        }
        
        // neutral
        for (auto & lv_neutral_ : lv_neutral) {
    
            double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_neutral_.Eta(), lv_neutral_.Phi());
            if(dr > dr_min && dr < dr_max) sumNeutral += lv_neutral_.P();
        }
        
        double sum = sumCharged + sumNeutral;
        double ratio= sum / lv_reco_.P();
        result.emplace_back(ratio);
    }
    return result;
}
 
 
 
// returns missing energy vector, based on reco particles
Vec_rp missingEnergy(float ecm, Vec_rp in, float p_cutoff = 0.0) {
    float px = 0, py = 0, pz = 0, e = 0;
    for(auto &p : in) {
        if (std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y*p.momentum.y) < p_cutoff) continue;
        px += -p.momentum.x;
        py += -p.momentum.y;
        pz += -p.momentum.z;
        e += p.energy;
    }
    
    Vec_rp ret;
    rp res;
    res.momentum.x = px;
    res.momentum.y = py;
    res.momentum.z = pz;
    res.energy = ecm-e;
    ret.emplace_back(res);
    return ret;
}

// calculate the cosine(theta) of the missing energy vector
float get_cosTheta_miss(Vec_rp met){
    
    float costheta = 0.;
    if(met.size() > 0) {
        
        TLorentzVector lv_met;
        lv_met.SetPxPyPzE(met[0].momentum.x, met[0].momentum.y, met[0].momentum.z, met[0].energy);
        costheta = fabs(std::cos(lv_met.Theta()));
    }
    return costheta;
}

// get objects from indexs
ROOT::VecOps::RVec<edm4hep::MCParticleData>
get_mc(ROOT::VecOps::RVec<int> indexes,
    ROOT::VecOps::RVec<edm4hep::MCParticleData> inParticles) {
  ROOT::VecOps::RVec<edm4hep::MCParticleData> result;

  for (const auto &index : indexes) {
    if (index > -1)
      result.push_back(inParticles.at(index));
  }

  return result;
}


struct TwoParticleGroups {
    int high_mass_idx = -1;
    int low_mass_idx  = -1;
};

// given 2 hard W objects and the full Particle collection, return the indices
// (in Particle) of the on-shell (higher mass) and off-shell (lower mass) W
TwoParticleGroups
get_on_and_off_shell_WW_160ecm(ROOT::VecOps::RVec<edm4hep::MCParticleData> hard_ws,
    ROOT::VecOps::RVec<edm4hep::MCParticleData> all_particles,
    double mass_threshold,
    double mass_width)
{
    const auto& first  = hard_ws[0];
    const auto& second = hard_ws[1];
    const auto& on_shell  = (first.mass >= second.mass) ? first : second;
    const auto& off_shell = (first.mass >= second.mass) ? second : first;

    TwoParticleGroups result;
    for (int i = 0; i < (int)all_particles.size(); i++) {
        const auto& p = all_particles[i];
        bool match = (p.PDG == on_shell.PDG &&
                      p.momentum.x == on_shell.momentum.x &&
                      p.momentum.y == on_shell.momentum.y &&
                      p.momentum.z == on_shell.momentum.z);
        if (match && result.high_mass_idx < 0) { result.high_mass_idx = i; continue; }

        match = (p.PDG == off_shell.PDG &&
                 p.momentum.x == off_shell.momentum.x &&
                 p.momentum.y == off_shell.momentum.y &&
                 p.momentum.z == off_shell.momentum.z);
        if (match && result.low_mass_idx < 0) result.low_mass_idx = i;

        if (result.high_mass_idx >= 0 && result.low_mass_idx >= 0) break;
    }
    return result;
}



ROOT::VecOps::RVec<double> get_pair_masses(
    const ROOT::VecOps::RVec<int>& quark_idxs,
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>>& pairs,
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles)
{
    ROOT::VecOps::RVec<double> masses;
    for (size_t i = 0; i < pairs[0].size(); i++) {
        // compute invariant mass of pair i

        const auto& a_pair = pairs[0][i];
        const auto& b_pair = pairs[1][i];

        int idx_a = quark_idxs[a_pair];
        int idx_b = quark_idxs[b_pair];

        auto pa = particles[idx_a];
        auto pb = particles[idx_b];

        // momentum components
        double px_a = pa.momentum.x;
        double py_a = pa.momentum.y;
        double pz_a = pa.momentum.z;

        double px_b = pb.momentum.x;
        double py_b = pb.momentum.y;
        double pz_b = pb.momentum.z;

        // energy from E^2 = |p|^2 + m^2
        double E_a = std::sqrt(px_a*px_a + py_a*py_a + pz_a*pz_a + pa.mass*pa.mass);
        double E_b = std::sqrt(px_b*px_b + py_b*py_b + pz_b*pz_b + pb.mass*pb.mass);

        double m_pair = std::sqrt(
            std::max(0.0, (E_a + E_b)*(E_a + E_b)
                - (px_a + px_b)*(px_a + px_b)
                - (py_a + py_b)*(py_a + py_b)
                - (pz_a + pz_b)*(pz_a + pz_b))
        );

        masses.push_back(m_pair);

    }
    return masses;
}

// Walk decay chain from w_idx, following W→W copies, until reaching the W
// whose daughters are not another W (i.e. the one that actually decays to quarks/leptons).
int get_decaying_W_idx(int w_idx,
                       const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles,
                       const ROOT::VecOps::RVec<int>& daughter_indices) {
    if (w_idx < 0) return w_idx;
    int current = w_idx;
    for (int depth = 0; depth < 20; depth++) {
        int dbegin = particles[current].daughters_begin;
        int dend   = particles[current].daughters_end;
        int next = -1;
        for (int d = dbegin; d < dend; d++) {
            int didx = daughter_indices[d];
            if (std::abs(particles[didx].PDG) == 24) { next = didx; break; }
        }
        if (next < 0) break;
        current = next;
    }
    return current;
}

// gets the best pair to match to a W and returns their particle indexs
// in ours, we compare to the on shell W

struct OnOffidx {
    ROOT::VecOps::RVec<int> on_shell_idx;
    ROOT::VecOps::RVec<int> off_shell_idx;
    double on_shell_mass = 0.0;
    double off_shell_mass = 0.0;
    ROOT::VecOps::RVec<int> on_shell_flavor;
    ROOT::VecOps::RVec<int> off_shell_flavor;
    int match_truth = 0;
};

// quark_idxs are passed as w_on shell, w_off shell, so the first two indicies should be the on shell, so we do a check
OnOffidx compare_pair_mass_to_w(
    const ROOT::VecOps::RVec<double>& masses,
    const ROOT::VecOps::RVec<int>& quark_idxs,
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>>& pairs,
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles,
    int w_idx)
{
    double diff = 999999;
    int best_idx = -1;
    double w_mass = particles[w_idx].mass;
    OnOffidx result;

    for (size_t i = 0; i < masses.size(); i++) {
        double d = std::abs(masses[i] - w_mass);
        if (d < diff) {
            diff = d;
            best_idx = i;
        }
    }

    size_t pos1 = pairs[0][best_idx];
    size_t pos2 = pairs[1][best_idx];

    result.on_shell_idx = ROOT::VecOps::RVec<int>{quark_idxs[pos1], quark_idxs[pos2]};
    result.on_shell_flavor = ROOT::VecOps::RVec<int>{particles[quark_idxs[pos1]].PDG, particles[quark_idxs[pos2]].PDG};
    result.on_shell_mass = masses[best_idx];

    // positions 0,1 in quark_idxs are the on-shell W's quarks (Concatenate order)
    result.match_truth = ((pos1 < 2) && (pos2 < 2)) ? 1 : 0;



    ROOT::VecOps::RVec<int> off_shell_idx;
    ROOT::VecOps::RVec<int> off_shell_flavor;
    for (size_t j = 0; j < quark_idxs.size(); j++) {
        if (j != pos1 && j != pos2) {
            off_shell_idx.push_back(quark_idxs[j]);
            off_shell_flavor.push_back(particles[quark_idxs[j]].PDG);
        }
    }
    result.off_shell_idx    = off_shell_idx;
    result.off_shell_flavor = off_shell_flavor;

    for (size_t i = 0; i < pairs[0].size(); i++) {
    if (pairs[0][i] != pos1 && pairs[0][i] != pos2 &&
        pairs[1][i] != pos1 && pairs[1][i] != pos2) {
        result.off_shell_mass = masses[i];
        break;
        }
    }

    return result;
}







// modify get_indices to not just take the first decay found since we are interested in ww process

// get_indices::get_indices_all_occurances( int pdg_mother, std::vector<int> pdg_daughters, bool stableDaughters, bool chargeConjugateMother, bool chargeConjugateDaughters, bool inclusiveDecay) {
//   m_pdg_mother = pdg_mother;
//   m_pdg_daughters = pdg_daughters;
//   m_stableDaughters = stableDaughters;
//   m_chargeConjugateMother = chargeConjugateMother;
//   m_chargeConjugateDaughters = chargeConjugateDaughters;
//   m_inclusiveDecay = inclusiveDecay;
// } ;

// ROOT::VecOps::RVec<int>  get_indices::operator() ( ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind) {

//    // Look for a specific decay specified by the mother PDG_id and
//    // the PDG_ids of the daughters
//    // Returns a vector with the indices, in the Particle block, of the mother and of
//    // the daughters - in the order defined by std::vector<int> pdg_daughters.
//    //
//    // In case there are several such decays in the event, keep only the first one.

//    ROOT::VecOps::RVec<int>  result;

//    for ( int imother =0; imother < in.size(); imother ++){
//      int pdg = in[imother].PDG ;
//      bool found_a_mother = false;
//      if ( ! m_chargeConjugateMother ) found_a_mother = ( pdg == m_pdg_mother );
//      if ( m_chargeConjugateMother )   found_a_mother = ( abs(pdg) == abs(m_pdg_mother) ) ;
//      if ( ! found_a_mother ) continue;

//      ROOT::VecOps::RVec<int> a = get_indices_MotherByIndex( imother, m_pdg_daughters, m_stableDaughters, m_chargeConjugateDaughters, m_inclusiveDecay, in, ind );
//      if ( a.size() != 0 ) {
//         result = a;
//         break;    // return the first decay found
//      }

//    }
//    return result;
// }


}}

#endif