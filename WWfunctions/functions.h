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
    double high_mass  = -1;
    double low_mass   = -1;
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

    result.high_mass = on_shell.mass;
    result.low_mass = off_shell.mass;

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
    double best_chi2 = 0.0;
    double second_best_chi2 = 0.0;
    double third_best_chi2 = 0.0;
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

    // for now, consider the distribution of only the mass we are choosing for
    // do a chi squared on only one variable, the mass of the w boson we are interested in
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



// quark_idxs are passed as w_on shell, w_off shell, so the first two indicies should be the on shell, so we do a check
OnOffidx chi2_compare_pair_mass_to_w(
    const ROOT::VecOps::RVec<double>& masses,
    const ROOT::VecOps::RVec<int>& quark_idxs,
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>>& pairs,
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles,
    double w_var,
    int w_idx)
{
    double best_chi2        = 999999;
    double second_best_chi2 = 999999;
    double third_best_chi2  = 999999;
    int best_idx = -1;
    //double w_mass = particles[w_idx].mass;
    double w_mass = 80.419;
    OnOffidx result;

    for (size_t i = 0; i < masses.size(); i++) {
        double chi2 = ((masses[i] - w_mass) * (masses[i] - w_mass)) / (w_var*w_var);

        if (chi2 < best_chi2) {
            third_best_chi2  = second_best_chi2;
            second_best_chi2 = best_chi2;
            best_chi2 = chi2;
            best_idx = i;
        } else if (chi2 < second_best_chi2) {
            third_best_chi2  = second_best_chi2;
            second_best_chi2 = chi2;
        } else if (chi2 < third_best_chi2) {
            third_best_chi2 = chi2;
        }
    }

    size_t pos1 = pairs[0][best_idx];
    size_t pos2 = pairs[1][best_idx];

    result.on_shell_idx = ROOT::VecOps::RVec<int>{quark_idxs[pos1], quark_idxs[pos2]};
    result.on_shell_flavor = ROOT::VecOps::RVec<int>{particles[quark_idxs[pos1]].PDG, particles[quark_idxs[pos2]].PDG};
    result.on_shell_mass = masses[best_idx];
    result.best_chi2        = best_chi2;
    result.second_best_chi2 = second_best_chi2;
    result.third_best_chi2  = third_best_chi2;

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

// create function that takes two collections of 4 vectors of the same size, 
// and then matches them based on the smallest deltaR value first, then second smallest of the remaining and so on,
// optionally, can input an overall minimum deltaR that is imposed on the last match, and therefor all of the
// matches before too
// in this case we match reco jets to mc quarks


//ROOT::VecOps::RVec<TLorentzVector>

struct JetToQuarkInfo {
    ROOT::VecOps::RVec<int> idx;
    ROOT::VecOps::RVec<double> delta_Rs;
    ROOT::VecOps::RVec<double> delta_etas;
    ROOT::VecOps::RVec<double> delta_phis;
    int under_min_delR = 0;
    
};



// delta_R vector is not ordered the same as idx vector, it returns in order of how we chose jets, so we always get the minimum quark-jet combo first 
// want to get delta_phi and delta_eta info too in order to understand their distributions
JetToQuarkInfo MatchJetsToQuarks(
    const ROOT::VecOps::RVec<TLorentzVector>& jets,   // changed type
    const ROOT::VecOps::RVec<TLorentzVector>& quarks,
    double delR_constraint)
{
    JetToQuarkInfo result;

    int n = jets.size();
    result.idx = ROOT::VecOps::RVec<int>(n, -1);
    
    std::vector<bool> quark_used(quarks.size(), false);

    for (int iter = 0; iter < n; iter++) {
        double min_dr = 999;
        double chosen_deta = 999;
        double chosen_dphi = 999;
        int best_jet = -1, best_quark = -1;
        for (int j = 0; j < n; j++) {
            if (result.idx[j] >= 0) continue;
            for (int q = 0; q < (int)quarks.size(); q++) {

                if (quark_used[q]) continue;
                double deta = quarks[q].Eta() - jets[j].Eta();
                double dphi = TVector2::Phi_mpi_pi(quarks[q].Phi() - jets[j].Phi());
                double dr   = std::sqrt(deta*deta + dphi*dphi);
                if (dr < min_dr) { 
                    min_dr = dr; 
                    chosen_deta = deta;
                    chosen_dphi = dphi;
                    best_jet = j; 
                    best_quark = q; }

            }
        }
        result.idx[best_jet] = best_quark;
        result.delta_Rs.push_back(min_dr);
        result.delta_etas.push_back(chosen_deta);
        result.delta_phis.push_back(chosen_dphi);

        quark_used[best_quark] = true;

        if (iter == n - 1 && min_dr < delR_constraint) {
            result.under_min_delR = 1;
        }
        
    }
    return result;
}


// modified dijet mass calculator given 4vecs of jets
// want to modify to give pair indexs as well



struct DijetInfo {
    ROOT::VecOps::RVec<double> masses;
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>> pair_idxs;
};

DijetInfo all_invariant_masses_and_pair_idxs(
    ROOT::VecOps::RVec<TLorentzVector> AllJets)
{
    TLorentzVector tlv1;
    TLorentzVector tlv2;
    double E, px, py, pz; 
    double invmass; 

    DijetInfo result;

    ROOT::VecOps::RVec<size_t> first_idxs, second_idxs;

    if(AllJets.size() < 2) {
        result.masses = ROOT::VecOps::RVec<double>{};
        result.pair_idxs = ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>>{};
        return result;
    }

    // For each jet, take its invariant mass with the remaining jets. Stop at last jet.
    for(int i = 0; i < AllJets.size()-1; ++i) {

    tlv1 = AllJets.at(i); 

        for(int j=i+1; j < AllJets.size(); ++j){ // go until end
            tlv2 = AllJets.at(j);
            E = tlv1.E() + tlv2.E();
            px = tlv1.Px() + tlv2.Px();
            py = tlv1.Py() + tlv2.Py();
            pz = tlv1.Pz() + tlv2.Pz();

            invmass = std::sqrt(E*E - px*px - py*py - pz*pz);
            result.masses.push_back(invmass);
            first_idxs.push_back(i);
            second_idxs.push_back(j);
            
        }
    }

    result.pair_idxs.push_back(first_idxs);
    result.pair_idxs.push_back(second_idxs);

    return result;
}    

// function used to map one type of index to another for a pair

ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>>
transform_pair_idxs(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>>& pair_idxs,
    const ROOT::VecOps::RVec<int>& jet_to_quark)
{
    ROOT::VecOps::RVec<size_t> a, b;
    for (size_t k = 0; k < pair_idxs[0].size(); k++) {
        a.push_back(jet_to_quark[pair_idxs[0][k]]);
        b.push_back(jet_to_quark[pair_idxs[1][k]]);
    }
    return {a, b};
}

// given an array of size n, put out all combinations of assingments, ie 
// config[0] = [0,1,2,3] 
// config[1] = [0,1,3,2] 
// config[2] = [0,2,1,3] 
// ...                       (24 total)
std::vector<std::vector<int>> GetAllPermutations(int n) {
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);  // fill [0, 1, ..., n-1]

    std::vector<std::vector<int>> all;
    do {
        all.push_back(perm);
    } while (std::next_permutation(perm.begin(), perm.end()));

    return all;  // n! entries, each of length n
}


// function that takes permutations and creates the quark jet pairs we need using their tlvs

// takes a list of lists of observed values grouped by event ({deta1, dphi1}, {deta2, dphi2})
// and a list for sigmas with each sigma corresponding to a type of observable
double ChiSquared(
    const ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>& observed_values,
    const ROOT::VecOps::RVec<double>& sigmas)
{
    double result = 0;
    for (int i = 0; i < (int)observed_values.size(); i++) {
        for (int k = 0; k < (int)observed_values[i].size(); k++) {
            result += (observed_values[i][k] / sigmas[k]) * (observed_values[i][k] / sigmas[k]);
        }
    }
    return result;
}





struct JtoQ_dR_Info {
    ROOT::VecOps::RVec<int> idx;
    ROOT::VecOps::RVec<double> delta_Rs;
    int under_min_delR = 0;
    double best_chi2 = 0.0;
    double second_best_chi2 = 0.0;
    double third_best_chi2 = 0.0;
};

struct JtoQ_etaphi_Info {
    ROOT::VecOps::RVec<int> idx;
    ROOT::VecOps::RVec<double> delta_Rs;
    ROOT::VecOps::RVec<double> delta_etas;
    ROOT::VecOps::RVec<double> delta_phis;
    int under_min_delR = 0;
    double best_chi2        = 0.0;
    double second_best_chi2 = 0.0;
    double third_best_chi2  = 0.0;
};

// only uses deltaR for chi squared
JtoQ_dR_Info JtoQ_ChiSquared_deltaR(
    const ROOT::VecOps::RVec<TLorentzVector>& jets,
    const ROOT::VecOps::RVec<TLorentzVector>& quarks,
    const ROOT::VecOps::RVec<double>& sigmas,
    double delR_constraint)
{
    JtoQ_dR_Info result;

    std::vector<std::vector<int>> all_permutation_idxs = GetAllPermutations(jets.size());

    double best_chi2        = std::numeric_limits<double>::max();
    double second_best_chi2 = std::numeric_limits<double>::max();
    double third_best_chi2 = std::numeric_limits<double>::max();


    for (const auto& perm : all_permutation_idxs) {

        ROOT::VecOps::RVec<double> delta_Rs;
        bool all_within = true;

        for (int j = 0; j < (int)jets.size(); j++) {
            int q = perm[j];
            double dR = jets[j].DeltaR(quarks[q]);
            delta_Rs.push_back(dR);
            if (dR > delR_constraint) all_within = false;
        }

        ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> obs;
        for (auto& v : delta_Rs) obs.push_back({v});
        double chi2 = ChiSquared(obs, sigmas);

        if (chi2 < best_chi2) {
            second_best_chi2        = best_chi2;
            best_chi2               = chi2;
            result.second_best_chi2 = second_best_chi2;
            result.best_chi2        = chi2;
            result.delta_Rs         = delta_Rs;
            result.under_min_delR   = all_within ? 1 : 0;
            result.idx.clear();
            for (int j = 0; j < (int)perm.size(); j++)
                result.idx.push_back(perm[j]);
        } else if (chi2 < second_best_chi2) {
            second_best_chi2        = chi2;
            result.second_best_chi2 = chi2;
        } else if (chi2 < third_best_chi2) {
            third_best_chi2         = chi2;
            result.third_best_chi2  = chi2;
        }
    }

    return result;
}


// uses delta_eta and delta_phi for chi squared
JtoQ_etaphi_Info JtoQ_ChiSquared_eta_phi(
    const ROOT::VecOps::RVec<TLorentzVector>& jets,
    const ROOT::VecOps::RVec<TLorentzVector>& quarks,
    const ROOT::VecOps::RVec<double>& sigmas,
    double delR_constraint)
{
    JtoQ_etaphi_Info result;

    std::vector<std::vector<int>> all_permutation_idxs = GetAllPermutations(jets.size());

    double best_chi2        = std::numeric_limits<double>::max();
    double second_best_chi2 = std::numeric_limits<double>::max();
    double third_best_chi2  = std::numeric_limits<double>::max();

    for (const auto& perm : all_permutation_idxs) {

        ROOT::VecOps::RVec<double> delta_Rs, delta_etas, delta_phis;
        bool all_within = true;

        for (int j = 0; j < (int)jets.size(); j++) {
            int q       = perm[j];
            double deta = quarks[q].Eta() - jets[j].Eta();
            double dphi = TVector2::Phi_mpi_pi(quarks[q].Phi() - jets[j].Phi());
            double dR   = std::sqrt(deta*deta + dphi*dphi);
            delta_etas.push_back(deta);
            delta_phis.push_back(dphi);
            delta_Rs.push_back(dR);
            if (dR > delR_constraint) all_within = false;
        }

        ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> obs;
        for (int j = 0; j < (int)jets.size(); j++)
            obs.push_back({delta_etas[j], delta_phis[j]});
        double chi2 = ChiSquared(obs, sigmas);

        if (chi2 < best_chi2) {
            third_best_chi2         = second_best_chi2;
            second_best_chi2        = best_chi2;
            best_chi2               = chi2;
            result.third_best_chi2  = third_best_chi2;
            result.second_best_chi2 = second_best_chi2;
            result.best_chi2        = chi2;
            result.delta_Rs         = delta_Rs;
            result.delta_etas       = delta_etas;
            result.delta_phis       = delta_phis;
            result.under_min_delR   = all_within ? 1 : 0;
            result.idx.clear();
            for (int j = 0; j < (int)perm.size(); j++)
                result.idx.push_back(perm[j]);
        } else if (chi2 < second_best_chi2) {
            third_best_chi2         = second_best_chi2;
            second_best_chi2        = chi2;
            result.third_best_chi2  = third_best_chi2;
            result.second_best_chi2 = chi2;
        } else if (chi2 < third_best_chi2) {
            third_best_chi2         = chi2;
            result.third_best_chi2  = chi2;
        }
    }

    return result;
}


// want to look at ""incorrect" pairs of Ws, so we want to mix the quarks that are clustered from W1 and W2 and get their permutations, and ultimately their masses

// ROOT::VecOps::RVec<double> get_pair_masses(
//     const ROOT::VecOps::RVec<int>& quark_idxs,
//     const ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>>& pairs,
//     const ROOT::VecOps::RVec<edm4hep::MCParticleData>& particles)


struct MixQuarkPairsInfo {
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::MCParticleData>> perm1;
    ROOT::VecOps::RVec<ROOT::VecOps::RVec<edm4hep::MCParticleData>> perm2;
    ROOT::VecOps::RVec<double> perm1_mass;
    ROOT::VecOps::RVec<double> perm2_mass;
};


//inputs are randomly mixed because the quarks from each W are sorted by how they are identified by their PDG,
//  so we are biased more towards certain quarks than others without random mixing
MixQuarkPairsInfo MixQuarkPairsAndGetMass(
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& quarks_from_W1,
    const ROOT::VecOps::RVec<edm4hep::MCParticleData>& quarks_from_W2)
{
    MixQuarkPairsInfo result;
    if (quarks_from_W1.size() < 2 || quarks_from_W2.size() < 2) return result;

    // randomly flip ordering of each W's quarks independently to remove flavor bias
    thread_local std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> coin(0, 1);
    bool flip_W1 = coin(rng);
    bool flip_W2 = coin(rng);

    const auto& w1_q0 = flip_W1 ? quarks_from_W1[1] : quarks_from_W1[0];
    const auto& w1_q1 = flip_W1 ? quarks_from_W1[0] : quarks_from_W1[1];
    const auto& w2_q0 = flip_W2 ? quarks_from_W2[1] : quarks_from_W2[0];
    const auto& w2_q1 = flip_W2 ? quarks_from_W2[0] : quarks_from_W2[1];

    auto pair_inv_mass = [](const edm4hep::MCParticleData& a, const edm4hep::MCParticleData& b) -> double {
        double E_a = std::sqrt(a.mass*a.mass + a.momentum.x*a.momentum.x + a.momentum.y*a.momentum.y + a.momentum.z*a.momentum.z);
        double E_b = std::sqrt(b.mass*b.mass + b.momentum.x*b.momentum.x + b.momentum.y*b.momentum.y + b.momentum.z*b.momentum.z);
        double E  = E_a + E_b;
        double px = a.momentum.x + b.momentum.x;
        double py = a.momentum.y + b.momentum.y;
        double pz = a.momentum.z + b.momentum.z;
        return std::sqrt(std::max(0.0, E*E - px*px - py*py - pz*pz));
    };

    // mixed pairing 1: w1_q0+w2_q0, w1_q1+w2_q1
    result.perm1.push_back({w1_q0, w2_q0});
    result.perm2.push_back({w1_q1, w2_q1});
    result.perm1_mass.push_back(pair_inv_mass(w1_q0, w2_q0));
    result.perm2_mass.push_back(pair_inv_mass(w1_q1, w2_q1));

    // mixed pairing 2: w1_q0+w2_q1, w1_q1+w2_q0
    result.perm1.push_back({w1_q0, w2_q1});
    result.perm2.push_back({w1_q1, w2_q0});
    result.perm1_mass.push_back(pair_inv_mass(w1_q0, w2_q1));
    result.perm2_mass.push_back(pair_inv_mass(w1_q1, w2_q0));

    return result;
}



struct DRandDTheta {
    ROOT::VecOps::RVec<double> DR;
    ROOT::VecOps::RVec<double> DTheta;
};

DRandDTheta GetAngularInfoFromDijet(
    const TLorentzVector& lv1,
    const TLorentzVector& lv2)
{
    DRandDTheta result;

    result.DR.push_back(lv1.DeltaR(lv2));
    result.DTheta.push_back(lv1.Angle(lv2.Vect()));

    return result;
}








struct WJetKinematics {
    double Wa_dR, Wa_dTheta, Wa_P, Wa_Pt, Wa_E;
    double Wb_dR, Wb_dTheta, Wb_P, Wb_Pt, Wb_E;
};

inline WJetKinematics compute_W_jet_kinematics(
    const TLorentzVector& j1, const TLorentzVector& j2,
    const TLorentzVector& j3, const TLorentzVector& j4,
    int wlab1, int wlab2, int wlab3, int wlab4)
{
    TLorentzVector ja1, ja2, jb1, jb2;
    int na = 0, nb = 0;
    auto fill = [](TLorentzVector& out1, TLorentzVector& out2, int& n, const TLorentzVector& j) {
        if (n++ == 0) out1 = j; else out2 = j;
    };
    if (wlab1 == 0) fill(ja1, ja2, na, j1); else fill(jb1, jb2, nb, j1);
    if (wlab2 == 0) fill(ja1, ja2, na, j2); else fill(jb1, jb2, nb, j2);
    if (wlab3 == 0) fill(ja1, ja2, na, j3); else fill(jb1, jb2, nb, j3);
    if (wlab4 == 0) fill(ja1, ja2, na, j4); else fill(jb1, jb2, nb, j4);

    WJetKinematics r;
    r.Wa_dR     = ja1.DeltaR(ja2);
    r.Wa_dTheta = ja1.Angle(ja2.Vect());
    r.Wa_P      = (ja1 + ja2).P();
    r.Wa_Pt     = (ja1 + ja2).Pt();
    r.Wa_E      = (ja1 + ja2).E();
    r.Wb_dR     = jb1.DeltaR(jb2);
    r.Wb_dTheta = jb1.Angle(jb2.Vect());
    r.Wb_P      = (jb1 + jb2).P();
    r.Wb_Pt     = (jb1 + jb2).Pt();
    r.Wb_E      = (jb1 + jb2).E();
    return r;
}

struct HJetKinematics {
    double H_dR, H_dTheta, H_P, H_Pt, H_E;
};

inline HJetKinematics compute_H_jet_kinematics(
    const TLorentzVector& j1, const TLorentzVector& j2)
{
    HJetKinematics r;
    r.H_dR     = j1.DeltaR(j2);
    r.H_dTheta = j1.Angle(j2.Vect());
    r.H_P      = (j1 + j2).P();
    r.H_Pt     = (j1 + j2).Pt();
    r.H_E      = (j1 + j2).E();
    return r;
}


struct sel_quarks_fromHiggs {
    int boson_pdg;
    explicit sel_quarks_fromHiggs(int pdg = 25) : boson_pdg(pdg) {}
    ROOT::VecOps::RVec<edm4hep::MCParticleData> operator()(
        ROOT::VecOps::RVec<edm4hep::MCParticleData> in,
        const ROOT::VecOps::RVec<int>& parents_relation) const {
        // Pythia8 can produce multiple copies of the boson with different indices.
        // Grouping by parent index breaks in that case, so just collect all quarks
        // that have any direct parent with the target PDG, then require exactly 2.
        ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
        for (size_t i = 0; i < in.size(); ++i) {
            const auto& p = in[i];
            if (std::abs(p.PDG) > 6 || p.PDG == 0) continue;
            for (unsigned j = p.parents_begin; j < p.parents_end; ++j) {
                if (j >= parents_relation.size()) break;
                int idx = parents_relation[j];
                if (idx < 0 || idx >= (int)in.size()) continue;
                if (std::abs(in[idx].PDG) == boson_pdg) { result.emplace_back(p); break; }
            }
        }
        if (result.size() != 2) result.clear();
        return result;                                  // size 2 or empty
    }
};


struct sel_nunu_fromHiggs {
    int boson_pdg;
    explicit sel_nunu_fromHiggs(int pdg = 25) : boson_pdg(pdg) {}
    ROOT::VecOps::RVec<edm4hep::MCParticleData> operator()(
        ROOT::VecOps::RVec<edm4hep::MCParticleData> in,
        const ROOT::VecOps::RVec<int>& parents_relation) const {
        // Pythia8 can produce multiple copies of the boson with different indices.
        // Grouping by parent index breaks in that case, so just collect all quarks
        // that have any direct parent with the target PDG, then require exactly 2.
        ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
        for (size_t i = 0; i < in.size(); ++i) {
            const auto& p = in[i];
            if ((std::abs(p.PDG) != 12) & (std::abs(p.PDG) != 14) & (std::abs(p.PDG) != 16)) continue;
            for (unsigned j = p.parents_begin; j < p.parents_end; ++j) {
                if (j >= parents_relation.size()) break;
                int idx = parents_relation[j];
                if (idx < 0 || idx >= (int)in.size()) continue;
                if (std::abs(in[idx].PDG) == boson_pdg) { result.emplace_back(p); break; }
            }
        }
        if (result.size() != 2) result.clear();
        return result;                                  // size 2 or empty
    }
};

struct sel_daughter_fromParent {
    int parent_pdg;
    int num_daughters;
    ROOT::VecOps::RVec<int> daughter_pdgs;

    explicit sel_daughter_fromParent(
        int parent_pdg,
        int num_daughters,
        ROOT::VecOps::RVec<int> daughter_pdgs
    ) : parent_pdg(parent_pdg), num_daughters(num_daughters), daughter_pdgs(daughter_pdgs) {}
    ROOT::VecOps::RVec<edm4hep::MCParticleData> operator()(
        ROOT::VecOps::RVec<edm4hep::MCParticleData> in,
        ROOT::VecOps::RVec<int>& parents_relation) const {
        // Pythia8 can produce multiple copies of the parent boson with different indices.
        // Grouping by parent index breaks in that case, so collect all particles whose
        // |PDG| is in daughter_pdgs and have any direct parent with |PDG| == parent_pdg,
        // then require exactly num_daughters matches.
        ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
        for (size_t i = 0; i < in.size(); ++i) {
            const auto& p = in[i];
            bool is_daughter = false;
            for (size_t k = 0; k < daughter_pdgs.size(); ++k) {
                if (std::abs(p.PDG) == daughter_pdgs[k]) { is_daughter = true; break; }
            }
            if (!is_daughter) continue;
            for (unsigned j = p.parents_begin; j < p.parents_end; ++j) {
                if (j >= parents_relation.size()) break;
                int idx = parents_relation[j];
                if (idx < 0 || idx >= (int)in.size()) continue;
                if (std::abs(in[idx].PDG) == parent_pdg) { result.emplace_back(p); break; }
            }
        }
        if (result.size() != num_daughters) result.clear();
        return result;  // size num_daughters or empty
    }
};





}
}

#endif
