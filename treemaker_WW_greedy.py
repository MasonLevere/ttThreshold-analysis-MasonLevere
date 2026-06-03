import os, copy
import urllib
from itertools import combinations

# list of processes
all_processes = {

#################################
#  PROCESS TESTED ON PREVIOUSLY
#
#    "p8_ee_WW_ecm345": {
#        "fraction": 1,
#    },
#
#################################
################################
# EXCLUSIVE SEMIHADRONIC SAMPLE
    # "wzp6_ee_munumuqq_noCut_ecm163": {
    #     "fraction": 1,
    # },
################################

# One on shell, one off shell WW events
    "p8_ee_WW_ecm160": {
        "fraction": 1,
    },

###FOR foll WW background samples with PS variations xsec values are wrong in the database!!
#   "p8_ee_WW_PSdown_ecm340":{ "fraction": 1,},
#   "p8_ee_WW_PSup_ecm340":{ "fraction": 1,},
#   "p8_ee_WW_PSdown_ecm345":{ "fraction": 1,},
#   "p8_ee_WW_PSup_ecm345":{ "fraction": 1,},
#   "p8_ee_WW_PSdown_ecm365":{ "fraction": 1,},
#   "p8_ee_WW_PSup_ecm365":{ "fraction": 1,},

###### following samples we don't use
##am     "wzp6_ee_WbWb_semihad_ecm345": {
##am         "fraction": 1,
##am     },
##am     "wzp6_ee_WbWb_had_ecm345": {
##am         "fraction": 1,
##am     },
##am    "wzp6_ee_WbWb_semihad_ecm350": {
##am       "fraction": 1,
##am    },
##am    "wzp6_ee_WbWb_had_ecm350": {
##am       "fraction": 1,
##am    },
##am    "wzp6_ee_WbWb_semihad_ecm355": {
##am        "fraction": 1,
##am    },
##am    "wzp6_ee_WbWb_had_ecm355": {
##am        "fraction": 1,
##am    },
##am    "wzp6_ee_WbWb_semihad_ecm340": {
##am        "fraction": 1,
##am    },
##am    "wzp6_ee_WbWb_had_ecm340": {
##am        "fraction": 1,
##am    },
##am    "wzp6_ee_WbWb_lep_ecm340": {
##am     "fraction": 1,
##am     },
##am    "wzp6_ee_WbWb_lep_ecm345": {
##am             "fraction": 1,
##am     },
##am    "wzp6_ee_WbWb_lep_ecm350": {
##am     "fraction": 1,
##am     },
##am    "wzp6_ee_WbWb_lep_ecm355": {
##am             "fraction": 1,
##am     },
##am    "wzp6_ee_WbWb_lep_ecm365": {
##am         "fraction": 1,
##am     },
##am    "wzp6_ee_WbWb_semihad_ecm365": {
##am        "fraction": 1,
##am     },
##am    "wzp6_ee_WbWb_had_ecm365": {
##am         "fraction": 1,
##am    },

}

# ### CONDOR ###
# runBatch = True
# batchQueue = "espresso"
# compGroup = "group_u_FCC.local_gen"




available_ecm = ['125', '160', '163', '340','345', '350', '355','365']

hadronic  = True
#semihad  = False
#lep      = False
ecm       = 160
print(ecm)

saveExclJets = True


### commented out since using new sample
if not str(ecm) in available_ecm:
    raise ValueError("ecm value not in available_ecm")

channel = "hadronic"

if  channel not in ["lep","semihad","had"]:
    print("using defa channel settings")
    #channel="semihad"
print(channel)    

processList={key: value for key, value in all_processes.items() if str(ecm) in available_ecm and str(ecm) in key } # (True if str('p8_ee_WW_ecm'+ecm) in key else str('wzp6_ee_WbWb_ecm'+ecm) in key)}  

print(processList)
# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
#prodTag     = "FCCee/winter2023/IDEA/"

inputDir    = "/eos/user/m/mlevere/ttThreshold-analysis/localSamples/"

#Optional: output directory, default is local running directoryp
outputDir   = "/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/{}_WW_greedy".format(channel)


# additional/costom C++ functions, defined in header files (optional)
includePaths = ["examples/functions.h"]

## latest particle transformer model, trained on 9M jets in winter2023 samples
model_name = "fccee_flavtagging_edm4hep_wc" #"fccee_flavtagging_edm4hep_wc_v1"

## model files locally stored on /eos
eos_dir ="/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
# new model file directory


### CHANGE DIRECTORY

model_dir = (
    "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_7classes_12_04_2023/"    # "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
)
local_preproc = "{}/{}.json".format(model_dir, model_name)
local_model = "{}/{}.onnx".format(model_dir, model_name)

url_model_dir = "https://fccsw.web.cern.ch/fccsw/testsamples/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
url_preproc = "{}/{}.json".format(url_model_dir, model_name)
url_model = "{}/{}.onnx".format(url_model_dir, model_name)

## get local file, else download from url
def get_file_path(url, filename):
    if os.path.exists(filename):
        return os.path.abspath(filename)
    else:
        urllib.request.urlretrieve(url, os.path.basename(url))
        return os.path.basename(url)

def get_files(eos_dir, proc):
    files=[]
    basepath=os.path.join(eos_dir,proc)
    if os.path.exists(basepath):
        files =  [os.path.join(basepath,x) for x in os.listdir(basepath) if os.path.isfile(os.path.join(basepath, x)) ]
    return files

weaver_preproc = get_file_path(url_preproc, local_preproc)
weaver_model = get_file_path(url_model, local_model)



from addons.ONNXRuntime.jetFlavourHelper import JetFlavourHelper
from addons.FastJet.jetClusteringHelper import (
    ExclusiveJetClusteringHelper,
    InclusiveJetClusteringHelper,
)

jetFlavourHelper = None
jetFlavourHelper_R5 = None
jetClusteringHelper = None
jetClusteringHelper_R5 = None



# create all combinations of quark pair names
quark_pdg = {
            'd': 1,
            'u': 2,
            's': 3,
            'c': 4,
            'b': 5,
            't': 6,
        }

lepton_pdg = {
    'e':  (11, 12),
    'mu': (13, 14),
    'tau': (15, 16),
}


w_hadron_decay_names = []
w_lepton_decay_names = []

# greedy matcher — no sigma parameters needed



all_branches = [

    # --- inspect_jets.py ---
    "jets_p4",
    "All_W_quarks_idx", "all_W_quarks_obj",
    "HardWs_all_mass",
    "matched_jets_to_q_idx", "matched_jets_to_q_under_min_delR",
    "dijet_masses", "dijet_pair_idx_a", "dijet_pair_idx_b",
    "Candidate_on_shell_W_qq_mass", "Candidate_off_shell_W_qq_mass",
    "Candidate_reco_on_shell_W_jj_mass", "Candidate_reco_off_shell_W_jj_mass",
    "reco_W_jj_match_truth",

    # --- delR_study.py ---
    "simple_jet_1_deltaR", "simple_jet_2_deltaR", "simple_jet_3_deltaR", "simple_jet_4_deltaR",
    "simple_jet_1_delta_eta", "simple_jet_2_delta_eta", "simple_jet_3_delta_eta", "simple_jet_4_delta_eta",
    "simple_jet_1_delta_phi", "simple_jet_2_delta_phi", "simple_jet_3_delta_phi", "simple_jet_4_delta_phi",

    # --- greedy matching ---
    "W_qq_match_truth",

    # TEMPORARY — not read by any current plotting or inspection script
    # "nlep", "lep_p", 'lep_theta', 'lep_phi',
    # "missing_p", "missing_p_theta", "missing_p_phi", "missing_pt",
    # "jet_p", "jet_e", "jet_mass", "jet_phi", "jet_theta", "jet_nconst", "event_njet",
    # "njets_R5",
    # "nbjets_R5_true", "ncjets_R5_true", "nljets_R5_true", "ngjets_R5_true",
    # "bjet1_R5_true_p", "ljet1_R5_true_p",
    # "jets_R5_p", "jets_R5_theta", "jets_R5_pflavor",
    # "bjets_R5_WPp5", "bjets_R5_WPp8", "bjets_R5_WPp85", "bjets_R5_WPp9",
    # "nbjets_R5_WPp5", "nbjets_R5_WPp8", "nbjets_R5_WPp85", "nbjets_R5_WPp9",
    # "recojet_isG_R5", "recojet_isU_R5", "recojet_isB_R5", "recojet_isS_R5",
    # "recojet_isC_R5", "recojet_isD_R5", "recojet_isTAU_R5",
    # "Ws_all", "HardWs_all", "HardWs_all_energy", "HardWs_all_p",
    # "W_on_shell_idx", "W_off_shell_idx",
    # "on_shell_quark_idxs", "off_shell_quark_idxs",
    # "all_W_quarks_tlv",
    # "matched_jets_to_q_delta_Rs", "matched_jets_to_q_delta_etas", "matched_jets_to_q_delta_phis",
    # "chi2_matched_jets_to_q",
    # "Candidate_reco_on_shell_W_jj_p_idxs", "Candidate_reco_off_shell_W_jj_p_idxs",
    # "Candidate_reco_on_shell_W_jj_flavor", "Candidate_reco_off_shell_W_jj_flavor",
    # "on_shell_quark_objs", "off_shell_quark_objs",
    # "Mass_qq_pairs",
    # "Candidate_on_shell_W_qq_p_idxs", "Candidate_off_shell_W_qq_p_idxs",
    # "Candidate_on_shell_W_qq_p_flavor", "Candidate_off_shell_W_qq_p_flavor",
    # 'chi2_matched_jets_to_q_delta_Rs', 'chi2_matched_jets_to_q_delta_etas', 'chi2_matched_jets_to_q_delta_phis',

] + w_hadron_decay_names + w_lepton_decay_names





#print('saving these branches',all_branches)
# Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis:

    # __________________________________________________________
    # Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):


    

        # __________________________________________________________
        # Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2

        # define some aliases to be used later on
        df = df.Alias("Muon0", "Muon#0.index")
        df = df.Alias("Electron0","Electron#0.index")

        # get all the leptons from the collection
        df = df.Define(
            "muons_all",
            "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)",
        )
        df = df.Define(
            "electrons_all",
            "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)",
        )

        # select leptons with momentum > 12 GeV
        df = df.Define(
            "muons_sel",
            "FCCAnalyses::ReconstructedParticle::sel_p(12)(muons_all)",
        )

        df = df.Define(
            "electrons_sel",
            "FCCAnalyses::ReconstructedParticle::sel_p(12)(electrons_all)",
        )

        # compute the muon isolation and store muons with an isolation cut of 0df = df.25 in a separate column muons_sel_iso
        df = df.Define(
            "muons_iso",
            "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.5)(muons_sel, ReconstructedParticles)",
        )
        # compute the electron isolation and store electrons with an isolation cut of 0df = df.25 in a separate column electrons_sel_iso
        df = df.Define(
            "electrons_iso",
            "FCCAnalyses::ZHfunctions::coneIsolation(0.01, 0.5)(electrons_sel, ReconstructedParticles)",
        )

        df = df.Define(
            "muons_sel_iso",
            "FCCAnalyses::ZHfunctions::sel_iso(0.25)(muons_sel, muons_iso)",
        )

        df = df.Define(
            "electrons_sel_iso",
            "FCCAnalyses::ZHfunctions::sel_iso(0.25)(electrons_sel, electrons_iso)",
        )

        if channel == "hadronic":
            #hadronic=True
            df = df.Filter("muons_sel_iso.size() + electrons_sel_iso.size() == 0")
        # elif  channel == "semihad":
        #     #semihad=True
        #     df = df.Filter("muons_sel_iso.size() + electrons_sel_iso.size() == 1")
        # elif channel == "leptonic":
        #     #lep=True
        #     df = df.Filter("muons_sel_iso.size() + electrons_sel_iso.size() == 2")

        # uncomment to restore channel seperating
        # if not (channel == "had"):
        #     df = df.Define(
        #         "muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_sel_iso)"
        #     )


        #     df = df.Define(
        #         "electrons_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_sel_iso)"
        #     )

        #     df = df.Define(
        #         "muons_theta",
        #         "FCCAnalyses::ReconstructedParticle::get_theta(muons_sel_iso)",
        #     )
        #     df = df.Define(
        #         "muons_phi",
        #         "FCCAnalyses::ReconstructedParticle::get_phi(muons_sel_iso)",
        #     )
        #     df = df.Define(
        #         "muons_q",
        #         "FCCAnalyses::ReconstructedParticle::get_charge(muons_sel_iso)",
        #     )
        #     df = df.Define(
        #         "muons_n", "FCCAnalyses::ReconstructedParticle::get_n(muons_sel_iso)",
        #     )

        #     df = df.Define(
        #         "electrons_theta",
        #         "FCCAnalyses::ReconstructedParticle::get_theta(electrons_sel_iso)",
        #     )
        #     df = df.Define(
        #         "electrons_phi",
        #         "FCCAnalyses::ReconstructedParticle::get_phi(electrons_sel_iso)",
        #         )
        #     df = df.Define(
        #         "electrons_q",
        #         "FCCAnalyses::ReconstructedParticle::get_charge(electrons_sel_iso)",
        #         )
        #     df = df.Define(
        #         "electrons_n", "FCCAnalyses::ReconstructedParticle::get_n(electrons_sel_iso)",
        #     )

        
        ## here cluster jets in the events but first remove muons from the list of
        ## reconstructed particles
        
        ## create a new collection of reconstructed particles removing muons with p>12
        df = df.Define(
            "ReconstructedParticlesNoMuons",
            "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles,muons_sel_iso)",
        )
        df = df.Define(
            "ReconstructedParticlesNoMuNoEl",
            "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticlesNoMuons,electrons_sel_iso)",
        )


        ## perform exclusive jet clustering
        global jetClusteringHelper
        global jetFlavourHelper
        global jetFlavourHelper_R5
        global jetClusteringHelper_R5
        
        ## define jet and run clustering parameters
        ## name of collections in EDM root files
        collections = {
            "GenParticles": "Particle",
            "PFParticles": "ReconstructedParticles",
            "PFTracks": "EFlowTrack",
            "PFPhotons": "EFlowPhoton",
            "PFNeutralHadrons": "EFlowNeutralHadron",
            "TrackState": "EFlowTrack_1",
            "TrackerHits": "TrackerHits",
            "CalorimeterHits": "CalorimeterHits",
            "dNdx": "EFlowTrack_2",
            "PathLength": "EFlowTrack_L",
            "Bz": "magFieldBz",
        }

        # want to get some gen level information too
            # use this function to get quarks of Ws get_indices_MotherByIndex

            # first, select Ws
        df = df.Define(
            "Ws_all",
            "FCCAnalyses::MCParticle::sel_pdgID(24, true)(Particle)"
        )

        df = df.Define(
            "HardWs_all",
            "FCCAnalyses::MCParticle::sel_genStatus(22)(Ws_all)"

        )

        df = df.Filter("HardWs_all.size() == 2", "exactly 2 hard Ws")

        df = df.Define(
            "HardWs_all_mass",
            "FCCAnalyses::MCParticle::get_mass(HardWs_all)"
        )


        df = df.Define(
            "HardWs_all_energy",
            "FCCAnalyses::MCParticle::get_e(HardWs_all)"
        )

        df = df.Define(
            "HardWs_all_p",
            "FCCAnalyses::MCParticle::get_p(HardWs_all)"
        )


        # peak center and width
        df = df.Define(
            "Ws_on_and_off_shell",
            "FCCAnalyses::ZHfunctions::get_on_and_off_shell_WW_160ecm(HardWs_all, Particle, 80.4, 2.1)"
        )

        df = df.Define("W_on_shell_idx",  "Ws_on_and_off_shell.high_mass_idx")
        df = df.Define("W_off_shell_idx", "Ws_on_and_off_shell.low_mass_idx")

         # quark_pdg = {
        #     'd': 1,
        #     'u': 2,
        #     's': 3,
        #     'c': 4,
        #     'b': 5,
        #     't': 6,
        # }


        # loop over all quarks pair decays from W and saves the events by object to "W_to{q1}{q2}"

        # start writing my python code to treemaker instead

        df = df.Alias("DaughterIndex", "Particle#1.index")

        # there are multiple copies of the W, and we look for the Ws with only status 22, but these copies
        # are not where the daughters are always stored
        # the daughters of status-22 can be other ws objects used in pythia8 that reference the same w,
        # some of which may store the quark daughters, thus we have to have a function that looks for the w with quark daughters
        
        df = df.Define("W_on_shell_decay_idx",
            "FCCAnalyses::ZHfunctions::get_decaying_W_idx(W_on_shell_idx, Particle, DaughterIndex)")
        df = df.Define("W_off_shell_decay_idx",
            "FCCAnalyses::ZHfunctions::get_decaying_W_idx(W_off_shell_idx, Particle, DaughterIndex)")

        dup_tracker = []

        for q1, id1 in quark_pdg.items():

            dup_tracker.append(q1)

            for q2, id2 in quark_pdg.items():

                if q2 in dup_tracker:
                    continue


                # goes from w idxs to w and q idxs to w and q particle objects
                df = df.Define(
                    f"W_on_shell_to_{q1}_{q2}_idxs",
                    f"FCCAnalyses::MCParticle::get_indices_MotherByIndex(W_on_shell_decay_idx, {{{id1},{id2}}}, false, true, true, Particle, DaughterIndex)"
                )

                df = df.Define(
                    f"W_off_shell_to_{q1}_{q2}_idxs",
                    f"FCCAnalyses::MCParticle::get_indices_MotherByIndex(W_off_shell_decay_idx, {{{id1},{id2}}}, false, true, true, Particle, DaughterIndex)"
                )


                # [W, q1, q2] as MCParticleData objects — truth data
                df = df.Define(
                    f"W_on_shell_to_{q1}_{q2}_objs",
                    f"W_on_shell_to_{q1}_{q2}_idxs.size() > 0 ? "
                    f"FCCAnalyses::ZHfunctions::get_mc(W_on_shell_to_{q1}_{q2}_idxs, Particle) : "
                    f"ROOT::VecOps::RVec<edm4hep::MCParticleData>{{}}"
                )
                df = df.Define(
                    f"W_off_shell_to_{q1}_{q2}_objs",
                    f"W_off_shell_to_{q1}_{q2}_idxs.size() > 0 ? "
                    f"FCCAnalyses::ZHfunctions::get_mc(W_off_shell_to_{q1}_{q2}_idxs, Particle) : "
                    f"ROOT::VecOps::RVec<edm4hep::MCParticleData>{{}}"
                )

                df = df.Define(
                    f"W_on_shell_to_{q1}_{q2}_quarks_idxs",
                    f"W_on_shell_to_{q1}_{q2}_idxs.size() > 2 ? "
                    f"ROOT::VecOps::Take(W_on_shell_to_{q1}_{q2}_idxs, {{1, 2}}) : "
                    f"ROOT::VecOps::RVec<int>{{}}"
                )

                df = df.Define(
                    f"W_off_shell_to_{q1}_{q2}_quarks_idxs",
                    f"W_off_shell_to_{q1}_{q2}_idxs.size() > 2 ? "
                    f"ROOT::VecOps::Take(W_off_shell_to_{q1}_{q2}_idxs, {{1, 2}}) : "
                    f"ROOT::VecOps::RVec<int>{{}}"
                )

                df = df.Define(
                    f"All_W_quarks_to_{q1}_{q2}_idx",
                    f"(W_on_shell_to_{q1}_{q2}_idxs.size() > 2) && (W_off_shell_to_{q1}_{q2}_idxs.size() > 2) ? "
                    f"ROOT::VecOps::Concatenate(W_on_shell_to_{q1}_{q2}_quarks_idxs, W_off_shell_to_{q1}_{q2}_quarks_idxs) : "
                    f"ROOT::VecOps::RVec<int>{{}}"
                )

                # df = df.Define(
                #     f"All_W_quarks_to_{q1}_{q2}_pairs_idx",
                #     f"(All_W_quarks_to_{q1}_{q2}_idx.size() == 4) ? "
                #     f"ROOT::VecOps::Combinations(All_W_quarks_to_{q1}_{q2}_idx, 2) : "
                #     f"ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>>{{}}"
                # )

                # df = df.Define(
                #     f"Mass_{q1}_{q2}_pairs",
                #     f"(All_W_quarks_to_{q1}_{q2}_pairs_idx.size() >= 2) ? "
                #     f"FCCAnalyses::ZHfunctions::get_pair_masses(All_W_quarks_to_{q1}_{q2}_idx, All_W_quarks_to_{q1}_{q2}_pairs_idx, Particle) : "
                #     f"ROOT::VecOps::RVec<double>{{}}"
                # )

                # df = df.Define(
                #     f"Candidate_On_Shell_W_{q1}_{q2}_pairs",
                #     f"(Mass_{q1}_{q2}_pairs.size() == 6) ? "
                #     f"FCCAnalyses::ZHfunctions::compare_pair_mass_to_w(Mass_{q1}_{q2}_pairs, All_W_quarks_to_{q1}_{q2}_idx, All_W_quarks_to_{q1}_{q2}_pairs_idx, Particle, W_on_shell_decay_idx) : "
                #     f"FCCAnalyses::ZHfunctions::OnOffidx{{}}"
                # )

                # df = df.Define(f"W_on_shell_qq_p_idxs_{q1}_{q2}",  f"Candidate_On_Shell_W_{q1}_{q2}_pairs.on_shell_idx")
                # df = df.Define(f"W_off_shell_qq_p_idxs_{q1}_{q2}", f"Candidate_On_Shell_W_{q1}_{q2}_pairs.off_shell_idx")

                


        ##### now want to define the two pairs invariant mass as well as if they are the true match

        on_shell_had  = " || ".join(f"W_on_shell_to_{q1}_{q2}_quarks_idxs.size() > 0"  for (q1, _), (q2, _) in combinations(quark_pdg.items(), 2))
        off_shell_had = " || ".join(f"W_off_shell_to_{q1}_{q2}_quarks_idxs.size() > 0" for (q1, _), (q2, _) in combinations(quark_pdg.items(), 2))
        df = df.Filter(f"({on_shell_had}) && ({off_shell_had})", "fully hadronic")


        # now we combine on and off shells quarks per event, 
        # so we filter across all names for the non-empty qq channel and combine

        channel_list = [(q1, q2) for (q1,_),(q2,_) in combinations(quark_pdg.items(), 2)]

        on_expr = "ROOT::VecOps::RVec<int>{}"
        for q1, q2 in reversed(channel_list):
            on_expr = f"(W_on_shell_to_{q1}_{q2}_quarks_idxs.size() > 0 ? W_on_shell_to_{q1}_{q2}_quarks_idxs : {on_expr})"

        off_expr = "ROOT::VecOps::RVec<int>{}"
        for q1, q2 in reversed(channel_list):
            off_expr = f"(W_off_shell_to_{q1}_{q2}_quarks_idxs.size() > 0 ? W_off_shell_to_{q1}_{q2}_quarks_idxs : {off_expr})"

        df = df.Define("on_shell_quark_idxs", on_expr)
        df = df.Define("off_shell_quark_idxs", off_expr)
        df = df.Define("All_W_quarks_idx", "ROOT::VecOps::Concatenate(on_shell_quark_idxs, off_shell_quark_idxs)")
        df = df.Define("on_shell_quark_objs",  "FCCAnalyses::ZHfunctions::get_mc(on_shell_quark_idxs, Particle)")
        df = df.Define("off_shell_quark_objs", "FCCAnalyses::ZHfunctions::get_mc(off_shell_quark_idxs, Particle)")

        # want the 4 vecs of all quarks in hadronic events
        df = df.Define("all_W_quarks_obj", "FCCAnalyses::ZHfunctions::get_mc(All_W_quarks_idx, Particle)")
        df = df.Define("all_W_quarks_tlv", "FCCAnalyses::MCParticle::get_tlv(all_W_quarks_obj)")

        
        df = df.Define(
            f"All_W_quarks_pairs_idx",
            f"(All_W_quarks_idx.size() == 4) ? "
            f"ROOT::VecOps::Combinations(All_W_quarks_idx, 2) : "
            f"ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>>{{}}"
        )

        df = df.Define(
            f"Mass_qq_pairs",
            f"(All_W_quarks_pairs_idx.size() >= 2) ? "
            f"FCCAnalyses::ZHfunctions::get_pair_masses(All_W_quarks_idx, All_W_quarks_pairs_idx, Particle) : "
            f"ROOT::VecOps::RVec<double>{{}}"
        )

        df = df.Define(
            f"Candidate_W_qq_pairs",
            f"(Mass_qq_pairs.size() == 6) ? "
            f"FCCAnalyses::ZHfunctions::compare_pair_mass_to_w(Mass_qq_pairs, All_W_quarks_idx, All_W_quarks_pairs_idx, Particle, W_on_shell_decay_idx) : "
            f"FCCAnalyses::ZHfunctions::OnOffidx{{}}"
        )

        df = df.Define(f"Candidate_on_shell_W_qq_p_idxs",  f"Candidate_W_qq_pairs.on_shell_idx")
        df = df.Define(f"Candidate_off_shell_W_qq_p_idxs", f"Candidate_W_qq_pairs.off_shell_idx")
        df = df.Define(f"Candidate_on_shell_W_qq_mass", f"Candidate_W_qq_pairs.on_shell_mass")
        df = df.Define(f"Candidate_off_shell_W_qq_mass", f"Candidate_W_qq_pairs.off_shell_mass")
        df = df.Define(f"Candidate_on_shell_W_qq_p_flavor", f"Candidate_W_qq_pairs.on_shell_flavor")
        df = df.Define(f"Candidate_off_shell_W_qq_p_flavor", f"Candidate_W_qq_pairs.off_shell_flavor")

        df = df.Define(f"W_qq_match_truth", f"Candidate_W_qq_pairs.match_truth")


        for lep, ids in lepton_pdg.items():
            lep_id, nu_id = ids

            df = df.Define(
                    f"W_on_shell_to_{lep}_nu",
                    f"FCCAnalyses::ZHfunctions::get_mc(FCCAnalyses::MCParticle::get_indices_MotherByIndex(W_on_shell_idx, {{{lep_id},{nu_id}}}, false, true, false, Particle, DaughterIndex), Particle)"
                )

            df = df.Define(
                    f"W_off_shell_to_{lep}_nu",
                    f"FCCAnalyses::ZHfunctions::get_mc(FCCAnalyses::MCParticle::get_indices_MotherByIndex(W_off_shell_idx, {{{lep_id},{nu_id}}}, false, true, false, Particle, DaughterIndex), Particle)"
                )

            df = df.Define(
                    f"W_on_shell_to_{lep}_nu_tlv",
                    f"FCCAnalyses::MCParticle::get_tlv(W_on_shell_to_{lep}_nu)"
                )

            df = df.Define(
                    f"W_off_shell_to_{lep}_nu_tlv",
                    f"FCCAnalyses::MCParticle::get_tlv(W_off_shell_to_{lep}_nu)"
                )





        nJets = 4  # WW → qqqq fully hadronic: 4 quarks → 4 jets

        collections_noleps = copy.deepcopy(collections)
        collections_noleps["PFParticles"] = "ReconstructedParticlesNoMuNoEl"
        if saveExclJets:
            jetClusteringHelper = ExclusiveJetClusteringHelper(
                collections_noleps["PFParticles"], nJets
            )
        
        jetClusteringHelper_R5  = InclusiveJetClusteringHelper(
            collections_noleps["PFParticles"], 0.5, 10, "R5"
        )
        if  saveExclJets:
            df = jetClusteringHelper.define(df)
        df = jetClusteringHelper_R5.define(df)
        ## define jet flavour tagging parameters
        if saveExclJets:
            jetFlavourHelper = JetFlavourHelper(
                collections_noleps,
                jetClusteringHelper.jets,
                jetClusteringHelper.constituents,
            )
        jetFlavourHelper_R5 = JetFlavourHelper(
            collections_noleps,
            jetClusteringHelper_R5.jets,
            jetClusteringHelper_R5.constituents,
            "R5",
        )
        if saveExclJets:    df = jetFlavourHelper.define(df)
        
        df = jetFlavourHelper_R5.define(df)
        ## tagger inference
        if  saveExclJets: df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)
        df = jetFlavourHelper_R5.inference(weaver_preproc, weaver_model,df)

        df = df.Define(
            "lep_p", "muons_sel_iso.size() >0 ? FCCAnalyses::ReconstructedParticle::get_p(muons_sel_iso)[0] : (electrons_sel_iso.size() > 0 ? FCCAnalyses::ReconstructedParticle::get_p(electrons_sel_iso)[0] : -999) "
        )
        df = df.Define(
            'lep_theta', 'muons_sel_iso.size() >0 ? FCCAnalyses::ReconstructedParticle::get_theta(muons_sel_iso)[0] : (electrons_sel_iso.size() > 0 ? FCCAnalyses::ReconstructedParticle::get_theta(electrons_sel_iso)[0] : -999) '
        )
        df = df.Define(
            'lep_phi', 'muons_sel_iso.size() >0 ? FCCAnalyses::ReconstructedParticle::get_phi(muons_sel_iso)[0] : (electrons_sel_iso.size() > 0 ? FCCAnalyses::ReconstructedParticle::get_phi(electrons_sel_iso)[0] : -999) '
        )
        df = df.Define(
            "nlep",
            "electrons_sel_iso.size()+muons_sel_iso.size()")
        


        df = df.Define(
            "missing_p",
            "FCCAnalyses::ReconstructedParticle::get_p(MissingET)[0]",
        )

        
        df = df.Define(
            'missing_p_theta', 'ReconstructedParticle::get_theta(MissingET)[0]',
        )

        df = df.Define(
            'missing_p_phi', 'ReconstructedParticle::get_phi(MissingET)[0]',
        )


        # add transverse momentum

        df = df.Define(
            'missing_pt', 'ReconstructedParticle::get_pt(MissingET)[0]',
        )

        if    saveExclJets:
            df = df.Define(
                "jets_p4",
                "JetConstituentsUtils::compute_tlv_jets({})".format(
                    jetClusteringHelper.jets
                ),
            )

        # start analyzing jets by first matching jets to quarks, so need 4vecs of quarks and jets
        # returns indexs that map from jets to quarks
        df = df.Define(
            'matched_jets_to_q',
            'FCCAnalyses::ZHfunctions::MatchJetsToQuarks(jets_p4, all_W_quarks_tlv, 0.1)'
        )


        df = df.Define(
            'matched_jets_to_q_idx',
            'matched_jets_to_q.idx'
        )

        df = df.Define(
            'matched_jets_to_q_under_min_delR',
            'matched_jets_to_q.under_min_delR'
        )

        df = df.Define(
            'matched_jets_to_q_delta_Rs',
            'matched_jets_to_q.delta_Rs'
        )

        df = df.Define(
            'matched_jets_to_q_delta_etas',
            'matched_jets_to_q.delta_etas'
        )

        df = df.Define(
            'matched_jets_to_q_delta_phis',
            'matched_jets_to_q.delta_phis'
        )

        # make seperate collections for each jet
        for i in range(nJets):
            df = df.Define(f"simple_jet_{i+1}_deltaR", f"matched_jets_to_q_delta_Rs[{i}]")
            df = df.Define(f"simple_jet_{i+1}_delta_phi", f"matched_jets_to_q_delta_phis[{i}]")
            df = df.Define(f"simple_jet_{i+1}_delta_eta", f"matched_jets_to_q_delta_etas[{i}]")



        # #matching jets to quarks using a chi-squared fit using only deltaR
        # df = df.Define(
        #     'chi2_matched_jets_to_q',
        #     'FCCAnalyses::ZHfunctions::JtoQ_ChiSquared_deltaR(jets_p4, all_W_quarks_tlv, ROOT::VecOps::RVec<double>({0.05, 0.05, 0.05, 0.05}), 0.1)'
        # )

        # matching jets to quarks using a chi-squared fit using dPhi and dEta
        # sigmas defined in order of {eta, phi}
        # sigmas estimated from 5,000 events when using old deltaR algorithim (no chi-squared fit), {s_eta = 0.2480, s_phi = 0.2389}

        # greedy variant — no chi2 matching

        



        # Filter out events below dR = 0.1
        ##df = df.Filter("chi2_matched_jets_to_q_under_min_delR == 1", "all jets matched within delR for new matching")

        # for i in range(nJets):
        #     df = df.Define(f"chi2_jet_{i+1}_deltaR", f"chi2_matched_jets_to_q_delta_Rs[{i}]")
        



        # now we want a function that will do the chi squared test and return the best permutation




        #df = df.Filter("matched_jets_to_q_under_min_delR == 1", "all jets matched within delR")

        # now we want to take our 4vecs of our jets and compute all dijet masses
        # function currently only gives the masses, not how they are paired

        df = df.Define(
            "dijet_info",
            "FCCAnalyses::ZHfunctions::all_invariant_masses_and_pair_idxs(jets_p4)"
        )

        df = df.Define(
            "dijet_masses",
            "dijet_info.masses"
        )

        # list of tuples that say how the pairs are built from the jets vector
        df = df.Define(
            "dijet_pairs_idxs",
            "dijet_info.pair_idxs"
        )

        # temporary
        df = df.Define("dijet_pair_idx_a", "dijet_info.pair_idxs[0]")  # first jet of each pair
        df = df.Define("dijet_pair_idx_b", "dijet_info.pair_idxs[1]")  # second jet of each pair


        # now try to do a similar thing to what we did with quarks, may have to manipulate some indexs around since function expects only quark info

        
        # need to get dijet_pairs_idxs to point

        df = df.Define(
            "dijet_pairs_as_quark_idx",
            "FCCAnalyses::ZHfunctions::transform_pair_idxs(dijet_pairs_idxs, matched_jets_to_q_idx)"
        )


        df = df.Define(
            f"Candidate_reco_W_jj_pairs",
            f"(dijet_masses.size() == 6) ? "  #dijet mass   # indexs of quarks (so need to point jet_to_quark_idxs to get their idxs first)    # modified dijet_pairs_idxs (need to check these should be pointing to jets which need to point to quarks)
            f"FCCAnalyses::ZHfunctions::compare_pair_mass_to_w(dijet_masses, All_W_quarks_idx, dijet_pairs_as_quark_idx, Particle, W_on_shell_decay_idx) : "
            f"FCCAnalyses::ZHfunctions::OnOffidx{{}}"
        )

        df = df.Define("Candidate_reco_on_shell_W_jj_p_idxs",  "Candidate_reco_W_jj_pairs.on_shell_idx")
        df = df.Define("Candidate_reco_off_shell_W_jj_p_idxs", "Candidate_reco_W_jj_pairs.off_shell_idx")
        df = df.Define("Candidate_reco_on_shell_W_jj_mass",    "Candidate_reco_W_jj_pairs.on_shell_mass")
        df = df.Define("Candidate_reco_off_shell_W_jj_mass",   "Candidate_reco_W_jj_pairs.off_shell_mass")
        df = df.Define("Candidate_reco_on_shell_W_jj_flavor",  "Candidate_reco_W_jj_pairs.on_shell_flavor")
        df = df.Define("Candidate_reco_off_shell_W_jj_flavor", "Candidate_reco_W_jj_pairs.off_shell_flavor")
        df = df.Define("reco_W_jj_match_truth",                "Candidate_reco_W_jj_pairs.match_truth")

        # greedy variant — no chi2 scores to split











        # exclusive at 4 jets
        
        df = df.Define(
            "jets_R5_p4",
            "JetConstituentsUtils::compute_tlv_jets({})".format(
                jetClusteringHelper_R5.jets
            ),
        )

        df = df.Define("jets_R5_p",           "JetClusteringUtils::get_p({})".format(jetClusteringHelper_R5.jets))
        df = df.Define("jets_R5_theta",       "JetClusteringUtils::get_theta({})".format(jetClusteringHelper_R5.jets))


        # df = df.Define("jet1_R5_p","jets_R5_p[0]")
        # df = df.Define("jet2_R5_p","jets_R5_p[1] ? jets_R5_p[2] : -999")
        # df = df.Define("jet3_R5_p","jets_R5_p.size()>2 ? jets_R5_p[2] : -999")
        # df = df.Define("jet4_R5_p","jets_R5_p.size()>3 ? jets_R5_p[3] : -999")
        # df = df.Define("jet5_R5_p","jets_R5_p.size()>4 ? jets_R5_p[4] : -999")
        # df = df.Define("jet6_R5_p","jets_R5_p.size()>5 ? jets_R5_p[5] : -999")

        # df = df.Define("jet1_R5_theta","jets_R5_theta[0]")
        # df = df.Define("jet2_R5_theta","jets_R5_theta[1]")
        # df = df.Define("jet3_R5_theta","jets_R5_theta.size()>2 ? jets_R5_theta[2] : -999")
        # df = df.Define("jet4_R5_theta","jets_R5_theta.size()>3 ? jets_R5_theta[3] : -999")
        # df = df.Define("jet5_R5_theta","jets_R5_theta.size()>4 ? jets_R5_theta[4] : -999")
        # df = df.Define("jet6_R5_theta","jets_R5_theta.size()>5 ? jets_R5_theta[5] : -999")
        
        df = df.Define("jets_R5_pflavor", "JetTaggingUtils::get_flavour({}, Particle)".format(jetClusteringHelper_R5.jets) )
        # df = df.Define("jet1_R5_pflavor","jets_R5_pflavor[0]")
        # df = df.Define("jet2_R5_pflavor","jets_R5_pflavor[1]")
        # df = df.Define("jet3_R5_pflavor","jets_R5_p.size()>2 ? jets_R5_pflavor[2] : -999")
        # df = df.Define("jet4_R5_pflavor","jets_R5_p.size()>3 ? jets_R5_pflavor[3] : -999")
        # df = df.Define("jet5_R5_pflavor","jets_R5_p.size()>4 ? jets_R5_pflavor[4] : -999")
        # df = df.Define("jet6_R5_pflavor","jets_R5_p.size()>5 ? jets_R5_pflavor[5] : -999")
        df = df.Define("njets_R5",       "return int(jets_R5_pflavor.size())")


        

        
        
        df = df.Define("jets_R5_btag_true", "JetTaggingUtils::get_btag({}, 1.0)".format('jets_R5_pflavor'))
        df = df.Define("jets_R5_ctag_true", "JetTaggingUtils::get_ctag({}, 1.0)".format('jets_R5_pflavor'))
        df = df.Define("jets_R5_ltag_true", "JetTaggingUtils::get_ltag({}, 1.0)".format('jets_R5_pflavor'))
        df = df.Define("jets_R5_gtag_true", "JetTaggingUtils::get_gtag({}, 1.0)".format('jets_R5_pflavor'))
        

        b_eff = .95
        c_eff = 10**-1.5
        l_eff = 10**-3
        g_eff = 10**-1.7
        uncert_b_eff = 0.01
        
        # df = df.Define("jets_R5_btag_eff_p9",  "JetTaggingUtils::get_btag({},{},{},{},{})".format('jets_R5_pflavor', b_eff, c_eff, l_eff, g_eff))
        # df = df.Define("jets_R5_btag_eff_p89", "JetTaggingUtils::get_btag({},{},{},{},{})".format('jets_R5_pflavor', b_eff-uncert_b_eff, c_eff, l_eff, g_eff))
        # df = df.Define("jets_R5_btag_eff_p91", "JetTaggingUtils::get_btag({},{},{},{},{})".format('jets_R5_pflavor', b_eff+uncert_b_eff, c_eff, l_eff, g_eff))
        
        df = df.Define("jets_R5_btagged_true", "JetTaggingUtils::sel_tag(true)(jets_R5_btag_true,{})".format(jetClusteringHelper_R5.jets))
        df = df.Define("jets_R5_ctagged_true", "JetTaggingUtils::sel_tag(true)(jets_R5_ctag_true,{})".format(jetClusteringHelper_R5.jets))
        df = df.Define("jets_R5_ltagged_true", "JetTaggingUtils::sel_tag(true)(jets_R5_ltag_true,{})".format(jetClusteringHelper_R5.jets))
        df = df.Define("jets_R5_gtagged_true", "JetTaggingUtils::sel_tag(true)(jets_R5_gtag_true,{})".format(jetClusteringHelper_R5.jets))

        # df = df.Define("jets_R5_btagged_eff_p9", "JetTaggingUtils::sel_tag(true)(jets_R5_btag_eff_p9,{})".format(jetClusteringHelper_R5.jets))
        # df = df.Define("jets_R5_btagged_eff_p89", "JetTaggingUtils::sel_tag(true)(jets_R5_btag_eff_p89,{})".format(jetClusteringHelper_R5.jets))
        # df = df.Define("jets_R5_btagged_eff_p91", "JetTaggingUtils::sel_tag(true)(jets_R5_btag_eff_p91,{})".format(jetClusteringHelper_R5.jets))

        df = df.Define("nbjets_R5_true", "return int(jets_R5_btagged_true.size())")
        df = df.Define("ncjets_R5_true", "return int(jets_R5_ctag_true.size())")
        df = df.Define("nljets_R5_true", "return int(jets_R5_ltag_true.size())")
        df = df.Define("ngjets_R5_true", "return int(jets_R5_gtag_true.size())")


        df = df.Define("bjets_R5_true",  "JetConstituentsUtils::compute_tlv_jets({})". format('jets_R5_btagged_true'))
        df = df.Define("ljets_R5_true",  "JetConstituentsUtils::compute_tlv_jets({})". format('jets_R5_ltagged_true'))
        df = df.Define("bjet1_R5_true_p","nbjets_R5_true > 0 ? bjets_R5_true[0].P() : -999")
        df = df.Define("ljet1_R5_true_p","nljets_R5_true > 0 ? ljets_R5_true[0].P() : -999")

        

        
        # df = df.Define("nbjets_R5_eff_p9", "return int(jets_R5_btagged_eff_p9.size())")
        # df = df.Define("nbjets_R5_eff_p89", "return int(jets_R5_btagged_eff_p89.size())")
        # df = df.Define("nbjets_R5_eff_p91", "return int(jets_R5_btagged_eff_p91.size())")

        # df = df.Define("bjet_R5_eff_p9_p4",  "JetConstituentsUtils::compute_tlv_jets({})". format('jets_R5_btagged_eff_p9'))
        # df = df.Define("bjet_R5_eff_p91_p4", "JetConstituentsUtils::compute_tlv_jets({})".format('jets_R5_btagged_eff_p91'))
        # df = df.Define("bjet_R5_eff_p89_p4", "JetConstituentsUtils::compute_tlv_jets({})".format('jets_R5_btagged_eff_p89'))
        # df = df.Define("mbbar_p9",  "nbjets_R5_eff_p9 >   1 ? JetConstituentsUtils::InvariantMass(bjet_R5_eff_p9_p4[0],  bjet_R5_eff_p9_p4[1])  : -999")
        # df = df.Define("mbbar_p91", "nbjets_R5_eff_p91 >  1 ? JetConstituentsUtils::InvariantMass(bjet_R5_eff_p91_p4[0], bjet_R5_eff_p91_p4[1]) : -999")
        # df = df.Define("mbbar_p89", "nbjets_R5_eff_p89 >  1 ? JetConstituentsUtils::InvariantMass(bjet_R5_eff_p89_p4[0], bjet_R5_eff_p89_p4[1]) : -999")


        # df = df.Define("jet1_R5_isG", "recojet_isG_R5[0]")
        # df = df.Define("jet2_R5_isG", "recojet_isG_R5[1]")
        # df = df.Define("jet3_R5_isG", "jets_R5_p.size()>2 ? recojet_isG_R5[2] : -999")
        # df = df.Define("jet4_R5_isG", "jets_R5_p.size()>3 ? recojet_isG_R5[3] : -999")
        # df = df.Define("jet5_R5_isG", "jets_R5_p.size()>4 ? recojet_isG_R5[4] : -999")
        # df = df.Define("jet6_R5_isG", "jets_R5_p.size()>5 ? recojet_isG_R5[5] : -999")

        # df = df.Define("jet1_R5_isU", "recojet_isU_R5[0]")
        # df = df.Define("jet2_R5_isU", "recojet_isU_R5[1]")
        # df = df.Define("jet3_R5_isU", "jets_R5_p.size()>2 ? recojet_isU_R5[2] : -999")
        # df = df.Define("jet4_R5_isU", "jets_R5_p.size()>3 ? recojet_isU_R5[3] : -999")
        # df = df.Define("jet5_R5_isU", "jets_R5_p.size()>4 ? recojet_isU_R5[4] : -999")
        # df = df.Define("jet6_R5_isU", "jets_R5_p.size()>5 ? recojet_isU_R5[5] : -999")

        # df = df.Define("jet1_R5_isB", "recojet_isB_R5[0]")
        # df = df.Define("jet2_R5_isB", "recojet_isB_R5[1]")
        # df = df.Define("jet3_R5_isB", "jets_R5_p.size()>2 ? recojet_isB_R5[2] : -999")
        # df = df.Define("jet4_R5_isB", "jets_R5_p.size()>3 ? recojet_isB_R5[3] : -999")
        # df = df.Define("jet5_R5_isB", "jets_R5_p.size()>4 ? recojet_isB_R5[4] : -999")
        # df = df.Define("jet6_R5_isB", "jets_R5_p.size()>5 ? recojet_isB_R5[5] : -999")
        
        # df = df.Define("jet1_R5_isS", "recojet_isS_R5[0]")
        # df = df.Define("jet2_R5_isS", "recojet_isS_R5[1]")
        # df = df.Define("jet3_R5_isS", "jets_R5_p.size()>2 ? recojet_isS_R5[2] : -999")
        # df = df.Define("jet4_R5_isS", "jets_R5_p.size()>3 ? recojet_isS_R5[3] : -999")
        # df = df.Define("jet5_R5_isS", "jets_R5_p.size()>4 ? recojet_isS_R5[4] : -999")
        # df = df.Define("jet6_R5_isS", "jets_R5_p.size()>5 ? recojet_isS_R5[5] : -999")

        # df = df.Define("jet1_R5_isC", "recojet_isC_R5[0]")
        # df = df.Define("jet2_R5_isC", "recojet_isC_R5[1]")
        # df = df.Define("jet3_R5_isC", "jets_R5_p.size()>2 ? recojet_isC_R5[2] : -999")
        # df = df.Define("jet4_R5_isC", "jets_R5_p.size()>3 ? recojet_isC_R5[3] : -999")
        # df = df.Define("jet5_R5_isC", "jets_R5_p.size()>4 ? recojet_isC_R5[4] : -999")
        # df = df.Define("jet6_R5_isC", "jets_R5_p.size()>5 ? recojet_isC_R5[5] : -999")

        # df = df.Define("jet1_R5_isD", "recojet_isD_R5[0]")
        # df = df.Define("jet2_R5_isD", "recojet_isD_R5[1]")
        # df = df.Define("jet3_R5_isD", "jets_R5_p.size()>2 ? recojet_isD_R5[2] : -999")
        # df = df.Define("jet4_R5_isD", "jets_R5_p.size()>3 ? recojet_isD_R5[3] : -999")
        # df = df.Define("jet5_R5_isD", "jets_R5_p.size()>4 ? recojet_isD_R5[4] : -999")
        # df = df.Define("jet6_R5_isD", "jets_R5_p.size()>5 ? recojet_isD_R5[5] : -999")

        # df = df.Define("jet1_R5_isTAU", "recojet_isTAU_R5[0]")
        # df = df.Define("jet2_R5_isTAU", "recojet_isTAU_R5[1]")
        # df = df.Define("jet3_R5_isTAU", "jets_R5_p.size()>2 ? recojet_isTAU_R5[2] : -999")
        # df = df.Define("jet4_R5_isTAU", "jets_R5_p.size()>3 ? recojet_isTAU_R5[3] : -999")
        # df = df.Define("jet5_R5_isTAU", "jets_R5_p.size()>4 ? recojet_isTAU_R5[4] : -999")
        # df = df.Define("jet6_R5_isTAU", "jets_R5_p.size()>5 ? recojet_isTAU_R5[5] : -999")

        df = df.Define("jets_R5_isB","recojet_isB_R5")
        
        # ADD BACK IN TO COLLECTIONS
        df = df.Define("bjets_R5_WPp5","ZHfunctions::sel_btag(0.5)(jets_R5_isB)")
        df = df.Define("bjets_R5_WPp8","ZHfunctions::sel_btag(0.8)(jets_R5_isB)")
        df = df.Define("bjets_R5_WPp85","ZHfunctions::sel_btag(0.85)(jets_R5_isB)")
        df = df.Define("bjets_R5_WPp9","ZHfunctions::sel_btag(0.9)(jets_R5_isB)")

        df = df.Define("nbjets_R5_WPp5","bjets_R5_WPp5.size()")
        df = df.Define("nbjets_R5_WPp8","bjets_R5_WPp8.size()")
        df = df.Define("nbjets_R5_WPp85","bjets_R5_WPp85.size()")
        df = df.Define("nbjets_R5_WPp9","bjets_R5_WPp9.size()")

        
        #['recojet_isG_R5', 'recojet_isU_R5', 'recojet_isS_R5', 'recojet_isC_R5', 'recojet_isB_R5', 'recojet_isTAU_R5', 'recojet_isD_R5', 'jet_nmu_R5', 'jet_nel_R5', 'jet_nchad_R5', 'jet_ngamma_R5', 'jet_nnhad_R5']
        if  saveExclJets:
            df = df.Define("jets_isB",   "JetFlavourUtils::get_weight(MVAVec_, 4)")
            df = df.Define("bjets_WPp5", "ZHfunctions::sel_btag(0.5)(jets_isB)")
            df = df.Define("bjets_WPp8", "ZHfunctions::sel_btag(0.8)(jets_isB)")
            df = df.Define("bjets_WPp85", "ZHfunctions::sel_btag(0.85)(jets_isB)")
            df = df.Define("bjets_WPp9", "ZHfunctions::sel_btag(0.9)(jets_isB)")
            df = df.Define("nbjets_WPp5", "bjets_WPp5.size()")
            df = df.Define("nbjets_WPp8", "bjets_WPp8.size()")
            df = df.Define("nbjets_WPp85", "bjets_WPp85.size()")
            df = df.Define("nbjets_WPp9", "bjets_WPp9.size()")
    
            # df = df.Define("jet1", "jets_p4[0]")
            # df = df.Define("jet2", "jets_p4[1]")
            # df = df.Define("jet3", "jets_p4[2]")
            # df = df.Define("jet4", "jets_p4[3]")
            # df = df.Define("jet1_p","jet1.P()")
            # df = df.Define("jet2_p","jet2.P()")
            # df = df.Define("jet3_p","jet3.P()")
            # df = df.Define("jet4_p","jet4.P()")
            # df = df.Define("jet5_p","jets_p4.size()>4 ? jets_p4[4].P() : -999")
            # df = df.Define("jet6_p","jets_p4.size()>5 ? jets_p4[5].P() : -999")
            df = df.Define("recojet_theta", "JetClusteringUtils::get_theta(jet)")
            # df = df.Define("jet1_theta","recojet_theta[0]")
            # df = df.Define("jet2_theta","recojet_theta[1]")
            # df = df.Define("jet3_theta","recojet_theta[2]")
            # df = df.Define("jet4_theta","recojet_theta[3]")
            # df = df.Define("jet5_theta","jets_p4.size()>4 ? recojet_theta[4] : -999")
            # df = df.Define("jet6_theta","jets_p4.size()>5 ? recojet_theta[5] : -999")
            
            df = df.Define("jet1_isTau","recojet_isTAU[0]")
            df = df.Define("jet2_isTau","recojet_isTAU[1]")
            df = df.Define("jet3_isTau","recojet_isTAU[2]")
            df = df.Define("jet4_isTau","recojet_isTAU[3]")
            df = df.Define("jet5_isTau","jets_p4.size()>4 ? recojet_isTAU[4] : -999")
            df = df.Define("jet6_isTau","jets_p4.size()>5 ? recojet_isTAU[5] : -999")
        
            df = df.Define("recojet_phi", "JetClusteringUtils::get_phi_std(jet)")
            df = df.Define("jet1_phi","recojet_phi[0]")
            df = df.Define("jet2_phi","recojet_phi[1]")
            df = df.Define("jet3_phi","recojet_phi[2]")
            df = df.Define("jet4_phi","recojet_phi[3]")
            df = df.Define("jet5_phi","jets_p4.size()>4 ? recojet_phi[4] : -999")
            df = df.Define("jet6_phi","jets_p4.size()>5 ? recojet_phi[5] : -999")
            df = df.Define("njets", "jets_p4.size()")
        
            df = df.Define("d_12", "JetClusteringUtils::get_exclusive_dmerge(_jet, 1)")
            df = df.Define("d_23", "JetClusteringUtils::get_exclusive_dmerge(_jet, 2)")
            df = df.Define("d_34", "JetClusteringUtils::get_exclusive_dmerge(_jet, 3)")
            df = df.Define("d_45", "jets_p4.size()>4 ? JetClusteringUtils::get_exclusive_dmerge(_jet, 4) : -999")
            df = df.Define("d_56", "jets_p4.size()>5 ? JetClusteringUtils::get_exclusive_dmerge(_jet, 5) : -999")


            df = df.Define("jet1_isG", "JetFlavourUtils::get_weight(MVAVec_, 0)[0]")
            df = df.Define("jet2_isG", "JetFlavourUtils::get_weight(MVAVec_, 0)[1]")
            df = df.Define("jet3_isG", "JetFlavourUtils::get_weight(MVAVec_, 0)[2]")
            df = df.Define("jet4_isG", "JetFlavourUtils::get_weight(MVAVec_, 0)[3]")
            df = df.Define("jet5_isG", "jets_p4.size()>4 ? JetFlavourUtils::get_weight(MVAVec_, 0)[4] : -999")
            df = df.Define("jet6_isG", "jets_p4.size()>5 ? JetFlavourUtils::get_weight(MVAVec_, 0)[5] : -999")
            
            df = df.Define("jet1_isQ", "JetFlavourUtils::get_weight(MVAVec_, 1)[0]")
            df = df.Define("jet2_isQ", "JetFlavourUtils::get_weight(MVAVec_, 1)[1]")
            df = df.Define("jet3_isQ", "JetFlavourUtils::get_weight(MVAVec_, 1)[2]")
            df = df.Define("jet4_isQ", "JetFlavourUtils::get_weight(MVAVec_, 1)[3]")
            df = df.Define("jet5_isQ", "jets_p4.size()>4 ? JetFlavourUtils::get_weight(MVAVec_, 1)[4] : -999")
            df = df.Define("jet6_isQ", "jets_p4.size()>5 ? JetFlavourUtils::get_weight(MVAVec_, 1)[5] : -999")
            
            df = df.Define("jet1_isS", "JetFlavourUtils::get_weight(MVAVec_, 2)[0]")
            df = df.Define("jet2_isS", "JetFlavourUtils::get_weight(MVAVec_, 2)[1]")
            df = df.Define("jet3_isS", "JetFlavourUtils::get_weight(MVAVec_, 2)[2]")
            df = df.Define("jet4_isS", "JetFlavourUtils::get_weight(MVAVec_, 2)[3]")
            df = df.Define("jet5_isS", "jets_p4.size()>4 ? JetFlavourUtils::get_weight(MVAVec_, 2)[4] : -999")
            df = df.Define("jet6_isS", "jets_p4.size()>5 ? JetFlavourUtils::get_weight(MVAVec_, 2)[5] : -999")
            
            df = df.Define("jet1_isC", "JetFlavourUtils::get_weight(MVAVec_, 3)[0]")
            df = df.Define("jet2_isC", "JetFlavourUtils::get_weight(MVAVec_, 3)[1]")
            df = df.Define("jet3_isC", "JetFlavourUtils::get_weight(MVAVec_, 3)[2]")
            df = df.Define("jet4_isC", "JetFlavourUtils::get_weight(MVAVec_, 3)[3]")
            df = df.Define("jet5_isC", "jets_p4.size()>4 ? JetFlavourUtils::get_weight(MVAVec_, 3)[4] : -999")
            df = df.Define("jet6_isC", "jets_p4.size()>5 ? JetFlavourUtils::get_weight(MVAVec_, 3)[5] : -999")
    
            df = df.Define("jet1_isB", "JetFlavourUtils::get_weight(MVAVec_, 4)[0]")
            df = df.Define("jet2_isB", "JetFlavourUtils::get_weight(MVAVec_, 4)[1]")
            df = df.Define("jet3_isB", "JetFlavourUtils::get_weight(MVAVec_, 4)[2]")
            df = df.Define("jet4_isB", "JetFlavourUtils::get_weight(MVAVec_, 4)[3]")
            df = df.Define("jet5_isB", "jets_p4.size()>4 ? JetFlavourUtils::get_weight(MVAVec_, 4)[4] : -999")
            df = df.Define("jet6_isB", "jets_p4.size()>5 ? JetFlavourUtils::get_weight(MVAVec_, 4)[5] : -999")


        return df

    # __________________________________________________________
    # Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        #print('incl jets',jetFlavourHelper_R5.outputBranches())
        #print('excl jets',jetFlavourHelper.outputBranches())
        #all_branches += jetFlavourHelper_R5.outputBranches()
        return all_branches

        
        #all_branches+= jetClusteringHelper.outputBranches()

        ## outputs jet scores and constituent breakdown
        #branchList += jetFlavourHelper.outputBranches()
    
        #return all_branches #branchList




        ##test command fccanalysis run --nevents=10 treemaker_WbWb_reco.py
