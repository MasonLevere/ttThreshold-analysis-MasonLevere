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

# ww sample just under threshold
# One on shell, one off shell WW events
    "p8_ee_WW_ecm160": {
        "fraction": 1,
    },


### ONES WE WILL USE SOON
# # semileptonic WW -> mu nu qq, split by CKM matrix element
#     "wzp6_ee_munumuqq_Vcd_ecm163": {"fraction": 1},
#     "wzp6_ee_munumuqq_Vcb_ecm163": {"fraction": 1},
#     "wzp6_ee_munumuqq_Vub_ecm163": {"fraction": 1},
#     "wzp6_ee_munumuqq_Vcs_ecm163": {"fraction": 1},
#     "wzp6_ee_munumuqq_Vud_ecm163": {"fraction": 1},
#     "wzp6_ee_munumuqq_Vus_ecm163": {"fraction": 1},



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
outputDir   = "/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/{}_WW_new_matching".format(channel)


# additional/costom C++ functions, defined in header files (optional)
#                   namespace is ZH         namespace is WWFunctions
includePaths = ["WWfunctions/functions.h", "WWfunctions/WWFunctions.h", "WWfunctions/BWPairing.h"]

## latest particle transformer model, trained on 9M jets in winter2023 samples

## model files locally stored on /eos
eos_dir ="/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
# new model file directory



# model_dir = (
#     "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_7classes_12_04_2023/"    # "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
# )

model_name = "fccee_flavtagging_edm4hep_wc" #"fccee_flavtagging_edm4hep_wc_v1"



model_dir = (
    "/afs/cern.ch/user/m/mlevere/private/FCCTutorial/ttThreshold-analysis/local_tagger"
)

model_name = "ZmWHmW_163_parT"

local_preproc = "{}/{}.json".format(model_dir, model_name)
local_model = "{}/{}.onnx".format(model_dir, model_name)

url_model_dir = "https://fccsw.web.cern.ch/fccsw/testsamples/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
url_preproc = "{}/{}.json".format(url_model_dir, model_name)
url_model = "{}/{}.onnx".format(url_model_dir, model_name)

## get local file, else download from url
def get_file_path(url, filename):
    if os.path.exists(filename):
        print('local filepath worked')
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

BW_BOSON   = os.environ.get("BW_BOSON", "W").strip().upper()
BW_SQRTD45 = float(os.environ.get("BW_SQRTD45_MAX", "7.0"))


_boson_pdg  = {"W": 24, "Z": 23}[BW_BOSON]



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


# --- sigma values for chi2 jet-quark matching ---
# Update these each iteration and rerun; fit_sigmas.py reads the output and reports new values.
SIGMA_ETA = 0.0505
SIGMA_PHI = 0.0529

SIGMA_W_ON_SHELL = 5.3665

# refers to detector resolution for gauss in voigt distribution
sigma = 3.6110261681321907


chi2_cut = 100



all_branches = [

    # --- gen truth ---
    "All_W_quarks_idx", "all_W_quarks_obj",
    "HardWs_all_mass",
    "W_on_shell_mass", "W_off_shell_mass",

    # --- jets ---
    "jets_p4",
    "n_reco_jets",

    # --- chi2 jet-quark matching ---
    'chi2_etaphi_best_chi2', 'chi2_etaphi_second_best_chi2', 'chi2_etaphi_third_best_chi2',
    'chi2_etaphi_delta_Rs', 'chi2_etaphi_delta_etas', 'chi2_etaphi_delta_phis',
    'chi2_dR_best_chi2', 'chi2_dR_delta_Rs',

    # --- truth matching ---
    "gen_pairing_true",
    "jet1_wlab", "jet2_wlab", "jet3_wlab", "jet4_wlab",
    "jet1_matched_q_dR", "jet2_matched_q_dR", "jet3_matched_q_dR", "jet4_matched_q_dR",
    "gen_Wa_mass", "gen_Wb_mass",
    "reco_matched_Wa_mass", "reco_matched_Wb_mass",

    # --- BW pairing (Voigt) ---
    "bwpair_pairing", "bwpair_gof_best", "bwpair_prob_best", "bwpair_dgof",
    "bwpair_gof0", "bwpair_gof1", "bwpair_gof2",
    "bwpair_prob0", "bwpair_prob1", "bwpair_prob2",
    "bwpair_ma0", "bwpair_ma1", "bwpair_ma2",
    "bwpair_mb0", "bwpair_mb1", "bwpair_mb2",
    "bwpair_correct",

    # --- BW pairing (pure BW, no detector smearing) ---
    "bwpair_bw_pairing", "bwpair_bw_gof_best", "bwpair_bw_prob_best", "bwpair_bw_dgof",
    "bwpair_bw_correct",

]





#print('saving these branches',all_branches)
# Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis:

    # __________________________________________________________
    # Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):


    

        # __________________________________________________________
        # Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2




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

        # change hard to production
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


        # splits Ws into on and off shell
        # instead we want to use approach where we just get daughter quarks sorted into W1 and W2 that are not ordered any way in particular

        # peak center and width
        df = df.Define(
            "Ws_on_and_off_shell",
            "FCCAnalyses::ZHfunctions::get_on_and_off_shell_WW_160ecm(HardWs_all, Particle, 80.4, 2.1)"
        )

        df = df.Define("W_on_shell_idx",  "Ws_on_and_off_shell.high_mass_idx")
        df = df.Define("W_off_shell_idx", "Ws_on_and_off_shell.low_mass_idx")
        df = df.Define("W_on_shell_mass", "Ws_on_and_off_shell.high_mass")
        df = df.Define("W_off_shell_mass", "Ws_on_and_off_shell.low_mass")


        df = df.Alias("Particle0", "Particle#0.index")


        # returns as particles
        df = df.Define("all_W_quarks_obj",
                    f"FCCAnalyses::WWFunctions::sel_quarks_fromBoson({_boson_pdg})(Particle, Particle0)")

        # hadronic cut
        df = df.Filter("all_W_quarks_obj.size() == 4", f"gen: 4 quarks from 2 {BW_BOSON} (signal)")

        nJets = 4  # WW → qqqq fully hadronic: 4 quarks → 4 jets

        collections_noleps = copy.deepcopy(collections)
        # collections_noleps["PFParticles"] = "ReconstructedParticlesNoMuNoEl"
        if saveExclJets:
            jetClusteringHelper = ExclusiveJetClusteringHelper(
                collections_noleps["PFParticles"], nJets
            )

        jetClusteringHelper_R5  = InclusiveJetClusteringHelper(
            collections_noleps["PFParticles"], 0.5, 10, "R5"
        )
        if saveExclJets:
            df = jetClusteringHelper.define(df)
            df = df.Define("d_45", "JetClusteringUtils::get_exclusive_dmerge(_jet, 4)")
            if BW_SQRTD45 > 0:
                df = df.Filter(f"d_45 < {BW_SQRTD45 * BW_SQRTD45}",
                               f"genuine 4-jet: sqrt(d_45) < {BW_SQRTD45:g} GeV")
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
        if saveExclJets: df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)
        df = jetFlavourHelper_R5.inference(weaver_preproc, weaver_model, df)



        # loop over all quarks pair decays from W and saves the events by object to "W_to{q1}{q2}"

        # start writing my python code to treemaker instead

        df = df.Alias("DaughterIndex", "Particle#1.index")


        dup_tracker = []




        # what we call gen_quarks now


        df = df.Define("All_W_quarks_idx", f"FCCAnalyses::WWFunctions::sel_quarks_fromBoson_idx({_boson_pdg})(Particle, Particle0)")
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




        if saveExclJets:
            df = df.Define(
                "jets_p4",
                "JetConstituentsUtils::compute_tlv_jets({})".format(
                    jetClusteringHelper.jets
                ),
            )
            for i in (1, 2, 3, 4):
                df = df.Define(f"jet{i}", f"jets_p4[{i-1}]")
            df = df.Define("n_reco_jets", "(int)jets_p4.size()")
            df = df.Filter("n_reco_jets == 4", "exactly 4 reco jets")   # safety (very rare to fail)

        # TRUTH MATCHING
        # BEST METHOD FOR MATCHING J TO Q, WE WANT TO SAVE THIS
        # matching jets to quarks using a chi-squared fit using dPhi and dEta
        # sigmas defined in order of {eta, phi}
        # sigmas estimated from 5,000 events when using old deltaR algorithim (no chi-squared fit), {s_eta = 0.2480, s_phi = 0.2389}

        # the indexs will be our perms from bw_pairing, matches jet to quark like
        # perm = {2, 0, 3, 1}
        # jet 0 → quark 2
        # jet 1 → quark 0
        # jet 2 → quark 3
        # jet 3 → quark 1

        df = df.Define(
            'chi2_matched_jets_to_q_etaphi',
            f'FCCAnalyses::ZHfunctions::JtoQ_ChiSquared_eta_phi(jets_p4, all_W_quarks_tlv, ROOT::VecOps::RVec<double>({{{SIGMA_ETA}, {SIGMA_PHI}}}), 0.1)'
        )

        df = df.Define('etaphi_perm',              'chi2_matched_jets_to_q_etaphi.idx')
        df = df.Define('chi2_etaphi_under_min_delR',   'chi2_matched_jets_to_q_etaphi.under_min_delR')
        df = df.Define('chi2_etaphi_delta_Rs',         'chi2_matched_jets_to_q_etaphi.delta_Rs')
        df = df.Define('chi2_etaphi_delta_etas',       'chi2_matched_jets_to_q_etaphi.delta_etas')
        df = df.Define('chi2_etaphi_delta_phis',       'chi2_matched_jets_to_q_etaphi.delta_phis')
        df = df.Define('chi2_etaphi_best_chi2',        'chi2_matched_jets_to_q_etaphi.best_chi2')
        df = df.Define('chi2_etaphi_second_best_chi2', 'chi2_matched_jets_to_q_etaphi.second_best_chi2')
        df = df.Define('chi2_etaphi_third_best_chi2',  'chi2_matched_jets_to_q_etaphi.third_best_chi2')

        # chi2 using ΔR only (sigma=1 → equivalent to minimising Σ ΔR²)
        df = df.Define(
            'chi2_matched_jets_to_q_dR',
            'FCCAnalyses::ZHfunctions::JtoQ_ChiSquared_deltaR(jets_p4, all_W_quarks_tlv, ROOT::VecOps::RVec<double>({1.0}), 0.1)'
        )
        df = df.Define('dR_perm',            'chi2_matched_jets_to_q_dR.idx')
        df = df.Define('chi2_dR_under_min_delR', 'chi2_matched_jets_to_q_dR.under_min_delR')
        df = df.Define('chi2_dR_delta_Rs',       'chi2_matched_jets_to_q_dR.delta_Rs')
        df = df.Define('chi2_dR_best_chi2',      'chi2_matched_jets_to_q_dR.best_chi2')


        # perm from previous script works approximately the same as the chi2 approach, this is almost identical
        # but takes minimum sum of dR instead of minimum sum of squares

        # NEED TO ATTEMPT MATCHING FIRST

        df = df.Define("bwpair", f"FCCAnalyses::WWFunctions::bwPairing(jet1, jet2, jet3, jet4, {SIGMA_W_ON_SHELL})")
        df = df.Define("bwpair_pairing",   "bwpair.pairing")
        df = df.Define("bwpair_gof_best",  "bwpair.gof_best")
        df = df.Define("bwpair_prob_best", "bwpair.prob_best")
        df = df.Define("bwpair_dgof",      "bwpair.dgof")
        for k in range(3):
            df = df.Define(f"bwpair_gof{k}",  f"bwpair.gof[{k}]")
            df = df.Define(f"bwpair_prob{k}", f"bwpair.prob[{k}]")
            df = df.Define(f"bwpair_ma{k}",   f"bwpair.m_a[{k}]")
            df = df.Define(f"bwpair_mb{k}",   f"bwpair.m_b[{k}]")

        
        # Pure BW pairing (no detector smearing) — saved alongside Voigt for comparison
        df = df.Define("bwpair_bw", "FCCAnalyses::WWFunctions::bwPairingBW(jet1, jet2, jet3, jet4)")
        df = df.Define("bwpair_bw_pairing",   "bwpair_bw.pairing")
        df = df.Define("bwpair_bw_gof_best",  "bwpair_bw.gof_best")
        df = df.Define("bwpair_bw_prob_best", "bwpair_bw.prob_best")
        df = df.Define("bwpair_bw_dgof",      "bwpair_bw.dgof")







        # TRUTH MATCHING

        for k in range(4):
            df = df.Define(f"gen_q{k}", f"all_W_quarks_tlv[{k}]")

        for i in (1, 2, 3, 4):
            k = i - 1
            df = df.Define(f"gen_quark{i}", f"all_W_quarks_tlv[dR_perm[{k}]]")
            # how close jet i is to its quark (used to call an event "matched")
            df = df.Define(f"jet{i}_matched_q_dR", f"(double)jet{i}.DeltaR(gen_quark{i})")
            # W-label of jet i: quarks 0,1 -> boson A (0); quarks 2,3 -> boson B (1)
            df = df.Define(f"jet{i}_wlab", f"(int)(dR_perm[{k}] >= 2)")
        # the true pairing index {0,1,2} from the 4 jet W-labels
        df = df.Define("gen_pairing_true",
            "FCCAnalyses::WWFunctions::pairing_index_from_groups("
            "jet1_wlab, jet2_wlab, jet3_wlab, jet4_wlab)")
        # did each tool pick the true pairing? (1/0)  <- efficiency numerator
        df = df.Define("bwpair_correct",
            "(int)(gen_pairing_true >= 0 && bwpair.pairing == gen_pairing_true)")
        df = df.Define("bwpair_bw_correct",
            "(int)(gen_pairing_true >= 0 && bwpair_bw.pairing == gen_pairing_true)")
        # gen-level W masses: quarks 0,1 -> boson A; quarks 2,3 -> boson B
        df = df.Define("gen_Wa_mass", "(float)(gen_q0 + gen_q1).M()")
        df = df.Define("gen_Wb_mass", "(float)(gen_q2 + gen_q3).M()")
        # reco di-jet masses of the truth-matched W pairs.
        # m_a in bwpair always contains jet1 (see order[3][4] in BWPairing.h),
        # so jet1_wlab tells us whether m_a corresponds to W_A or W_B.
        df = df.Define("reco_matched_Wa_mass",
            "gen_pairing_true >= 0 ? "
            "(jet1_wlab == 0 ? bwpair.m_a[gen_pairing_true] : bwpair.m_b[gen_pairing_true]) "
            ": -1.0f")
        df = df.Define("reco_matched_Wb_mass",
            "gen_pairing_true >= 0 ? "
            "(jet1_wlab == 0 ? bwpair.m_b[gen_pairing_true] : bwpair.m_a[gen_pairing_true]) "
            ": -1.0f")







        # exclusive at 4 jets
        
        df = df.Define(
            "jets_R5_p4",
            "JetConstituentsUtils::compute_tlv_jets({})".format(
                jetClusteringHelper_R5.jets
            ),
        )

        df = df.Define("jets_R5_p",           "JetClusteringUtils::get_p({})".format(jetClusteringHelper_R5.jets))
        df = df.Define("jets_R5_theta",       "JetClusteringUtils::get_theta({})".format(jetClusteringHelper_R5.jets))


        df = df.Define("jets_R5_pflavor", "JetTaggingUtils::get_flavour({}, Particle)".format(jetClusteringHelper_R5.jets) )
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
        
        df = df.Define("jets_R5_btagged_true", "JetTaggingUtils::sel_tag(true)(jets_R5_btag_true,{})".format(jetClusteringHelper_R5.jets))
        df = df.Define("jets_R5_ctagged_true", "JetTaggingUtils::sel_tag(true)(jets_R5_ctag_true,{})".format(jetClusteringHelper_R5.jets))
        df = df.Define("jets_R5_ltagged_true", "JetTaggingUtils::sel_tag(true)(jets_R5_ltag_true,{})".format(jetClusteringHelper_R5.jets))
        df = df.Define("jets_R5_gtagged_true", "JetTaggingUtils::sel_tag(true)(jets_R5_gtag_true,{})".format(jetClusteringHelper_R5.jets))

        df = df.Define("nbjets_R5_true", "return int(jets_R5_btagged_true.size())")
        df = df.Define("ncjets_R5_true", "return int(jets_R5_ctag_true.size())")
        df = df.Define("nljets_R5_true", "return int(jets_R5_ltag_true.size())")
        df = df.Define("ngjets_R5_true", "return int(jets_R5_gtag_true.size())")


        df = df.Define("bjets_R5_true",  "JetConstituentsUtils::compute_tlv_jets({})". format('jets_R5_btagged_true'))
        df = df.Define("ljets_R5_true",  "JetConstituentsUtils::compute_tlv_jets({})". format('jets_R5_ltagged_true'))
        df = df.Define("bjet1_R5_true_p","nbjets_R5_true > 0 ? bjets_R5_true[0].P() : -999")
        df = df.Define("ljet1_R5_true_p","nljets_R5_true > 0 ? ljets_R5_true[0].P() : -999")

        

        
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
        
            # df = df.Define("d_12", "JetClusteringUtils::get_exclusive_dmerge(_jet, 1)")
            # df = df.Define("d_23", "JetClusteringUtils::get_exclusive_dmerge(_jet, 2)")
            # df = df.Define("d_34", "JetClusteringUtils::get_exclusive_dmerge(_jet, 3)")
            # df = df.Define("d_45", "jets_p4.size()>4 ? JetClusteringUtils::get_exclusive_dmerge(_jet, 4) : -999")
            # df = df.Define("d_56", "jets_p4.size()>5 ? JetClusteringUtils::get_exclusive_dmerge(_jet, 5) : -999")


            # df = df.Define("jet1_isG", "JetFlavourUtils::get_weight(MVAVec_, 0)[0]")
            # df = df.Define("jet2_isG", "JetFlavourUtils::get_weight(MVAVec_, 0)[1]")
            # df = df.Define("jet3_isG", "JetFlavourUtils::get_weight(MVAVec_, 0)[2]")
            # df = df.Define("jet4_isG", "JetFlavourUtils::get_weight(MVAVec_, 0)[3]")
            # df = df.Define("jet5_isG", "jets_p4.size()>4 ? JetFlavourUtils::get_weight(MVAVec_, 0)[4] : -999")
            # df = df.Define("jet6_isG", "jets_p4.size()>5 ? JetFlavourUtils::get_weight(MVAVec_, 0)[5] : -999")
            
            # df = df.Define("jet1_isQ", "JetFlavourUtils::get_weight(MVAVec_, 1)[0]")
            # df = df.Define("jet2_isQ", "JetFlavourUtils::get_weight(MVAVec_, 1)[1]")
            # df = df.Define("jet3_isQ", "JetFlavourUtils::get_weight(MVAVec_, 1)[2]")
            # df = df.Define("jet4_isQ", "JetFlavourUtils::get_weight(MVAVec_, 1)[3]")
            # df = df.Define("jet5_isQ", "jets_p4.size()>4 ? JetFlavourUtils::get_weight(MVAVec_, 1)[4] : -999")
            # df = df.Define("jet6_isQ", "jets_p4.size()>5 ? JetFlavourUtils::get_weight(MVAVec_, 1)[5] : -999")
            
            # df = df.Define("jet1_isS", "JetFlavourUtils::get_weight(MVAVec_, 2)[0]")
            # df = df.Define("jet2_isS", "JetFlavourUtils::get_weight(MVAVec_, 2)[1]")
            # df = df.Define("jet3_isS", "JetFlavourUtils::get_weight(MVAVec_, 2)[2]")
            # df = df.Define("jet4_isS", "JetFlavourUtils::get_weight(MVAVec_, 2)[3]")
            # df = df.Define("jet5_isS", "jets_p4.size()>4 ? JetFlavourUtils::get_weight(MVAVec_, 2)[4] : -999")
            # df = df.Define("jet6_isS", "jets_p4.size()>5 ? JetFlavourUtils::get_weight(MVAVec_, 2)[5] : -999")
            
            # df = df.Define("jet1_isC", "JetFlavourUtils::get_weight(MVAVec_, 3)[0]")
            # df = df.Define("jet2_isC", "JetFlavourUtils::get_weight(MVAVec_, 3)[1]")
            # df = df.Define("jet3_isC", "JetFlavourUtils::get_weight(MVAVec_, 3)[2]")
            # df = df.Define("jet4_isC", "JetFlavourUtils::get_weight(MVAVec_, 3)[3]")
            # df = df.Define("jet5_isC", "jets_p4.size()>4 ? JetFlavourUtils::get_weight(MVAVec_, 3)[4] : -999")
            # df = df.Define("jet6_isC", "jets_p4.size()>5 ? JetFlavourUtils::get_weight(MVAVec_, 3)[5] : -999")
    
            # df = df.Define("jet1_isB", "JetFlavourUtils::get_weight(MVAVec_, 4)[0]")
            # df = df.Define("jet2_isB", "JetFlavourUtils::get_weight(MVAVec_, 4)[1]")
            # df = df.Define("jet3_isB", "JetFlavourUtils::get_weight(MVAVec_, 4)[2]")
            # df = df.Define("jet4_isB", "JetFlavourUtils::get_weight(MVAVec_, 4)[3]")
            # df = df.Define("jet5_isB", "jets_p4.size()>4 ? JetFlavourUtils::get_weight(MVAVec_, 4)[4] : -999")
            # df = df.Define("jet6_isB", "jets_p4.size()>5 ? JetFlavourUtils::get_weight(MVAVec_, 4)[5] : -999")


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
