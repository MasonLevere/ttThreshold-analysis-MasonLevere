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

# # One on shell, one off shell WW events
#     "p8_ee_WW_ecm160": {
#         "fraction": 1,
#     },

# # both very on shell
#     "p8_ee_WW_ecm240": {
#         "fraction": 1,
#     }

    "mgp8_ee_zh_ecm240": {
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




available_ecm = ['125', '160', '163', '240', '340','345', '350', '355','365']

hadronic  = True
#semihad  = False
#lep      = False
ecm = int(os.environ.get("ecm", list(all_processes.keys())[-1][-3:]))
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

processList={key: value for key, value in all_processes.items() if str(ecm) in available_ecm and str(ecm) in key }

# override sample via env var (e.g. BW_SAMPLE=p8_ee_ZZ_ecm160 for ZZ control sample)
_BW_SAMPLE = os.environ.get("BW_SAMPLE", "").strip()
_sample_tag = _BW_SAMPLE if _BW_SAMPLE else f"mgp8_ee_zh_ecm{ecm}"

inputDir    = "/eos/user/m/mlevere/ttThreshold-analysis/localSamples/"

# compute fraction to limit events: count available files, assume 100k events each
_n_events = int(os.environ.get("N_EVENTS", 0))
import glob as _glob
_n_files = len(_glob.glob(os.path.join(inputDir, _sample_tag, "*.root")))
if _n_events > 0 and _n_files > 0:
    _total_events = _n_files * 100000
    _fraction = min(1.0, _n_events / _total_events)
else:
    _fraction = 1.0

processList = {_sample_tag: {"fraction": _fraction}}

print(f"processList: {processList}  (fraction={_fraction:.4f}, {_n_files} files, ~{int(_fraction*_n_files*100000)} events)")

# Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics (mandatory)
#prodTag     = "FCCee/winter2023/IDEA/"

_sigma_str  = os.environ.get("BW_SIGMA", "").strip()
_sigma_tag  = f"_sig{float(_sigma_str):.4f}" if _sigma_str else ""
outputDir   = "/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/{}{}_new_matching".format(_sample_tag, _sigma_tag)


# additional/costom C++ functions, defined in header files (optional)
#                   namespace is ZH         namespace is WWFunctions
includePaths = ["WWfunctions/functions.h", "WWfunctions/WWFunctions.h", "WWfunctions/BWPairing.h"]

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

BW_BOSON   = os.environ.get("BW_BOSON", "W").strip().upper()
BW_SQRTD45 = float(os.environ.get("BW_SQRTD45_MAX", "7.0"))


_boson_pdg   = {"W": 24,   "Z": 23}[BW_BOSON]
_boson_mass  = {"W": 80.4, "Z": 91.2}[BW_BOSON]
_boson_width = {"W": 2.1,  "Z": 2.495}[BW_BOSON]



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

# refers to detector resolution for gauss in voigt distribution; override with BW_SIGMA=<val>
sigma = float(os.environ.get("BW_SIGMA", "3.6110261681321907"))

_bw2d_table_dir = "/afs/cern.ch/user/m/mlevere/private/FCCTutorial/ttThreshold-analysis/bw2d_tables"
_bw2d_table_path = f"{_bw2d_table_dir}/bw2d_mWW{float(ecm):.1f}_mw80.419_gw2.049_sig{sigma:.4f}_dm0.05.bin"
os.environ.setdefault("BW2D_TABLE_PATH", _bw2d_table_path)
print(f"[treemaker] BW2D_TABLE_PATH = {os.environ['BW2D_TABLE_PATH']}")

chi2_cut = 100



all_branches = [

    # --- jets ---
    "jets_p4",
    "n_reco_jets",
    "d_45",

    # --- chi2 jet-quark matching ---
    "chi2_dR_best_chi2", "chi2_dR_delta_Rs",

    # --- truth matching ---
    "jet1_matched_q_dR", "jet2_matched_q_dR",

    # --- H jet kinematics ---
    "H_dR_between_jets", "H_dTheta_between_jets",
    "H_P_from_jets", "H_Pt_from_jets", "H_E_from_jets",

]





#print('saving these branches',all_branches)
# Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis:

    # __________________________________________________________
    # Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        import ROOT
        ROOT.FCCAnalyses.WWFunctions.init_bw2d_table()

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

        df = df.Alias("DaughterIndex", "Particle#1.index")
        df = df.Alias("Particle0", "Particle#0.index")

        # # returns as particles
        # # get quarks from Z decays
        # df = df.Define("all_W_quarks_obj",
        #             "FCCAnalyses::ZHfunctions::sel_daughter_fromParent(23, 2, {1, 2, 3, 4, 5, 6})(Particle, Particle0)")

        # get nunu decays from Z decay
        # gets antiparticles too
        df = df.Define("Z_to_nunu",
                    "FCCAnalyses::ZHfunctions::sel_daughter_fromParent(23, 2, {12, 14, 16})(Particle, Particle0)")

        # hadronic cut (keeps nunu and ll)
        df = df.Filter("Z_to_nunu.size() == 2", f"gen: keep nunu Z decays (signal)")

        # add a cut for H to qq
        # doesnt include gg emission
        df = df.Define("quarks_from_higgs",
            f"FCCAnalyses::ZHfunctions::sel_quarks_fromHiggs({25})(Particle, Particle0)")
        df = df.Filter("quarks_from_higgs.size() == 2", "keep only hadronic higgs decay")


        nJets = 2  

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
            df = df.Define("d_45", f"JetClusteringUtils::get_exclusive_dmerge(_jet, {nJets})")
            # if BW_SQRTD45 > 0:
                # df = df.Filter(f"d_45 < {BW_SQRTD45 * BW_SQRTD45}",
                #                f"genuine 4-jet: sqrt(d_45) < {BW_SQRTD45:g} GeV")
                # instead of a filter, just create a collection
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



        dup_tracker = []



        df = df.Define("q_tlv", "FCCAnalyses::MCParticle::get_tlv(quarks_from_higgs)")



        if saveExclJets:
            df = df.Define(
                "jets_p4",
                "JetConstituentsUtils::compute_tlv_jets({})".format(
                    jetClusteringHelper.jets
                ),
            )
            for i in (1, 2):
                df = df.Define(f"jet{i}", f"jets_p4[{i-1}]")
            df = df.Define("n_reco_jets", "(int)jets_p4.size()")
            df = df.Filter("n_reco_jets == 2", "exactly 2 reco jets")   # safety (very rare to fail)

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

        # removed for now since etaphi require knowledge of uncertainty in both

        # df = df.Define(
        #     'chi2_matched_jets_to_q_etaphi',
        #     f'FCCAnalyses::ZHfunctions::JtoQ_ChiSquared_eta_phi(jets_p4, all_W_quarks_tlv, ROOT::VecOps::RVec<double>({{{SIGMA_ETA}, {SIGMA_PHI}}}), 0.1)'
        # )

        # df = df.Define('etaphi_perm',              'chi2_matched_jets_to_q_etaphi.idx')
        # df = df.Define('chi2_etaphi_under_min_delR',   'chi2_matched_jets_to_q_etaphi.under_min_delR')
        # df = df.Define('chi2_etaphi_delta_Rs',         'chi2_matched_jets_to_q_etaphi.delta_Rs')
        # df = df.Define('chi2_etaphi_delta_etas',       'chi2_matched_jets_to_q_etaphi.delta_etas')
        # df = df.Define('chi2_etaphi_delta_phis',       'chi2_matched_jets_to_q_etaphi.delta_phis')
        # df = df.Define('chi2_etaphi_best_chi2',        'chi2_matched_jets_to_q_etaphi.best_chi2')
        # df = df.Define('chi2_etaphi_second_best_chi2', 'chi2_matched_jets_to_q_etaphi.second_best_chi2')
        # df = df.Define('chi2_etaphi_third_best_chi2',  'chi2_matched_jets_to_q_etaphi.third_best_chi2')

        # for i in range(1, 5):
        #     df = df.Define(f"simple_jet_{i}_deltaR", f"chi2_etaphi_delta_Rs[{i-1}]")
        #     df = df.Define(f"simple_jet_{i}_delta_eta", f"chi2_etaphi_delta_etas[{i-1}]")
        #     df = df.Define(f"simple_jet_{i}_delta_phi", f"chi2_etaphi_delta_phis[{i-1}]")




        # chi2 using ΔR only (sigma=1 → equivalent to minimising Σ ΔR²)
        df = df.Define(
            'chi2_matched_jets_to_q_dR',
            'FCCAnalyses::ZHfunctions::JtoQ_ChiSquared_deltaR(jets_p4, q_tlv, ROOT::VecOps::RVec<double>({1.0}), 0.1)'
        )

        df = df.Define('dR_perm',            'chi2_matched_jets_to_q_dR.idx')
        df = df.Define('chi2_dR_under_min_delR', 'chi2_matched_jets_to_q_dR.under_min_delR')
        df = df.Define('chi2_dR_delta_Rs',       'chi2_matched_jets_to_q_dR.delta_Rs')
        df = df.Define('chi2_dR_best_chi2',      'chi2_matched_jets_to_q_dR.best_chi2')


        # perm from previous script works approximately the same as the chi2 approach, this is almost identical
        # but takes minimum sum of dR instead of minimum sum of squares

        # matching not needed for now, just looking at kinematics
        # NEED TO ATTEMPT MATCHING FIRST

        # df = df.Define("bwpair", f"FCCAnalyses::WWFunctions::bwPairing(jet1, jet2, jet3, jet4, {SIGMA_W_ON_SHELL})")
        # df = df.Define("bwpair_pairing",   "bwpair.pairing")
        # df = df.Define("bwpair_gof_best",  "bwpair.gof_best")
        # df = df.Define("bwpair_prob_best", "bwpair.prob_best")
        # df = df.Define("bwpair_dgof",      "bwpair.dgof")
        # for k in range(3):
        #     df = df.Define(f"bwpair_gof{k}",  f"bwpair.gof[{k}]")
        #     df = df.Define(f"bwpair_prob{k}", f"bwpair.prob[{k}]")
        #     df = df.Define(f"bwpair_ma{k}",   f"bwpair.m_a[{k}]")
        #     df = df.Define(f"bwpair_mb{k}",   f"bwpair.m_b[{k}]")

        
        # # Pure BW pairing (no detector smearing) — saved alongside Voigt for comparison
        # df = df.Define("bwpair_bw", "FCCAnalyses::WWFunctions::bwPairingBW(jet1, jet2, jet3, jet4)")
        # df = df.Define("bwpair_bw_pairing",   "bwpair_bw.pairing")
        # df = df.Define("bwpair_bw_gof_best",  "bwpair_bw.gof_best")
        # df = df.Define("bwpair_bw_prob_best", "bwpair_bw.prob_best")
        # df = df.Define("bwpair_bw_dgof",      "bwpair_bw.dgof")


        # # now we want to try to use 2D BW

        # # 2D BW pairing (phase-space weighted double BW, no detector smearing)
        # # doubleBWPairing takes m_WW (GeV) and squares it internally -> pass ecm, NOT ecm*ecm
        # df = df.Define("double_bwpair", f"FCCAnalyses::WWFunctions::doubleBWPairing(jet1, jet2, jet3, jet4, {ecm})")
        # df = df.Define("double_bwpair_pairing",   "double_bwpair.pairing")
        # df = df.Define("double_bwpair_gof_best",  "double_bwpair.gof_best")
        # df = df.Define("double_bwpair_prob_best", "double_bwpair.prob_best")
        # df = df.Define("double_bwpair_dgof",      "double_bwpair.dgof")

        # # 2D BW pairing with detector smearing (precomputed mmap table, init'd above)
        # df = df.Define("double_smeared_bwpair", "FCCAnalyses::WWFunctions::doubleBWPairingSmeared(jet1, jet2, jet3, jet4)")
        # df = df.Define("double_smeared_bwpair_pairing",   "double_smeared_bwpair.pairing")
        # df = df.Define("double_smeared_bwpair_gof_best",  "double_smeared_bwpair.gof_best")
        # df = df.Define("double_smeared_bwpair_prob_best", "double_smeared_bwpair.prob_best")
        # df = df.Define("double_smeared_bwpair_dgof",      "double_smeared_bwpair.dgof")
        # for k in range(3):
        #     df = df.Define(f"double_smeared_bwpair_gof{k}",  f"double_smeared_bwpair.gof[{k}]")
        #     df = df.Define(f"double_smeared_bwpair_prob{k}", f"double_smeared_bwpair.prob[{k}]")
        #     df = df.Define(f"double_smeared_bwpair_ma{k}",   f"double_smeared_bwpair.m_a[{k}]")
        #     df = df.Define(f"double_smeared_bwpair_mb{k}",   f"double_smeared_bwpair.m_b[{k}]")


        # TRUTH MATCHING

        for k in range(2):
            df = df.Define(f"gen_q{k}", f"q_tlv[{k}]")

        for i in (1, 2):
            k = i - 1
            df = df.Define(f"gen_quark{i}", f"q_tlv[dR_perm[{k}]]")
            # how close jet i is to its quark (used to call an event "matched")
            df = df.Define(f"jet{i}_matched_q_dR", f"(double)jet{i}.DeltaR(gen_quark{i})")
            # W-label of jet i: quarks 0,1 -> boson A (0); quarks 2,3 -> boson B (1)
            df = df.Define(f"jet{i}_wlab", f"(int)(dR_perm[{k}] >= 2)")

        # sometimes want to look at tightly cut events
        # df = df.Filter("jet1_matched_q_dR < 0.1 && jet2_matched_q_dR < 0.1 && jet3_matched_q_dR < 0.1 && jet4_matched_q_dR < 0.1", "all jets matched dR < 0.1")

        df = df.Define("H_jet_kin",
            "FCCAnalyses::ZHfunctions::compute_H_jet_kinematics("
            "jet1, jet2)")
        df = df.Define("H_dR_between_jets",     "H_jet_kin.H_dR")
        df = df.Define("H_dTheta_between_jets", "H_jet_kin.H_dTheta")
        df = df.Define("H_P_from_jets",         "H_jet_kin.H_P")
        df = df.Define("H_Pt_from_jets",        "H_jet_kin.H_Pt")
        df = df.Define("H_E_from_jets",         "H_jet_kin.H_E")


        # more truth matching we dont use for now

        # # the true pairing index {0,1,2} from the 4 jet W-labes
        # df = df.Define("gen_pairing_true",
        #     "FCCAnalyses::WWFunctions::pairing_index_from_groups("
        #     "jet1_wlab, jet2_wlab, jet3_wlab, jet4_wlab)")
        # # did each tool pick the true pairing? (1/0)  <- efficiency numerator
        # df = df.Define("bwpair_correct",
        #     "(int)(gen_pairing_true >= 0 && bwpair.pairing == gen_pairing_true)")
        # df = df.Define("bwpair_bw_correct",
        #     "(int)(gen_pairing_true >= 0 && bwpair_bw.pairing == gen_pairing_true)")
        # df = df.Define("double_bwpair_correct",
        #     "(int)(gen_pairing_true >= 0 && double_bwpair.pairing == gen_pairing_true)")
        # df = df.Define("double_smeared_bwpair_correct",
        #     "(int)(gen_pairing_true >= 0 && double_smeared_bwpair.pairing == gen_pairing_true)")
        # # gen-level W masses: quarks 0,1 -> boson A; quarks 2,3 -> boson B
        # df = df.Define("gen_Wa_mass", "(float)(gen_q0 + gen_q1).M()")
        # df = df.Define("gen_Wb_mass", "(float)(gen_q2 + gen_q3).M()")
        # # reco di-jet masses of the truth-matched W pairs.
        # # m_a in bwpair always contains jet1 (see order[3][4] in BWPairing.h),
        # # so jet1_wlab tells us whether m_a corresponds to W_A or W_B.
        # df = df.Define("reco_matched_Wa_mass",
        #     "gen_pairing_true >= 0 ? "
        #     "(jet1_wlab == 0 ? bwpair.m_a[gen_pairing_true] : bwpair.m_b[gen_pairing_true]) "
        #     ": -1.0f")
        # df = df.Define("reco_matched_Wb_mass",
        #     "gen_pairing_true >= 0 ? "
        #     "(jet1_wlab == 0 ? bwpair.m_b[gen_pairing_true] : bwpair.m_a[gen_pairing_true]) "
        #     ": -1.0f")







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
