#!/usr/bin/env python3
"""Print chi2_etaphi delta_R, delta_eta, delta_phi values for each event."""
import uproot
import numpy as np

WW_FN = ("/eos/user/m/mlevere/ttThreshold-analysis/outputs/treemaker/WbWb/"
         "p8_ee_WW_ecm240_new_matching/p8_ee_WW_ecm240.root")

t = uproot.open(WW_FN)["events"]
C = t.arrays(["chi2_etaphi_delta_Rs", "chi2_etaphi_delta_etas", "chi2_etaphi_delta_phis"],
             library="np")

dRs   = C["chi2_etaphi_delta_Rs"]
detas = C["chi2_etaphi_delta_etas"]
dphis = C["chi2_etaphi_delta_phis"]

print(f"{'event':>6}  {'delta_Rs':>40}  {'delta_etas':>40}  {'delta_phis':>40}")
print("-" * 135)
for i in range(min(20, len(dRs))):
    print(f"{i:>6}  {str(dRs[i]):>40}  {str(detas[i]):>40}  {str(dphis[i]):>40}")
