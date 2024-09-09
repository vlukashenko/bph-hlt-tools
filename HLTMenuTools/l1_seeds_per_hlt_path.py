import ROOT
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

# Directories containing ROOT files
#directories = ["/eos/user/v/valukash/TriggerContact/bph-hlt-tools/HLTMenuTools/"]
directories = ["/eos/cms/store/group/phys_bphys/trigger/Run2024/Muon0/BPHTriggerTuples_Muon0-Run2024C-PromptReco-v1HLTStudy/240508_050810/0000/"]#
#"/eos/cms/store/group/phys_bphys/trigger/Run2024/Muon0/BPHTriggerTuples_Muon0-Run2024C-PromptReco-v1HLTStudy/240508_050810/0001/", ]
# List of L1 seeds and HLT paths
l1_seeds = [

    "L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p2",
"L1_DoubleMu4er2p0_SQ_OS_dR_Max1p6", 
"L1_DoubleMu5_SQ_OS_dR_Max1p6", 
"L1_DoubleMu3er2p0_SQ_OS_dR_Max1p6", 
"L1_DoubleMu0er1p5_SQ_OS_dEta_Max1p2", 
"L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", 
"L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6", 
"L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5", 
"L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", 
"L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", 
"L1_DoubleMu4p5_SQ_OS_dR_Max1p2", 
"L1_DoubleMu4_SQ_OS_dR_Max1p2"
]

hlt_paths = [
    "HLT_DoubleMu4_3_LowMass_v",
    "HLT_DoubleMu4_3_Bs_v",
    "HLT_DoubleMu4_LowMass_Displaced_v",
    #"HLT_DoubleMu4_3_Photon4_BsToMMG_v",
    #"HLT_DoubleMu3_TkMu_DsTau3Mu_v"
    
]

# Dictionary to store counts for each HLT path under different L1 seeds
results = {hlt: {l1: 0 for l1 in l1_seeds} for hlt in hlt_paths}

# Function to process a single ROOT file using RDataFrame
def process_root_file(file_path):
    df = ROOT.RDataFrame("rootuple/ntuple", file_path)
    #print(df.GetColumnNames())
    for hlt in hlt_paths:
        # Apply a filter for each HLT path
        df_hlt = df.Filter(f"{hlt} == 1")
        print(f"Total {hlt}: ", df_hlt.Count().GetValue())   
        for l1 in l1_seeds:
            # Count events that pass both the HLT path and the L1 seed
            count = df_hlt.Filter(f"{l1} == 1").Count().GetValue()
            print(l1, " : ", count)
            results[hlt][l1] += count

def percentage_uncertainty(p, n, N):
    if N == 0:
        return 0
    return p * np.sqrt(1 / n + 1/ N)

count = 100
count_idx = 0
# Iterate through each directory and process each ROOT file
for directory in directories:
    for filename in os.listdir(directory):
        count_idx = count_idx + 1
        print(count_idx)
        if count_idx > count: break
        if filename.endswith(".root"):
            file_path = os.path.join(directory, filename)
            print(f"Processing file: {file_path}")
            process_root_file(file_path)

# Plot pie charts for each HLT path
for hlt, l1_counts in results.items():
    labels = []
    sizes = []
    uncertainties = []
    for l1, count in l1_counts.items():
        if count > 0:
            percentage = 100 * count / sum(l1_counts.values())
            uncertainty = percentage_uncertainty(percentage, count, sum(l1_counts.values()))
            labels.append(f"{l1} ({percentage:.1f}%)")
            sizes.append(percentage)
            uncertainties.append(uncertainty)

    fig, ax = plt.subplots(figsize=(10, 10))
    wedges, texts, autotexts = ax.pie(
        sizes, labels=None, autopct='', startangle=140, colors=plt.cm.tab20.colors)

    # Customizing the pie chart
    bbox_props = dict(boxstyle="round,pad=0.3", edgecolor="w", facecolor="none")
    kw = dict(xycoords='data', textcoords='data', arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")
    
    for i, (wedge, label, size, uncertainty) in enumerate(zip(wedges, labels, sizes, uncertainties)):
        angle = (wedge.theta2 - wedge.theta1) / 2. + wedge.theta1 + i*0.5
        x = np.cos(np.deg2rad(angle))
        y = np.sin(np.deg2rad(angle))
        
        if size < 2:  # Adjust this threshold for when to plot outside
            # Annotate outside with arrow
            connectionstyle = f"angle,angleA=0,angleB={angle}"
            kw["arrowprops"].update({"connectionstyle": connectionstyle, "color": wedge.get_facecolor()})
            ax.annotate(f"{size:.1f}%", xy=(x, y), xytext=(1.4*x, 1.4*y), ha="center", **kw)
            #ax.annotate(f"{size:.1f} ± {uncertainty:.1f}%", xy=(x, y), xytext=(1.4*x, 1.4*y), ha="center", **kw)
        else:
            # Annotate inside the wedge
            ax.text(1.1*x, 1.1*y, f"{size:.1f}%", ha="center", va="center")
            #ax.text(1.1*x, 1.1*y, f"{size:.1f} ± {uncertainty:.1f}%", ha="center", va="center")
    
    # Add legend with labels
    ax.legend(wedges, labels, title="L1 Seeds",loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    plt.title(f"Percentage of events passing L1 seeds for {hlt}", ha="center")

    # Adjust the layout so everything fits
    plt.tight_layout()
    plt.savefig(hlt+"_l1_seeds_chart.png")

df_results = pd.DataFrame(columns=l1_seeds, index=hlt_paths)

for hlt, l1_counts in results.items():
    total_events = sum(l1_counts.values())
    for l1, count in l1_counts.items():
        if count > 0:
            percentage = 100 * count / total_events
            uncertainty = percentage_uncertainty(percentage, count, total_events)
            #df_results.loc[hlt, l1] = f"{percentage:.1f} ± {uncertainty:.1f}"
            df_results.loc[hlt, l1] = f"{percentage:.1f}"
 
        else:
            df_results.loc[hlt, l1] = "-"

# Save the DataFrame as a LaTeX table
df_results.to_latex("l1_seeds_per_hlt_path.tex", caption="Percentage of events passing L1 seeds for each HLT path.", label="tab:l1_seeds_per_hlt_path")

print("LaTeX table saved as l1_seeds_per_hlt_path.tex")
