import pandas as pd
import matplotlib.pyplot as plt
import mplhep as hep
import uproot3 as uproot
import seaborn as sns
import numpy as np
import json
import os
from scipy import stats
plt.style.use(hep.style.CMS)

def clopper_pearson(x, n, alpha=0.32, return_errors=True):
    """Estimate the confidence interval for binomial distributions.
    `x` is the number of successes and `n` is the number trials (x <=
    n). `alpha` is the confidence level (i.e., the true probability is
    inside the confidence interval with probability 1-alpha). The
    function returns a `(low, high)` pair of numbers indicating the
    interval on the probability.
    https://root.cern.ch/doc/master/classTEfficiency.html#ae80c3189bac22b7ad15f57a1476ef75b
    """
    b = stats.beta.ppf
    #if isinstance(x, np.ndarray):
    #    alpha = alpha*np.ones_like(x)
    ratio = x/n
    ratio = np.nan_to_num(ratio, nan=0)
    
    lo = b(alpha / 2, x, n - x + 1)
    hi = b(1 - alpha / 2, x + 1, n - x)
    if isinstance(x, np.ndarray):
        lo = np.nan_to_num(lo, nan=0)
        hi = np.nan_to_num(hi, nan=1)
        if return_errors:
            lo = ratio - lo
            hi = hi -ratio
            hi = np.where((ratio==0) & (hi>0.5), 0, hi )
        to_return = np.array([lo, hi])
    else:
        to_return = [0.0 if math.isnan(lo) else lo, 
                    1.0 if math.isnan(hi) else hi]
        if return_errors:
            to_return[0] = ratio - to_return[0]
            to_return[1] = to_return[1] -ratio
    return to_return


class npEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
    

Bins1d = dict(
    muProbe_pt = [0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,40,50],
    muProbe_eta = np.linspace(-2.4, 2.4, 24),
    DiMu_mass = [2.9,2.95, 3, 3.05, 3.1, 3.15 ,3.2, 3.25, 3.3],
    muProbe_phi = np.linspace(-np.pi, np.pi, 20),
    dR_muons = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,  0.7],
    lxySig = [0, 0.5, 1, 1.5, 2, 2.5, 3,  3.5, 4, 4.5, 5, 5.5]
)

pretty_name=dict(
    muProbe_pt = r'Offline $\mu_{pT}$',
    muProbe_eta = r'Offline $\mu_{\eta}$',
    DiMu_mass = r'$m_{\mu\mu}$',
    muProbe_phi = r'Offline $\mu_{\phi}$',
    dR_muons="$\Delta R(\mu_{tag}, \mu_{probe})$",
)


if __name__== '__main__':
    
    tagPath = 'HLT_Mu8_v'
    probePath = 'HLT_Mu0_L1DoubleMu_v'
    input_file = '/eos/cms/store/group/phys_bphys/trigger/Run2024/TagTuples/Muon_Run2024C.root'
    run  = '2024C'
    Name = 'L1Efficiency_v0'

    tagQuery = "2.9<DiMu_mass<3.3"\
             +"& DiMu_Prob>0.005 "\
             +"& abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 "\
             +"& muTag_pt>8 "\
             +"& muTag_L1_match==1  & muProbe_L1_match==1 "\
             +"& muTag_charge+muProbe_charge==0"\
    


    tagPath = 'HLT_Mu8_v'
    probePath = 'HLT_Mu0_L1DoubleMu_v'
    input_file = '/eos/cms/store/group/phys_bphys/trigger/Run2024/TagTuples/Muon_Run2024C_v1.root'
    run  = '2024C'
    Name = 'L1Efficiency_STEAM_May'

    tagQuery = "2.9<DiMu_mass<3.3"\
             +"& DiMu_Prob>0.005 "\
             +"& abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 "\
             +"& muTag_pt>8 "\
             +"& muTag_L1_match==1  & muProbe_L1_match==1 "\
             +"& muTag_charge+muProbe_charge==0"\



    # tagPath = 'HLT_Mu4_L1DoubleMu_v'
    # probePath = 'HLT_DoubleMu4_3_LowMass_v'
    # input_file = '/eos/cms/store/group/phys_bphys/trigger/Run2024/TagTuples/Muon_Run2024C.root'
    # run  = '2024C'
    # Name = 'HLTEfficiency_v0'

    # tagQuery = "2.9<DiMu_mass<3.3"\
    #          +"& DiMu_Prob>0.005 "\
    #          +"& abs(muTag_eta)<2.4 & abs(muProbe_eta)<2.4 "\
    #          +"& muTag_charge+muProbe_charge==0"\
    #          #+"& muTag_L1_match==1  & muProbe_L1_match==1 "\             




    arrays = ['DiMu_mass', 'DiMu_Prob', 'event', '*HLT_*' , 'mu*match', '*lxy*', '*charge*', 'L1*', '*dR*'] #+ ['DiMu_Prob']
    arrays+= 'muProbe_pt,muProbe_eta,muProbe_phi,muTag_pt,muTag_eta,muTag_phi'.split(',')
    arrays+= 'L3_muProbe_pt,L3_muProbe_eta,L3_muProbe_phi,L3_muTag_pt,L3_muTag_eta,L3_muTag_phi'.split(',')



    outputdir = os.path.join('Run'+run, Name)
    os.makedirs(outputdir, exist_ok=True)
    file = uproot.open(input_file)
    data = file[tagPath].arrays( outputtype=pd.DataFrame)
    
    
    if 'HLT_Mu8' in tagPath and 'HLT_Mu0_L1' in probePath:
        efficiencyText = 'L1 Efficiency'
    elif 'HLT_Mu4_L1' in tagPath and 'HLT_Double' in probePath:
        efficiencyText = 'HLT Efficiency'
    else:
        efficiencyText = 'Efficiency'

    for var1 in ['muProbe_pt', 'muProbe_eta', 'muProbe_phi', 'dR_muons','lxySig']:
    
    
        h_all      = np.histogram(data.query(tagQuery)[var1], bins=Bins1d[var1])
        h_passprob = np.histogram(data.query(tagQuery+f' and muProbe_{probePath}==1')[var1], 
                                 bins=Bins1d[var1])        
        ratio = h_passprob[0]/h_all[0]
        err  = clopper_pearson(h_passprob[0], h_all[0])

        
        efficiency_output = dict(
            Bins    = h_all[1],
            Passing = h_passprob[0],
            All     = h_all[0],
            Ratio   = ratio,
            Error   = err,
            Probe   = probePath,
            Tag     = tagPath,
            Query   = tagQuery,
            Var     = var1,
            input   = input_file
        )
        
        with open(f'{outputdir}/Efficiency1D_{var1}.json', 'w+') as jj:
            json.dump(efficiency_output, jj, indent=4, cls=npEncoder)
        
        
        bin_mean = (h_all[1][1:] + h_all[1][:-1])/2
        bin_size = (h_all[1][1:] - h_all[1][:-1])/2

        
        fig,ax = plt.subplots(figsize=[15,10])
        ax.errorbar(bin_mean, ratio, err, xerr=bin_size, ls='none', marker='o', capsize=2 )    
        ax.set_ylabel(efficiencyText)
        ax.set_xlabel(pretty_name.get(var1, var1))
        #ax.legend(frameon=True, title='Run 2024D')
        hep.cms.label(data=True, label=run, com=13.6)
        ax.set_ylim(0, 1.2)
        ax.grid(True)
        plt.savefig(f'{outputdir}/Tag{tagPath}_Probe{probePath}_{var1}.pdf', bbox_inches='tight')
        plt.close()