import ROOT
#import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt
import os
import glob
import time
np_pi = 3.141592653589793
ROOT.EnableImplicitMT() # Tell ROOT you want to go parallel

SingleMuon = ["HLT_Mu0_Barrel_v",
            "HLT_Mu0_Barrel_L1HP10_v",
            "HLT_Mu0_Barrel_L1HP11_v",
            "HLT_Mu9_Barrel_L1HP10_IP6_v",
            "HLT_Mu10_Barrel_L1HP11_IP6_v",]

DoubleMuon = ["HLT_DoubleMu4_3_LowMass_v",
            "HLT_DoubleMu4_3_LowMass_SS_v",
            "HLT_DoubleMu4_LowMass_Displaced_v",
            "HLT_DoubleMu4_MuMuTrk_Displaced_v",
            "HLT_DoubleMu2_Jpsi_LowPt",
            "HLT_Dimuon10_Upsilon_y1p4_v",]

HLTPaths =  SingleMuon + DoubleMuon




Histos = dict(

    dimuonSpectrum = {"var_name" : "DimuonMass", 
                     "bins"      : [20000, 0, 20], 
                     "var"       : "DiMu_mass"},

    dimuonpT       = {"var_name" : "DimuonpT",   
                     "bins"      : [10000, 0, 100], 
                     "var"       : "DiMu_pt"},                     
)








def include_new_histo(name, var_name, bins, var, rdataframe):
    #Define the model
    histo_model = ROOT.RDF.TH1DModel(name, var_name, *bins)
    #Histos
    histo = rdataframe.Histo1D(histo_model, var)
    #Write
    histo.Write()

if __name__ == '__main__':


    # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    import argparse

    my_parser = argparse.ArgumentParser(prog='Histogram producer',
                                        description='This scripts produces histograms for events that were collected by a given trigger on a given variable',
                                       )
    my_parser.add_argument('--input',
                       type = str,
                       default='Rootuple_MuMu_2023C-MiniAOD-DoubleMu_*.root',
                       help='path of the input files')
    my_parser.add_argument('--outname',
                       type = str,
                       default='tuple',
                       help='path of the input files')
    my_parser.add_argument('--makePlots',                       
                       action='store_true',
                       help='Use this flag to produce default plots') 


    args  = my_parser.parse_args()
    print(args.input, args.outname)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    
    
    
    ### READ FILES ####    
    datafiles_ = glob.glob(args.input)
    print('TOTAL FILES  : :  ', len(datafiles_))
    start = time.time()
    file__ = ROOT.RDataFrame('rootuple/ntuple',
                       datafiles_)


    ### DEFINE NEW VARIABLES #####
    file_ = file__.Define("DiMu_ptLead", "DiMu_mu1_pt > DiMu_mu2_pt ? DiMu_mu1_pt : DiMu_mu2_pt")  
    file_ = file_.Define("DiMu_ptTrail", "DiMu_mu1_pt > DiMu_mu2_pt ? DiMu_mu2_pt : DiMu_mu1_pt")  
    file_ = file_.Define("Dimuon_charge", "mu1_charge+mu2_charge")

   
    #### FILTER EVENTS

    ## Same and opposite sign matching
    fileSS_ = file_.Filter('Dimuon_charge!=0') 
    file_ = file_.Filter('Dimuon_charge==0') 
    
    ## Select muons matching L3 objects 
    file_matched = file_.Filter('(mu2_L3_match==1) && (mu1_L3_match==1)')
    ## Select muons with SoftMuon flag 
    file_matched_soft = file_matched.Filter('(mu1soft==1) && (mu2soft==1)')

    # Trigger matching
    HLTPathsFilters = dict()
    for path in HLTPaths:
        HLTPathsFilters[path] = file_matched.Filter(f'{path}==1')





    #########################################
    # Create Histograms
    outputfile_ = ROOT.TFile(f"{args.outname.replace('.root', '')}.root", "RECREATE")
    
    for directory, options in Histos.items():
        outputfile_.cd()
        outputfile_.mkdir(directory)  
        root_directory = outputfile_.GetDirectory(directory)
        root_directory.cd()
        for path in HLTPathsFilters:
            include_new_histo(path, **options, rdataframe=HLTPathsFilters[path])        

    outputfile_.Close()
    print('Close File : ', round((time.time()-start)/60,2), ' min')





    ####################################
    # Plot dimuons
    if args.makePlots:

        import commons
        import uproot3 as uproot
        import matplotlib.pyplot as plt
        import matplotlib.colors as mplcolors 
        import numpy as np

        file = uproot.open(f"{args.outname.replace('.root', '')}.root")

        Histos_mumu = dict()
        for key in file['dimuonSpectrum'].keys():
            Histos_mumu[key.decode()[:-2]] = file['dimuonSpectrum'][key]
            
        # Histos_pT = dict()
        # for key in file['dimuonpT'].keys():
        #     Histos_pT[key.decode()[:-2]] = file['dimuonpT'][key]  


        ## Dimuon Spectrum
        max_value = 8.5
        _bins     = Histos_mumu['HLT_DoubleMu4_3_LowMass_v'].alledges
        mask      = _bins<=8.5
        mask2     = _bins<11.5
        mask_ups  = np.bitwise_not(mask) & (mask2)
        bins    = _bins[mask][1:]
        bins_ups=_bins[mask_ups]

        bins_sum = 50

        counts_low = dict()
        counts_high= dict()
        for key in Histos_mumu.keys():
            c, bins = commons.reduce_histogram(Histos_mumu[key].allvalues[1:-1],  
                                                    _bins[1:-1], bins_sum)
            mask      = bins<=8.5
            mask2     = bins<11.5
            mask_ups  = np.bitwise_not(mask) & (mask2)
            
            bins_dim = bins[mask][1:]
            bins_ups = bins[mask_ups]
            
            counts_low[key] = c[mask[1:]][1:]
            counts_high[key] = c[mask_ups[1:]][1:]





        frac=12*bins_sum*0.005
        fig  = plt.figure(figsize=[20,7]) 
        ax, axu = commons.create_axes(fig, space_between=0.5)

        color1 = mplcolors.to_rgba('tab:blue')
        color1 = list(color1)
        color1[0]+=0.15
        color1[1]+=0.15
        color1[2]+=0.15
        color1[-1] = 0.9
        #print('color1 : ', color1)

        color2 = mplcolors.to_rgba('navy')
        color2 = list(color2)
        color2[-2] = 1
        color2[0]+=0.2
        color2[1]+=0.2
        #print('color2 : ', color2)

        color4='purple'
        commons.hist_from_heights(counts_low['HLT_DoubleMu4_3_LowMass_v'], bins_dim, axis=ax, 
                                    color=color1, label='Inclusive low mass dimuon trigger',
                                    rasterized=True)
        commons.hist_from_heights(counts_low['HLT_DoubleMu4_LowMass_Displaced_v'], bins_dim, axis=ax, 
                                    color=color2, label='Displaced low mass dimuon trigger',
                                    rasterized=True)
        commons.hist_from_heights(counts_low['HLT_DoubleMu2_Jpsi_LowPt'], bins_dim, axis=ax, 
                                    color=color4, label=R'J/$\Psi$ low pT', histtype='step', linewidth=3,
                                    rasterized=True)
        commons.hist_from_heights(counts_low['HLT_DoubleMu4_3_LowMass_SS_v'], bins_dim, axis=ax, 
                                    color='lime', label=R'Same Sign', histtype='step', linewidth=3,
                                    rasterized=True)


        ax.set_yscale('log')
        ax.set_xlim(0.2, 8.5)
        #ax.set_ylim(30*bins_sum, 6e7*bins_sum)
        ax.legend(frameon=True, fontsize=17)
        ax.set_ylabel('Candidates / '+ str(int(1000*np.mean(np.diff(bins))))+' MeV')

        #ax.text(0.5,  7e3*frac, r'$\eta$')
        #ax.text(0.7,  2e4*frac, r'$\omega$')
        #ax.text(0.98, 2e4*frac, r'$\phi$')
        #ax.text(2.9,  1e5*frac, r'$J/\psi$')
        #ax.text(3.6,  7e3*frac, r"$\psi '$")

        commons.hep.cms.label(data=True, loc=1, label='Preliminary',                         
                                com=None, ax=ax, 
                                year=None, lumi=None)

        commons.hep.cms.label(data=True, exp='',
                                    year=None, ax=axu,
                                    #@lumi=3.2, 
                                    com=13.6)
        color3 = mplcolors.to_rgba("#df979e")
        color3 = list(color3)
        color3[-1] = 0.8
        #print('color3 : ', color3)

        #axu.text(9.9, 3e2*frac, r"$\Upsilon$")
        commons.hist_from_heights(counts_high['HLT_Dimuon10_Upsilon_y1p4_v'], bins_ups, axis=axu, 
                                    color=color3, label=r'$\Upsilon(nS)$ trigger',
                                    rasterized=True)
        axu.set_xlabel('$m_{\mu \mu}$ [GeV]', fontsize=28)
        axu.set_yscale('log')
        axu.yaxis.tick_left()
        axu.yaxis.tick_right()
        axu.set_xlim(8.6, 11.4)

        #axu.set_ylim(1e2*bins_sum, 2e4*bins_sum)
        axu.legend(frameon=True, fontsize=17)

        for i, ch in enumerate(ax.get_children()):
            try:
                txt=ch.get_text()
                if '13 TeV' in txt: ax.get_children()[i].set_text('')            
            except Exception as e:
                continue
                
        plt.savefig(f'Plots/DoubleMuonPaths_{args.outname}_{bins_sum}MeV.pdf', bbox_inches='tight')
        plt.close()











        frac=12*bins_sum*0.005
        fig, ax  = plt.subplots(figsize=[15,7]) 

        color1 = mplcolors.to_rgba('tab:blue')
        color1 = list(color1)
        color1[0]+=0.15
        color1[1]+=0.15
        color1[2]+=0.15
        color1[-1] = 0.9
        #print('color1 : ', color1)

        color2 = mplcolors.to_rgba('navy')
        color2 = list(color2)
        color2[-2] = 1
        color2[0]+=0.2
        color2[1]+=0.2
        #print('color2 : ', color2)

        commons.hist_from_heights(counts_low['HLT_Mu10_Barrel_L1HP11_IP6_v'], bins_dim, axis=ax, 
                                    color=color1, label='Mu10_Barrel_L1HP11_IP6',
                                    rasterized=True)
        commons.hist_from_heights(counts_low['HLT_Mu0_Barrel_v'], bins_dim, axis=ax, 
                                    color=color2, label='Mu0_Barrel',
                                    rasterized=True)
        commons.hist_from_heights(counts_low['HLT_Mu0_Barrel_L1HP11_v'], bins_dim, axis=ax, 
                                    color=color4, label='Mu0_Barrel_L1HP11', histtype='step', linewidth=3,
                                    rasterized=True)



        ax.set_yscale('log')
        ax.set_xlim(0.2, 8.5)
        ax.legend(frameon=True)
        ax.set_ylabel('Candidates / '+ str(int(1000*np.mean(np.diff(bins))))+' MeV')

        commons.hep.cms.label(data=True, label='Preliminary',                         
                                com=13.6, ax=ax, 
                                year=None, lumi=None)

        # commons.hep.cms.label(data=True, exp='',
        #                             year=None, ax=ax,
        #                             com=13.6)
        color3 = mplcolors.to_rgba("#df979e")
        color3 = list(color3)
        color3[-1] = 0.8
        #print('color3 : ', color3)


        ax.set_xlabel('$m_{\mu \mu}$ [GeV]', fontsize=28)
        plt.savefig(f'Plots/SingleMuonPaths_{args.outname}_SingleMuonPaths_{bins_sum}MeV.pdf', bbox_inches='tight')
        plt.close()