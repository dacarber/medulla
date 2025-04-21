# imports
import uproot
import ast
from ROOT import TFile, TEfficiency, TH1D, TGraphAsymmErrors, RDataFrame, TCanvas
import numpy as np
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
plt.style.use('pi0ana.mplstyle')
import pandas as pd
import toml
import argparse
from array import array
from evaluation import *

def main(args):
    
    # Load variable configuration
    #var_cfg = toml.load(args.var_cfg)

    # Load ROOT file into pandas dataframe

    # Selected
    rf = uproot.open(args.in_file)
    pot = rf['events/mc/POT'].to_numpy()[0][0]
    sel_nu_tree = rf['events/mc/SelectedNu_PhaseCuts']
    sel_nu_df = sel_nu_tree.arrays(library='pd')
    sel_cos_tree = rf['events/mc/SelectedCos_PhaseCuts']
    sel_cos_df = sel_cos_tree.arrays(library='pd')
    sel_df = pd.concat([sel_nu_df, sel_cos_df])

    # Purity
    pur_tree = rf['events/mc/Purity_PhaseCuts']
    pur_df = pur_tree.arrays(library="pd")
    # Efficiency
    sig_tree = rf['events/mc/Signal_PhaseCuts']
    sig_df = sig_tree.arrays(library="pd")
    sig_df = sig_df[sig_df.category == 0] # only fiducialized signal events

    # Confusion
    '''
    eff_pid_confusion_tree = rf['events/mc/EffPIDConfusion_PhaseCuts']
    eff_pid_confusion_df = eff_pid_confusion_tree.arrays(library="pd")
    pur_pid_confusion_tree = rf['events/mc/PurPIDConfusion_PhaseCuts']
    pur_pid_confusion_df = pur_pid_confusion_tree.arrays(library="pd")
    eff_primary_confusion_tree = rf['events/mc/EffPrimaryConfusion_PhaseCuts']
    eff_primary_confusion_df = eff_primary_confusion_tree.arrays(library="pd")
    '''
    
    # Purity
    purity = calc_flat_purity(sel_df)
    print(f'Purity: {round(purity,2)}%')

    # Efficiency
    efficiency = calc_flat_efficiency(sig_df)
    print(f'Efficiency: {round(efficiency,2)}%')

    # As a function of each cut...
    cuts = {'Flash Cut':['flash_cut == 1'], 
            "Fiducial Cut":['flash_cut == 1','fiducial_cut == 1'], 
            "Topology Cut":['flash_cut == 1','fiducial_cut == 1','category_topology == 0'],
            
    calc_purity_by_cut(pur_df, cuts)
    calc_efficiency_by_cut(sig_df, cuts)
    
    # Calculate efficiency as function of...
    #plot_diff_pur_eff(sig_df, pot, 'muon_momentum_mag', var_cfg, [0,3200], 40, 'eff') # muon momentum
    #plot_diff_pur_eff(sig_df, pot, 'pi0_momentum_mag', var_cfg, [0,1000], 30, 'eff') # pi0 momentum

    # Confusion
    #plot_pid_confusion(eff_pid_confusion_df, pot, pur_or_eff='eff')
    #plot_pid_confusion(pur_pid_confusion_df, pot, pur_or_eff='pur')
    #plot_eff_primary_confusion(eff_primary_confusion_df, 2)    

def add_plot_labels(ax,pot, adj_y=0.045, title=str()):
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.025, #0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = xrange[0] + 0.025*(xrange[1] - xrange[0]) #0.025
    prelim_label = r'$\bf{ICARUS}$ Simulation Work-in-Progress'
    ax.text(x=usex, y=usey, s=prelim_label, fontsize=14, color='#d67a11')
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.02, 0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = (xrange[1] - xrange[0]) / 2
    ax.text(x=usex, y=usey, s=f'{title}', fontsize=14, color='black')
    yrange = ax.get_ylim()
    usey = yrange[1] + adj_y*(yrange[1] - yrange[0]) #0.02, 0.045 for confusion matrix
    xrange = ax.get_xlim()
    usex = xrange[1] - 0.025*(xrange[1] - xrange[0])
    mag = int(np.floor(np.log10(pot)))
    usepot = pot/10**mag
    s = f'{usepot:.2f}'+f'$\\times 10^{{{mag}}}$ POT'
    ax.text(x=usex, y=usey, s=s, fontsize=13, color='black', horizontalalignment='right')

def calc_flat_purity(sel_df):
    total_selected_events = len(sel_df)
    matched_selected_events = len(sel_df[sel_df.category == 0])
    return 100 * matched_selected_events / total_selected_events

def calc_flat_efficiency(sig_df):
    total_signal_events = len(sig_df)
    selected_signal_events = len(sig_df[sig_df['all_cut'] == 1])
    return 100 * selected_signal_events / total_signal_events

def plot_diff_pur_eff(df, pot, var, var_cfg, var_range, nbins, metric):
    bcs = []
    bxerr0s = []
    bxerr1s = []
    byerr0s = []
    byerr1s = []
    vals = []

    # filtering/cuts                                                                                      
    if var_range[0] == var_range[-1]:
        df = df
    else:
        df = df[(df[var] > var_range[0]) & (df[var] < var_range[1])]

    #df['var_q'] = pd.qcut(df[var], q=nbins)
    df['var_q'] = pd.cut(df[var], np.linspace(var_range[0], var_range[1], nbins+1))

    hedges = sorted([i.left for i in df.var_q.unique().tolist()])
    lastedge = max([i.right for i in df.var_q.unique().tolist()])
    hedges.append(lastedge)
    hedges = np.array(hedges)
    hpass = TH1D('hpass', '', nbins, hedges)
    htotal = TH1D('htotal', '', nbins, hedges)
    counts = []

    for i,(name,group) in enumerate(df.groupby('var_q')):

        bxerr0s.append(name.left)
        bcs.append(name.mid)
        bxerr1s.append(name.right)

        if metric == 'pur':
            vals.append(calc_flat_purity(group) / 100)
            hpass.SetBinContent(i+1, len(group[group['category'] == 0]))
        elif metric == 'eff':
            vals.append(calc_flat_efficiency(group) / 100)
            hpass.SetBinContent(i+1, len(group[group['all_cut'] == 1]))
        htotal.SetBinContent(i+1, len(group))
        counts.append(len(group))

    gr = TGraphAsymmErrors()
    gr.Divide(hpass, htotal, 'cl=0.683 b(1,1) mode')                                             
    for i in range(nbins):
        byerr0s.append(gr.GetErrorYlow(i))
        byerr1s.append(gr.GetErrorYhigh(i))

    fig, ax1 = plt.subplots(figsize=(10,6))
    #ax1.bar(x=hedges[:-1], height=counts, width=np.diff(hedges), align='edge', fc='C1', ec='black')
    ax1.hist(df[var], range=(var_range[0], var_range[1]), bins=nbins, color='C1', edgecolor='black')
    ax1.set_xlabel(var_cfg['variables'][var]['xlabel'])
    ax1.set_ylabel('Entries')
    ax2 = ax1.twinx()
    ax2.set_ylim([0,1])
    ax2.errorbar(bcs, vals, xerr=[np.array(bcs) - np.array(bxerr0s), np.array(bxerr1s) - np.array(bcs)], yerr=[np.array(byerr0s), np.array(byerr1s)], fmt='o', capsize=2)
    if metric == 'pur':
        ylabel = 'Purity'
    elif metric == 'eff':
        ylabel = 'Efficiency'
    ax2.set_ylabel(ylabel)
    ax2.spines['right'].set_color('C0')
    ax2.yaxis.label.set_color('C0')
    ax2.tick_params(axis='y', colors='C0')
    #plt.title(r'Signal $\nu_{\mu}$ CC $\pi^{0}$', fontsize=16)
    plottile=r'Signal $\nu_{\mu}$ CC $\pi^{0}$'
    add_plot_labels(ax1,pot, adj_y=0.030, title=plottile)
    plt.savefig(metric + '_vs_' + var + '.png')

def calc_purity_by_cut(pur_df, cuts):
    no_cut_purity = len(pur_df[pur_df.category == 0]) / len(pur_df)
    print(f'No Cut Purity: {round(100*no_cut_purity, 2)}%')
    for key, value in cuts.items():
        cond_string = str()
        for i,c in enumerate(value):
            cond_string += c
            if i != len(value) - 1: cond_string += ' & '

        sel_by_cut = pur_df.query(cond_string)
        cut_purity = len(sel_by_cut[sel_by_cut.category == 0]) / len(sel_by_cut)
        print(f'{key} Purity: {round(100*cut_purity, 2)}%')
        
def calc_efficiency_by_cut(sig_df, cuts):
    print(f'No Cut Efficiency: 100%')
    for key, value in cuts.items():
        cond_string = str()
        for i,c in enumerate(value):
            cond_string += c
            if i != len(value) - 1: cond_string += ' & '

        sel_by_cut = sig_df.query(cond_string)
        cut_efficiency = len(sel_by_cut[sel_by_cut.category == 0]) / len(sig_df)
        print(f'{key} Efficiency: {round(100*cut_efficiency, 2)}%')

def plot_pid_confusion(pid_confusion_df, pot, pur_or_eff='eff'):
    # Start from true particles for efficiency; reco particles for purity
    if(pur_or_eff == 'eff'): 
        _from = 'true'
        _to = 'reco'
    if(pur_or_eff == 'pur'): 
        _from = 'reco'
        _to = 'true'

    # Get particles
    photons = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 0]
    photon_pids_from = [0] * len(photons)
    photon_pids_to = photons[f'{_to}_pid'].to_list()

    electrons = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 1]
    electron_pids_from = [1] * len(electrons)
    electron_pids_to = electrons[f'{_to}_pid'].to_list()

    muons = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 2]
    muon_pids_from = [2] * len(muons)
    muon_pids_to = muons[f'{_to}_pid'].to_list()

    pions = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 3]
    pion_pids_from = [3] * len(pions)
    pion_pids_to = pions[f'{_to}_pid'].to_list()

    protons = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 4]
    proton_pids_from = [4] * len(protons)
    proton_pids_to = protons[f'{_to}_pid'].to_list()
    
    kaons = pid_confusion_df[pid_confusion_df[f'{_from}_pid'] == 5]
    kaon_pids_from = [5] * len(kaons)
    kaon_pids_to = kaons[f'{_to}_pid'].to_list()

    # One flat list for all particles
    pids_from = photon_pids_from + electron_pids_from + muon_pids_from + pion_pids_from + proton_pids_from + kaon_pids_from
    pids_to = photon_pids_to + electron_pids_to + muon_pids_to + pion_pids_to + proton_pids_to + kaon_pids_to

    # Set labels
    label_keys = np.unique(pids_from)
    label_keys = np.append(label_keys, 3) # show pions
    label_keys = sorted(label_keys)
    label_values = []
    if 0 in label_keys: label_values.append('Photon')
    if 1 in label_keys: label_values.append('Electron')
    if 2 in label_keys: label_values.append('Muon')
    #if 3 in label_keys: label_values.append('Pion')
    label_values.append('Pion') # show pions
    if 4 in label_keys: label_values.append('Proton')
    if 5 in label_keys: label_values.append('Kaon')

    # Plot
    cm = confusion_matrix(np.array(pids_from), np.array(pids_to), normalize='true', labels=label_keys)
    cm_counts = confusion_matrix(np.array(pids_from), np.array(pids_to), labels=label_keys)
    fig, ax = plt.subplots(figsize=(8, 8))

    im = heatmap(cm.T * 100,
                 label_values,
                 label_values,
                 ax=ax, cmap="Blues")
    annotate_heatmap(im, unc=cm_counts.T, valfmt="{:.2f}% \n({})", fontsize=14)
    plt.setp(ax.get_yticklabels(), rotation=90, ha="center",rotation_mode="anchor")
    ax.tick_params(axis='both', which='major', pad=15, bottom=False, top=False, left=False, right=False, labelsize=14)
    ax.tick_params(axis="x", pad=1)
    if(pur_or_eff == 'eff'): ax.set_xlabel('True Signal Particles', fontsize=15, labelpad=10)
    if(pur_or_eff == 'pur'): ax.set_xlabel('Selected Reco Particles', fontsize=15, labelpad=10)
    ax.xaxis.set_label_position('bottom')
    if(pur_or_eff == 'eff'): ax.set_ylabel('Matched Reco Particles', fontsize=15)
    if(pur_or_eff == 'pur'): ax.set_ylabel('Matched True Particles', fontsize=15)
    add_plot_labels(ax,pot, adj_y=0.045, title=str())
    plt.savefig(f'{pur_or_eff}_pid_confusion.png')

def plot_eff_primary_confusion(eff_primary_confusion_df, true_particle_pid):
    true_particle_df = eff_primary_confusion_df[eff_primary_confusion_df.true_pid == true_particle_pid]
    true_particle_true_primaries = true_particle_df.true_primary.to_list()
    true_particle_reco_primaries = true_particle_df.reco_primary.to_list()

    cm = confusion_matrix(np.array(true_particle_true_primaries), np.array(true_particle_reco_primaries), normalize='true', labels=[0,1])
    cm_counts = confusion_matrix(np.array(true_particle_true_primaries), np.array(true_particle_reco_primaries), labels=[0,1])
    fig, ax = plt.subplots(figsize=(8, 8))

    im = heatmap(cm.T * 100,
                 ['Nonprimary', 'Primary'],
                 ['Nonprimary', 'Primary'],
                 ax=ax, cmap="Blues", interpolation='nearest')
    annotate_heatmap(im, unc=cm_counts.T, valfmt="{:.2f}% \n({})", fontsize=14)

    plt.setp(ax.get_yticklabels(), rotation=90, ha="center",rotation_mode="anchor")
    ax.tick_params(axis='both', which='major', pad=15, bottom=False, top=False, left=False, right=False, labelsize=14)
    ax.tick_params(axis="x", pad=1)
    ax.set_xlabel('True Signal Particles', fontsize=15, labelpad=10)
    ax.xaxis.set_label_position('top')
    ax.set_ylabel('Matched Reco Particles', fontsize=15)
    plt.show()
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', required=True)
    parser.add_argument('--var_cfg', required=False)
    args = parser.parse_args()
    main(args)