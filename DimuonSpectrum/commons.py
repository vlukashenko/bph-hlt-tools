import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
plt.style.use(hep.style.CMS)

def reduce_histogram(counts, edges, bins_to_sum):
    if len(counts) != len(edges) - 1:
        print(len(counts), len(edges))
        raise ValueError("Counts and edges arrays must have compatible lengths")

    if len(counts) % bins_to_sum != 0:
        print(len(counts), bins_to_sum)
        raise ValueError("Number of bins to sum must evenly divide the total number of bins")

    reduced_counts = []
    reduced_edges = []

    for i in range(0, len(counts), bins_to_sum):
        sum_counts = np.sum(counts[i:i + bins_to_sum])
        reduced_counts.append(sum_counts)
        reduced_edges.append(edges[i])
    
    reduced_edges.append(edges[-1])  # Add the last edge

    return np.array(reduced_counts), np.array(reduced_edges)


def create_axes(fig, split = 70, space_between = 2):
    size = 200
    split *= size/100
    split = int(split)
    
    space_between *= size/100
    space_between = int(space_between)
    
    ax  = plt.subplot2grid(shape = (1,size), loc = (0,0),
                           colspan = split, fig = fig)
    axp = plt.subplot2grid(shape = (1,size), loc = (0, split+space_between),
                           colspan = size-(split+space_between), fig = fig)
    
    #axp.get_shared_x_axes().join(axp, ax)
    #ax.set_xticklabels([])
    return ax, axp


def hist_from_heights(heights, bin_edges, axis=None, histtype='bar', join_n_bins=1, how='right', yerr=[], **kwargs):

    if join_n_bins>1:
        heights_tmp, bin_edges_tmp = rebin_data(heights, bin_edges, join_n_bins=join_n_bins, how='right')
    else: 
        heights_tmp, bin_edges_tmp = heights, bin_edges
        
    if histtype=='errorbar':
        if len(yerr)==0:
            yerr = np.sqrt(heights_tmp)
        bin_mean = (np.array(bin_edges_tmp[1:])+np.array(bin_edges_tmp[:-1]))/2
        width = np.array(bin_edges_tmp[1:])-np.array(bin_edges_tmp[:-1])
        if axis:
            fig = axis.errorbar(bin_mean, heights_tmp, xerr=width/2, yerr=yerr,**kwargs)
        else:
            print(kwargs)
            fig = plt.errorbar(bin_mean, heights_tmp, xerr=width/2, yerr=yerr, **kwargs)       

    if histtype=='bar':
        bin_mean = (bin_edges_tmp[1:]+bin_edges_tmp[:-1])/2
        width = np.array(bin_edges_tmp[1:])-np.array(bin_edges_tmp[:-1])
        if axis:
            fig = axis.bar(bin_mean, heights_tmp, width=width, **kwargs)
        else:
            fig = plt.bar(bin_mean, heights_tmp, width=width, **kwargs)
    
    elif histtype in ['step', 'stepfilled'] :
        x_pos = np.append(bin_edges_tmp, bin_edges_tmp[-1])
        y_pos = np.append(np.append(0.1, heights_tmp),0.1)
        if axis:
            if histtype=='step':
                fig = axis.step(x_pos, y_pos, where="pre" ,**kwargs)
            elif histtype=='stepfilled':
                fig = axis.fill_between(x_pos, y_pos, step="pre" ,**kwargs)
        else:
            if histtype=='step':
                fig = plt.step(x_pos, y_pos,where="pre" ,**kwargs)
            elif histtype=='stepfilled':
                fig = plt.fill_between(x_pos, y_pos,step="pre" ,**kwargs)
    return fig
