import pandas as pd
import numpy as np
import haloutils as htils

def restrict_to_subs(hpath,data):
    filename = hpath+'/analysis/minihalos/bg_subidsofhost.npy'
    subids = np.load(filename)
    mask = np.in1d(data['base_rsid'],subids)
    return data[mask]

def restrict_to_host(hpath,data):
    zoomid = htils.load_zoomid(hpath)
    mask = (data['base_rsid'] == zoomid)
    return data[mask]

def restrict_to_host_and_subs(hpath,data):
    zoomid = htils.load_zoomid(hpath)
    filename = hpath+'/analysis/minihalos/bg_subidsofhost.npy'
    subids = np.load(filename)
    m1 = np.in1d(data['base_rsid'],subids)
    m2 = (data['base_rsid'] == zoomid)
    mask = np.logical_or(m1,m2)
    return data[mask]

def restrict_to_outside_host(hpath,data):
    zoomid = htils.load_zoomid(hpath)
    filename = hpath+'/analysis/minihalos/bg_subidsofhost.npy'
    subids = np.load(filename)
    m1 = np.in1d(data['base_rsid'],subids)
    m2 = (data['base_rsid'] == zoomid)
    mask = np.logical_or(m1,m2)
    return data[~mask]

def get_processed_data(hpath,fileout,pandas):
    filename = hpath + "/analysis/" + fileout
    data = np.load(filename)
    #print "found:",len(data)
    if pandas:
        return pd.DataFrame(data)
    else:
        return data

def load_minihalos_imf(hpath,imf='kroupa',pandas=True):
    fileout = "all_lw_minihalo_array_"+imf+".npy"
    return get_processed_data(hpath,fileout,pandas)

def load_minihalos(hpath,pandas=True):
    fileout = "all_minihalo_array.npy"
    return get_processed_data(hpath,fileout,pandas)

def load_ricotti_minihalos(hpath,imf,pandas=True):
    fileout = "modified_ricotti_minihalo_array_"+str(imf)+".npy"
    print(fileout)
    return get_processed_data(hpath,fileout,pandas)

def load_firstgal(hpath,pandas=True):
    fileout = "all_firstgal_array.npy"
    return get_processed_data(hpath,fileout,pandas)

def load_minihalo_descendants(hpath,pandas=True):
    data = np.load(hpath+"/analysis/minihalo_descendants.npy")
    if pandas:
        return pd.DataFrame(data)
    else:
        return data

# include Tegmark 
def tegmark_mcut(z,T=2000.):
    h = 0.6711
    omega_m = 0.3125
    M = 1e8/h * (T/(1+z))**1.5 * (0.6*10/1.22/1.98e4)**1.5 * (18*3.14*3.14/178/omega_m)**0.5 #in solar masses
    return M


def plot_meannumber(ax,nbins,data,color,linestyle,label,plot_line_only=False):
    data = np.fliplr(data)
    median = np.median(data,axis=0)
    cumsum = np.cumsum(data,axis=1)
    stds = np.std(cumsum,axis=0)
    bin_red = (nbins[1:] + nbins[:-1])/2.
    bin_red = bin_red[::-1]
    #print bin_red
    binned_redshift = np.hstack([0,bin_red])
    y = np.median(cumsum,axis=0)    
    scales = np.linspace(0,1,25)
    fmin = y-stds; fmax = y+stds
    bin_red = np.hstack([bin_red,0])
    y = np.hstack([y,np.max(y)])
    stds = np.hstack([stds,stds[-1]])
    if not plot_line_only:
        #bin_red.insert(0,0)
        #y.insert(0,np.max(y))
        #stds.insert(0,stds[0])
    #    for fillmin,fillmax in zip(vec_res*np.max(fmin),vec_res*np.max(fmax)):
        for scales in np.flipud(scales):
            #ax.fill_between(bin_red,fmin,fmax,alpha=0.5,color=color, edgecolor=color)
            ax.fill_between(bin_red,y-scales*stds,y+scales*stds,alpha=0.05,color=color, edgecolor=None, linewidth=0.0)

    #print bin_red[-1],y[-1],stds[-1]
    zlast,nlast,stdlast = bin_red[-1],y[-1],stds[-1]
    ax.plot(bin_red,y,color=color,linestyle=linestyle,linewidth=3,label=label)
    return zlast,nlast,stdlast

def plot_scatter(ax,nbins,data,color):
    for row in data:
        bin_red = (nbins[1:] + nbins[:-1])/2.
        ax.plot(np.flipud(bin_red),np.cumsum(np.flipud(row)),linestyle='-',color=color,linewidth=1)

def get_number_of_halos_per_redshift(hpath,data,nbins):
    
    groupedbysnap = data.groupby('snap')
    #print groupedbysnap.indices.keys()
    redshifts = htils.get_z_snap(hpath,groupedbysnap.indices.keys())
    n_minihalos = [len(row) for row in groupedbysnap.indices.values()]
    x = np.array(redshifts,'float32')
    y = np.array(n_minihalos,'float32')
    
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    return nbins,sy