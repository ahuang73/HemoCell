import sys
sys.path.append("C:/Users/jdabo/Dropbox/STAGE/code/jon/")  # 'default'
import matplotlib
from matplotlib import pyplot as plt        
import numpy as np
import pickle
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches



SAVE_FIGS = False

#%%

def plotDifference(diff1, geo_slice, N_t):
    diff2 = -diff1
    # First create two color maps. One for the geometry contour and one with 
    # transparency for the overlaid RBC heat map

    # Color map for geometry
    N = 256
    geometry_color = (200.0, 200.0, 200.0)
    geometry_color = [c/N for c in geometry_color]
    geometry_color.append(1.0)  # no transparency
    white_c = [1.0, 1.0, 1.0, 1.0]
    red_c = [1.0, 0.0, 0.0, 1.0]
    blue_c = [0.0, 0.0, 1.0, 1.0]
    
    
    geometry_color_array = np.array([geometry_color, white_c])
    map_object = LinearSegmentedColormap.from_list(name='geo_map', colors=geometry_color_array)
    plt.register_cmap(cmap=map_object)
    
    # Create color map for normal RBCs
    c1 = np.ones(N)
    c1 = np.linspace(geometry_color[0], red_c[0], N)
    c2 = np.zeros(N)
    c2 = np.linspace(geometry_color[1], red_c[1], N)
    c3 = np.zeros(N)
    c3 = np.linspace(geometry_color[2], red_c[2], N)
    c4 = np.ones(N)
    c4 = np.linspace(0.0, 1.0, N)   # transparency
    color_array = np.array([c1, c2, c3, c4]).transpose()
    
    # create a colormap object and register
    map_object = LinearSegmentedColormap.from_list(name='RBC_colors', colors=color_array)
    plt.register_cmap(cmap=map_object)

    # Create color map for normal RBCs
    c1 = np.ones(N)
    c1 = np.linspace(geometry_color[0], blue_c[0], N)
    c2 = np.zeros(N)
    c2 = np.linspace(geometry_color[1], blue_c[1], N)
    c3 = np.zeros(N)
    c3 = np.linspace(geometry_color[2], blue_c[2], N)
    c4 = np.ones(N)
    c4 = np.linspace(0.0, 1.0, N)   # transparency
    color_array = np.array([c1, c2, c3, c4]).transpose()
    
    # create a colormap object and register
    map_object = LinearSegmentedColormap.from_list(name='RBC_stiff_colors', colors=color_array)
    plt.register_cmap(cmap=map_object)
    
    # Mirror both matrices
    geo_slice = np.array(geo_slice[::-1])
    diff1  = np.array(diff1[::-1])
    diff2  = np.array(diff2[::-1])

    # Divide through time to get fraction of time
    diff1 /= N_t;   diff2 /= N_t
    
    # Transform sideview data with log for log scaling
    ThreePlusLog = lambda v : 3 + np.log10(v) if v > 0.0 else 0.0
    diff1 = np.array([[ThreePlusLog(v) for v in row] for row in diff1])
    diff2 = np.array([[ThreePlusLog(v) for v in row] for row in diff2])
    
    # Plot
    FONTSIZE = 36
    
    matplotlib.rc('xtick', labelsize=FONTSIZE) 
    matplotlib.rc('ytick', labelsize=FONTSIZE) 
    fig, ax = plt.subplots(1, 1, figsize=(17,15))
    ax.imshow(geo_slice.transpose(), cmap='geo_map')
    fig2 = ax.imshow(diff1.transpose(), cmap='RBC_colors', label="Healthy")
    fig3 = ax.imshow(diff2.transpose(), cmap='RBC_stiff_colors', label="Stiffened")
    for fig in [fig3, fig2]:
        colorbar = plt.colorbar(fig, shrink=0.285)
        colorbar.set_ticks([0.0, 1.0, 2.0, 3.0])
        if fig == fig3:
            colorbar.ax.set_ylabel(r"Difference in avg density ($\frac{RBC \ count}{N_{measure}})$", rotation=270, labelpad=60, fontsize=FONTSIZE)        
        colorbar.set_ticklabels(["-0.001", "-0.01", "-0.1", "-1.0"])
    ax = plt.gca()
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    spatial_conv_factor = 0.5
    xlabels = ['%d' % (tick * spatial_conv_factor) for tick in xticks]
    ylabels = ['%d' % (tick * spatial_conv_factor) for tick in yticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)    
    ax.set_xlim([10, 400])
    ax.set_ylim([210, 0])
    plt.grid(False)
    # MANUAL LEGEND
    red_patch = mpatches.Patch(color='red', label='Healthy')
    blue_patch = mpatches.Patch(color='blue', label='Stiffened')
    plt.legend(handles=[red_patch, blue_patch], fontsize=FONTSIZE, loc='lower right')
    plt.xlabel(r"Position in $\mu$m", fontsize=FONTSIZE)
    plt.ylabel(r"Position in $\mu$m", fontsize=FONTSIZE)
    if SAVE_FIGS:
        matplotlib.use('Agg')
        plt.savefig("RBCs_pos_over_time")
    else:
        matplotlib.use('Qt5Agg')
        plt.show()


#%%        
    
CASE = "AR2"
ROOT = "L:/no_backup/"
RESULTS_PATH = ROOT + CASE + "/output/heatmap/"

boundaries_path = RESULTS_PATH + CASE + "_spatial_boundaries.pickle"
geo_slice_path = RESULTS_PATH + CASE + "_geometry_slice.pickle"
        
xmin, xmax, ymin, ymax, zmin, zmax = pickle.load(open(boundaries_path, "rb"), encoding="latin1")
#xmin, xmax, ymin, ymax, zmin, zmax = -1.5, 410.5, -1.5, 251.5, -1.5, 410.5

geo_slice = pickle.load(open(geo_slice_path, "rb"))


#%%
N_t = 1000

C1 = "AR2"
C2 = "AR2_stiff"
  
sideview_AR2  = np.loadtxt(open(ROOT + C1 + "/output/heatmap/" + C1 + "_sideview.csv", "rb"), \
                       delimiter=" ", skiprows=1)
sideview_AR2s = np.loadtxt(open(ROOT + C2 + "/output/heatmap/" + C2 + "_sideview.csv", "rb"), \
                       delimiter=" ", skiprows=1)

diff = np.zeros(sideview_AR2.shape)
I,J = sideview_AR2.shape
for i in range(I):
    for j in range(J):
        diff[i][j] = sideview_AR2[i][j] - sideview_AR2s[i][j]

plotDifference(diff, geo_slice, N_t)

#%%