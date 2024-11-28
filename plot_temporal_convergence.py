from glob import glob
import matplotlib as mpl
import  matplotlib.ticker as mplticker
import matplotlib.pyplot as plt
import matplotlib.transforms as mpltr
import numpy as np

## TO DISPLAY MATPLOTLIBRC FILE: mpl.matplotlib_fname()
#mpl.rcParams["text.usetex"]=True
mpl.rcParams['font.family']='serif'
#mpl.rcParams['mathtext.rm']='sans'
#mpl.rcParams['mathtext.it']='sans:italic'
#mpl.rcParams['mathtext.bf']='sans:bold'
#mpl.rcParams['mathtext.default']='sf'
mpl.rcParams['font.size']=12
mpl.rcParams['axes.titlesize']='large'
mpl.rcParams['axes.labelsize']='large'
mpl.rcParams['xtick.labelsize']='small'
mpl.rcParams['ytick.labelsize']='small'
mpl.rcParams['legend.fontsize']='small'
mpl.rcParams['figure.figsize']=(8,8)
mpl.rcParams['figure.dpi']=80

def _list_monitor_files():

    return glob('./conv_*.plt')

def read_monitor_file (filename):

    coords = np.loadtxt(filename, delimiter=',', skiprows=1, max_rows=1)
    types = _read_line(filename, 3)
    data = np.loadtxt(filename, delimiter=',', skiprows=4)

    return coords, data, types

def read_all_monitor_files():

    filelist = _list_monitor_files()
    nfiles = len(filelist)

    coords = np.zeros((nfiles, 3))

    # read first file
    coords[0,:], data, point_label = read_monitor_file(filelist[0])
    ntime, nvar = np.shape(data)
    data_temp = np.zeros((ntime, nvar, nfiles))
    lbl = []
    data_temp[:,:,0] = data
    lbl.append(point_label)

    for i, filepath in enumerate(filelist):
        if i > 0:
            coords[i,:], data_temp[:,:,i], point_label = \
                    read_monitor_file(filepath)
            lbl.append(point_label)

    return coords, data_temp, lbl

def _read_line(filename, line_number):
    with open(filename, 'r') as file:
        for i, line in enumerate(file,start=1):
            if i == line_number:
                return line.strip()

def plot(max_variation=1.e-2, plot_variation_relative=True):

    coords, data_temp, lbl = read_all_monitor_files()
    ntime, nvar, nfiles = np.shape(data_temp)

    data_full = data_temp.copy()

    if plot_variation_relative:
        data_full[1:,3:-1,:] = np.abs(1. - data_temp[1:,3:-1,:] \
                                         / data_temp[:-1,3:-1,:]) * 100.
    else:
        data_full[1:,3:-1,:] = np.abs(data_temp[1:,3:-1,:] \
                                      - data_temp[:-1,3:-1,:])

    data_full[1,:,:] = np.nan


    itr = data_full[:,0,:]
    time = data_full[:,1,:]
    dt = data_full[:,2,:]
    u = data_full[:,3,:]
    v = data_full[:,4,:]
    w = data_full[:,5,:]
    p = data_full[:,6,:]
    max_div = data_full[:,7,:]

    target_variation = max_variation * np.ones(np.size(time[:,0]))

    fields = (time, u, v, w, p)
    _figure_all_points(nfiles, fields, lbl,
                       target_variation, plot_variation_relative)
    _figure_turbines(nfiles, fields, lbl,
                     target_variation, plot_variation_relative)
    _figure_downstream(nfiles, fields, lbl,
                       target_variation, plot_variation_relative)
    _figure_outlet(nfiles, fields, lbl,
                   target_variation, plot_variation_relative)
    _figure_max_div(time, max_div)

def _figure_generic(plot_function, nfiles, fields_tuple, lbl,
                    target_variation, plot_variation_relative):

    time, u, v, w, p = fields_tuple

    fig = plt.figure(1, figsize=[8,6])
    fig.subplots_adjust(left=None,bottom=None,right=None,top=None,\
    	                wspace=0.25,hspace=0.25)
    ax1=fig.add_subplot(4,1,1)
    plot_function(ax1, time, u, lbl, log=plot_variation_relative)
    ax1.plot(time,target_variation,color='k')

    ax2=fig.add_subplot(4,1,2)
    plot_function(ax2, time, v, lbl, log=plot_variation_relative)
    ax2.plot(time,target_variation,color='k')

    ax3=fig.add_subplot(4,1,3)
    plot_function(ax3, time, w, lbl, log=plot_variation_relative)
    ax3.plot(time,target_variation,color='k')
    leg3=ax3.legend(loc=0,ncol=nfiles)

    ax4=fig.add_subplot(4,1,4)
    plot_function(ax4, time, p, lbl, log=plot_variation_relative)
    ax4.plot(time,target_variation,color='k')

    ax1.set_ylim((1.e-12,1.e2))
    ax2.set_ylim(ax1.get_ylim())
    ax3.set_ylim(ax1.get_ylim())
    ax4.set_ylim(ax1.get_ylim())
    ax4.set_xlim(0)
    ax4.xaxis.set_minor_formatter( mplticker.NullFormatter() )
    ax1.set_xlim(ax4.get_xlim())
    ax2.set_xlim(ax4.get_xlim())
    ax3.set_xlim(ax4.get_xlim())
    # ax1.set_xticks(ax4.get_xticks())
    # ax2.set_xticks(ax4.get_xticks())
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax3.set_xticklabels([])
    ax1.set_xlabel("")
    ax2.set_xlabel("")
    ax3.set_xlabel("")
    ax4.set_xlabel(r"Time")
    ax1.set_ylabel(r"u")
    ax2.set_ylabel(r"v")
    ax3.set_ylabel(r"w")
    ax4.set_ylabel(r"p")

    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(left=.1,bottom=.1,right=.95,top=.95)

    return ax1

def _figure_all_points(nfiles, fields, lbl,
                       target_variation, plot_variation_relative):
    ax1 = _figure_generic(_plot_all_monitor_points,
                               nfiles, fields, lbl,
                               target_variation, plot_variation_relative)
    if plot_variation_relative:
        ax1.set_title("All monitoring points - Relative variation (%)")
    else:
        ax1.set_title("All monitoring points - Absolute variation")

    plt.savefig(("convergence_full.png"),format="png")
    plt.savefig(("convergence_full.pdf"),format="pdf")
    plt.show()

def _figure_turbines(nfiles, fields, lbl,
                     target_variation, plot_variation_relative):
    ax1 = _figure_generic(_plot_turbines,
                               nfiles, fields, lbl,
                               target_variation, plot_variation_relative)
    if plot_variation_relative:
        ax1.set_title("Rotor centres - Relative variation (%)")
    else:
        ax1.set_title("Rotor centres - Absolute variation")

    plt.savefig(("convergence_turbines.png"),format="png")
    plt.savefig(("convergence_turbines.pdf"),format="pdf")
    plt.show()

def _figure_downstream(nfiles, fields, lbl,
                       target_variation, plot_variation_relative):
    ax1 = _figure_generic(_plot_downstream,
                               nfiles, fields, lbl,
                               target_variation, plot_variation_relative)
    if plot_variation_relative:
        ax1.set_title("2D Downstream - Relative variation (%)")
    else:
        ax1.set_title("2D Downstream - Absolute variation")

    plt.savefig(("convergence_downstream.png"),format="png")
    plt.savefig(("convergence_downstream.pdf"),format="pdf")
    plt.show()

def _figure_outlet(nfiles, fields, lbl,
                   target_variation, plot_variation_relative):
    ax1 = _figure_generic(_plot_outlet,
                               nfiles, fields, lbl,
                               target_variation, plot_variation_relative)
    if plot_variation_relative:
        ax1.set_title("Outlet - Relative variation (%)")
    else:
        ax1.set_title("Outlet - Absolute variation")

    plt.savefig(("convergence_outlet.png"),format="png")
    plt.savefig(("convergence_outlet.pdf"),format="pdf")
    plt.show()

def _figure_max_div(time, max_div):
    fig = plt.figure(1, figsize=[8,5])
    fig.subplots_adjust(left=None,bottom=None,right=None,top=None,\
    	                wspace=0.25,hspace=0.25)
    ax1=fig.add_subplot(1,1,1)
    ax1.semilogy(time[:,0], max_div[:,0], c='k',marker='None', ls='-')
    ax1.set_xlabel(r"Time")
    ax1.set_ylabel(r"Max DIV")
    # ax1.set_ylim((1.e-4, 1.e1))
    ax1.set_xlim(0)
    plt.savefig(("max_div.png"),format="png")
    plt.savefig(("max_div.pdf"),format="pdf")
    plt.show()

def _plot_monitor_class(ax, x, y, lbl, monitor_type='T', linestyle='-', 
                        log=False):
    colours = [
               [.8,0,0],
               [0,.8,0],
               [0,0,.6],
               [.6,.6,.6],
               [0,.9,.9],
               [.8,.8,0],
               [.8,0,.8]]

    _, nfiles = np.shape(x)

    for k in range(nfiles):
        point_type = lbl[k][0]
        point_number = int(lbl[k][1])
        if point_type == monitor_type:
            ax.semilogy(x[:,k],y[:,k],
                        c=colours[point_number],marker='None',
                        ls=linestyle,label=lbl[k])

def _plot_turbines(ax, x, y, lbl, linestyle='-', log=False):
    _plot_monitor_class(ax, x, y, lbl, 'T', linestyle, log)

def _plot_downstream(ax, x, y, lbl, linestyle='-', log=False):
    _plot_monitor_class(ax, x, y, lbl, 'D', linestyle, log)

def _plot_outlet(ax, x, y, lbl, linestyle='-', log=False):
    _plot_monitor_class(ax, x, y, lbl, 'O', linestyle, log)

def _plot_all_monitor_points(ax, x, y, lbl, log=False):

    _plot_turbines(ax, x, y, lbl, '-', log)
    _plot_downstream(ax, x, y, lbl, '-.', log)
    _plot_outlet(ax, x, y, lbl, '--', log)

if __name__ == "__main__":
    plot(1.e-2, True)
