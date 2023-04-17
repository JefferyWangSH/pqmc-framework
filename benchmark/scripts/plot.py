import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import os


if "__main__":

    # font convention
    label_font = {'family': 'Times New Roman', 'weight': 'regular', 'size': 20,}
    title_font = {'family': 'Times New Roman', 'weight': 'regular', 'size': 20,}
    legend_font = {'family': 'Times New Roman', 'weight': 'regular', 'size': 18,}

    mpl.rcParams.update({ "text.usetex": True, "font.family": "serif", "font.serif": ["Times"],})

    dataset_double_occu = np.loadtxt( "../double_occu.dat", usecols=(0,1,2,3,4), skiprows=1 )
    dataset_kinetic_energy = np.loadtxt( "../kinetic_energy.dat", usecols=(0,1,2,3,4), skiprows=1 )
    kinetic_energy_u, kinetic_energy_pqmc_mean, kinetic_energy_pqmc_err, kinetic_energy_dqmc_mean, kinetic_energy_dqmc_err = dataset_kinetic_energy.transpose()
    double_occu_u, double_occu_pqmc_mean, double_occu_pqmc_err, double_occu_dqmc_mean, double_occu_dqmc_err = dataset_double_occu.transpose()
    
    fig, ax = plt.subplots( 1, 1, figsize=(8,6) )

    ax.tick_params(axis='both', which='major', direction='in', length=5.0, width=1.5, labelsize=20.0)
    frame_size = 1.5
    ax.spines['left'].set_linewidth(frame_size)
    ax.spines['right'].set_linewidth(frame_size)
    ax.spines['top'].set_linewidth(frame_size)
    ax.spines['bottom'].set_linewidth(frame_size)

    # ax.errorbar( double_occu_u, double_occu_pqmc_mean, double_occu_pqmc_err, fmt="-o", linewidth=2, markersize=8, markeredgewidth=2, markerfacecolor="none", label=r"PQMC" )
    # ax.errorbar( double_occu_u, double_occu_dqmc_mean, double_occu_dqmc_err, fmt="-o", linewidth=2, markersize=8, markeredgewidth=2, markerfacecolor="none", label=r"DQMC" )
    # ax.set_ylabel( r"$Double\ Occupation\ D$", label_font, labelpad=10 )

    # ax.errorbar( kinetic_energy_u, kinetic_energy_pqmc_mean, kinetic_energy_pqmc_err, fmt="-o", linewidth=2, markersize=8, markeredgewidth=2, markerfacecolor="none", label=r"PQMC" )
    # ax.errorbar( kinetic_energy_u, kinetic_energy_dqmc_mean, kinetic_energy_dqmc_err, fmt="-o", linewidth=2, markersize=8, markeredgewidth=2, markerfacecolor="none", label=r"DQMC" )
    # ax.set_ylabel( r"$E_k/(Nt)$", label_font, labelpad=10 )

    u = double_occu_u
    egs_pqmc_mean = u * double_occu_pqmc_mean + kinetic_energy_pqmc_mean
    egs_pqmc_err = u * double_occu_pqmc_err + kinetic_energy_pqmc_err
    egs_dqmc_mean = u * double_occu_dqmc_mean + kinetic_energy_dqmc_mean
    egs_dqmc_err = u * double_occu_dqmc_err + kinetic_energy_dqmc_err
    ax.errorbar( u, egs_pqmc_mean, egs_pqmc_err, fmt="-o", linewidth=2, markersize=8, markeredgewidth=2, markerfacecolor="none", label=r"PQMC" )
    ax.errorbar( u, egs_dqmc_mean, egs_dqmc_err, fmt="-o", linewidth=2, markersize=8, markeredgewidth=2, markerfacecolor="none", label=r"DQMC" )
    ax.set_ylabel( r"$E/(Nt)$", label_font, labelpad=10 )

    ax.set_xlabel( r"$U$", label_font, labelpad=10 )
    ax.set_xticks( [1.0, 2.0, 3.0, 4.0] )
    ax.xaxis.set_major_formatter( FormatStrFormatter('%.1f') )
    ax.set_title( r"$L=4,\quad \beta=2\Theta=10.0, \quad \Delta\tau=0.05$", title_font, pad=10 )

    plt.legend( loc='best',
                labelspacing=0.4, markerfirst=True,
                frameon=True, fancybox=False, edgecolor='black', handlelength=1.4, 
                ncol=1, columnspacing=0.8, prop=legend_font )
    fig.tight_layout()
    plt.savefig( "../egs.pdf", format="pdf", bbox_inches="tight", dpi=240 )
    plt.show()