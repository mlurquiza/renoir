
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors


HA2EV = 27.211396641308


def plot_eps (omegas, epsM_list):

    print("Ploting spectra...")
    print()

    tickfont = {'fontname':'DejaVu Sans'}
    labelfont = {'fontname':'DejaVu Sans'}

    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.xlabel('Energy [eV]', size=16, **labelfont)
    plt.ylabel(r'Im $\mathregular{\epsilon}_M$', size=16)

    #plt.xlim([min(omegas*HA2EV), max(omegas*HA2EV)])
    plt.xlim(173, 193)
    #plt.ylim(0, 3)

    plt.gcf().subplots_adjust(left=0.20)
    plt.gcf().subplots_adjust(bottom=0.15)

    idx=1
    for epsM in epsM_list:
        plt.plot(omegas*HA2EV, epsM.imag, label=idx, linewidth=2)
        idx+=1

    plt.tight_layout()
    plt.legend()

    plt.savefig("plot_eps.pdf", format='pdf')

    plt.show()




def plot_t1t2t3 (outfile_list, lo_list):

    for code, outfile in enumerate(outfile_list):

        print("Ploting spectra...")
        print()


        labelfont = {'fontname': 'DejaVu Sans'}
        fig, axs = plt.subplots(3, sharex=True, gridspec_kw={'hspace': 0})
        color_list = ['mediumblue', 'coral', 'teal', 'red', 'olive', 'darkmagenta', 'navy', 'yellow', 'orange', 'darkgray', 'yellowgreen', 'mediumaquamarine', 'black']

        nl_o = len(lo_list)
        evals_o = []
        evals_c = []
        t1  = []
        t2  = []
        t3  = []

        outfile.seek(0)
        for idx in range(nl_o):

            line = outfile.readline()
            line = outfile.readline()

            line = outfile.readline()
            tok  = line.strip().split()
            evals_o.append(float(tok[3])*HA2EV)

            line = outfile.readline()
            tok  = line.strip().split()
            nl_c = int(tok[3])

            line = outfile.readline()

            for lc in range(nl_c):
                line = outfile.readline()
                tok  = line.strip().split()
                if idx == 0:
                    evals_c.append(float(tok[0]))
                    t1.append(float(tok[1]))
                t2.append(float(tok[2]))
                t3.append(float(tok[3]))


        # Convert to arrays
        evals_c = np.array(evals_c)
        evals_o = np.array(evals_o)

        t1 = np.array(t1)
        t2 = np.array(t2).reshape((nl_o, nl_c))
        t3 = np.array(t3).reshape((nl_o, nl_c))


        axs[0].bar (evals_c, t1, width=0.1, align='center', color='silver')
        for idx in range(nl_o):
            name  = round(evals_o[idx], 2)
            color = color_list[idx]
            #color = plt.get_cmap('twilight')(idx / (len(lo_list) - 1))
            axs[1].scatter (evals_c, t2[idx,:], s=5, zorder=10-idx, color=color, alpha=0.8, label=str(name)+' eV')
            axs[2].plot    (evals_c, t3[idx,:], color=color, linewidth='1.0', label=str(name)+' eV')

        axs[1].legend(loc="upper left", fontsize=8)
        axs[2].legend  (loc="upper left", fontsize=8)

        axs[0].set_ylabel("|$t^{(1)}$|", size=14, **labelfont)
        axs[0].axis(ymin=0.0, ymax=np.max(t1)+np.max(t1)/10)
        axs[0].tick_params(axis='both', which='major', labelsize=12)

        axs[1].set_ylabel("|$t^{(2)}$|", size=14, **labelfont)
        axs[1].axis(ymin=0-np.max(t2)/10, ymax=np.max(t2)+np.max(t2)/10)
        axs[1].tick_params(axis='both', which='major', labelsize=12)

        axs[2].set_ylabel("|$t^{(3)}$|", size=14)
        axs[2].axis(ymin=0, ymax=np.max(t3)+np.max(t3)/10)
        axs[2].tick_params(axis='both', which='major', labelsize=12)

        print(np.max(t3))

        for ax in axs:
            ax.label_outer()

        plt.xlabel('Energy [eV]', size=12, **labelfont)
        plt.xlim([min(evals_c)-2, max(evals_c)+2])

        plt.tight_layout()
        plt.savefig(f"plot_t1t2t3_{code}.pdf", format='pdf', bbox_inches='tight')



def plot_rixs_vs_loss (wloss_list, win_list, rixs_list):

    print("Ploting spectra...")
    print()

    axisfont  = {'fontname':'DejaVu Sans'}
    fig = plt.figure(figsize=(20, len(win_list)*1.5))

    # one subplot per win in win_list:
    gs  = fig.add_gridspec(len(win_list), 1, hspace=0, wspace=0)
    axs = gs.subplots(sharex='col', sharey='row')

    colors = ['mediumblue', 'coral', 'firebrick', 'seagreen', 'darkorange', 'purple']
    ncolors = len(colors)

    for code, rixs  in enumerate(rixs_list):

        color = colors[code]
        n_win = len(win_list)
        idx = n_win
        for i in range(n_win):
            idx -= 1
            label = str(code) + r' -  $\omega_{in}$ = ' + "{:5.1f}".format(win_list[idx]+20.8)
            #label = str(code) + r' -  $\omega_{in}$ = ' + str(win_list[idx]) #str(win)
            axs[i].plot(wloss_list+2.4, rixs[idx,:].imag/np.max(rixs.imag), linewidth=2.0, color=color, label=label)
            #axs[i].plot(win_list[idx]+20.8-(wloss_list+2.4), rixs[idx,:].imag/np.max(rixs.imag), linewidth=2.0, color=color, label=label)
            axs[i].legend(loc="upper left", fontsize=16)


    ## Hide x labels and tick labels for all but bottom plot.
    for ax in axs:
        ax.label_outer()

    #plt.axis([ min(wloss_list), max(wloss_list), 0, 1.2 ])
    #plt.xlim([0, 16])
    plt.xlabel(r'$\mathregular{Energy \ loss \ [eV]}$', size=20, **axisfont)
    plt.tick_params(axis='x', which='major', labelsize=20)

    plt.savefig("plot_rixs_vs_loss.pdf", format='pdf')




def plot_rixs_map (wloss_list, win_list, rixs_list):

    for code, rixs in enumerate(rixs_list):

        print("Ploting spectra...")
        print()


        plt.figure(figsize=(14, 10))
        plt.gcf().subplots_adjust(left=0.15, bottom=0.15)

        axisfont = {'fontname': 'DejaVu Sans'}

        plt.pcolormesh(wloss_list, win_list, rixs.imag, norm=colors.LogNorm( vmin=np.max(rixs.imag)/1000, vmax=np.max(rixs.imag) ), cmap=plt.get_cmap('GnBu'))

        print(np.max(rixs.imag))

        cbar = plt.colorbar(extend='neither', spacing='proportional', orientation='vertical', shrink=1)
        cbar.ax.tick_params(labelsize=26)

        plt.xlim([min(wloss_list), max(wloss_list)])
        plt.ylim([min(win_list), max(win_list)])

        plt.xlabel(r'$\mathregular{Energy \ loss \ [eV]}$', size=26, **axisfont)
        plt.ylabel(r"$\mathregular{Excitation \ energy \ [eV]}$", size=26, **axisfont)
        plt.tick_params(axis='both', which='major', labelsize=26)

        plt.savefig(f"plot_rixs_map_{code}.png", format='png')


