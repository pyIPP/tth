import numpy as np
import fconf
import Tkinter as tk
import ttk
import sfh_20170404 # leave this older version
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

sfhr = sfh_20170404.SFH_READ

col = 2*['#c00000', '#00c000', '#0000c0', '#b0b000', '#b000b0', '#00b0b0']
fsize = 8
titsize=10
lblsize=10

tth_laws = ( \
    'ITERL-89P(tot)', \
    'ITERL-96P(th)', \
    'ITERH92PY:tau_E', \
    'ITERH-92P(y, tot)', \
    'ITERH-92P(y, th)', \
    'ITERH-93P(th)', \
    'ITERH-98P(y, th)', \
    'ITERH-98P(y, th, 2)', \
    'ESGB, McDonald', \
    'CORDEY05, NF 45', \
    'KARDAUN, IAEA 2006' \
#                'KARDAUN-LANG, IAEA 2012, EX/P4-01' \
                 )

def fig_plot(frame, dic, nshot, sig_list, y_ticks=[], col='#c00000', n_rows=5, n_cols=7):

    figframe  = tk.Frame(frame)
    figframe.pack(side=tk.TOP, fill=tk.X)
    figframe.pack_propagate(0)

    fig = Figure(figsize=(8., 8.7), dpi=100)
    can = FigureCanvasTkAgg(fig, master=frame)
    can._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    fig.subplots_adjust(left=0.1, bottom=0.06, right=0.95, top=0.92, hspace=0)
    fig.text(.5, .95, '#%d' %nshot, ha='center')

    xmajorLocator = MultipleLocator(2)
    xmajorFormatter = FormatStrFormatter('%d')
    ymajorFormatter = FormatStrFormatter('%3.0e')

    n_plots = len(sig_list)
    n_rows = int(np.sqrt(n_plots)) + 1
    n_cols = int(n_plots/n_rows + 0.5) + 1
    jrow = 0
    for jx, sig in enumerate(sig_list):
        ax = fig.add_subplot(n_rows, n_cols, jx+1)
        ax.plot(dic['time'], dic[sig], color=col)
        ax.text(.5, .9, sig, ha='center', transform=ax.transAxes, fontsize=titsize)

# x
        ax.xaxis.set_major_locator(xmajorLocator)
        ax.xaxis.set_major_formatter(xmajorFormatter)
# y
        ax.yaxis.set_major_formatter(ymajorFormatter)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fsize) 

        if (jx % n_cols == 0):
            jrow += 1
        if jrow < n_rows:
            ax.set_xticklabels([])
        else:
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(fsize) 
            ax.set_xlabel('Time [s]', fontsize=lblsize)

    ax.tick_params(which='major', length=4, width=0.5)

    can.mpl_connect('button_press_event', fconf.on_click)
    toolbar = NavigationToolbar2TkAgg(can, frame)
    toolbar.update()


def fig_plot2(frame, time, sgr, nshot, col='#c00000', n_rows=5, n_cols=7, ymax=2):

    figframe  = tk.Frame(frame)
    figframe.pack(side=tk.TOP, fill=tk.X)

    fig = Figure(figsize=(8., 8.7), dpi=100)
    can = FigureCanvasTkAgg(fig, master=frame)
    can._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    fig.subplots_adjust(left=0.1, bottom=0.06, right=0.95, top=0.92, hspace=0)
    fig.text(.5, .95, '#%d' %nshot, ha='center')

    print sgr.shape
    n_sig = sgr.shape[1]
    for jx in range(n_sig):
        ax = fig.add_subplot(n_rows, n_cols, jx+1)
        ax.plot(time, sgr[:, jx], color=col)
        ax.plot([0, 10], [1, 1], 'k-')
        ax.set_ylabel(tth_laws[jx], fontsize=fsize)
        ax.set_ylim([0, ymax])
        if int(jx/n_cols) == n_rows-1:
            ax.set_xlabel('Time [s]', fontsize=fsize)
    can.mpl_connect('button_press_event', fconf.on_click)
    toolbar = NavigationToolbar2TkAgg(can, frame)
    toolbar.update()


def plot_toth(nshot, toth):

        if __name__ == '__main__':
            myframe = tk.Tk()
        else:
            myframe = tk.Toplevel()

        myframe.geometry('1200x960')
        myframe.title('Time traces')
        nb = ttk.Notebook(myframe, name='nb')
        nb.pack(side=tk.TOP, fill=tk.X)

        totframe = tk.Frame(nb)
        tthframe = tk.Frame(nb)
        hf_frame = tk.Frame(nb)
        lh_frame = tk.Frame(nb)
        txt = {totframe: 'TOT signals', tthframe: 'TTH signals', \
               hf_frame: 'H factors'  , lh_frame: 'L2H factors'}
        for frame in totframe, tthframe, hf_frame, lh_frame:
            tk.Label(frame, text=txt[frame]).pack(side=tk.TOP)
            nb.add(frame, text=txt[frame])

# TOT signals

        tot_list = sfhr('TOT').sig
        fig_plot(totframe, toth.tot, nshot, tot_list, col=col[0])

# TTH signals

        tth_list = sfhr('TTH').sig
        fig_plot(tthframe, toth.tth, nshot, tth_list, col=col[0])

# H-factors

        fig_plot2(hf_frame, toth.tth['time'], toth.tth['H/L-facs'], nshot, col=col[0], n_rows=4, n_cols=3)

# L2H-factors

        fig_plot2(lh_frame, toth.tth['time'], toth.tth['L2H_facs'], nshot, col=col[0], n_rows=4, n_cols=3, ymax=3)

        myframe.mainloop()
