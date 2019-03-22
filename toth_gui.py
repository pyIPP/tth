""" TRGUI

This GUI-based program prepares the input for transport codes

"""

__author__  = 'Giovanni Tardini (Tel. 1898)'
__version__ = '2.0'
__date__    = '06.03.2012'

import os, sys
sys.path.append('/afs/ipp/home/g/git/python/repository')
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
try:
    import Tkinter as tk
    import ttk, tkFont
except:
    import tkinter as tk
    from tkinter import ttk
import numpy as np
import matplotlib.pylab as plt
import tkhyper, sfdiff
import exec_toth, plot_toth
import dd_20180130, ww_20180130

sf  = dd_20180130.shotfile()

frc = '#b0d0b0'  # Frame, Notebook, Checkbutton
tbc = '#eaeaea'  # Toolbar, Button
figc = '#70a0c0' # Figure


class TOTH:


    def __init__(self):

        tothframe = tk.Tk()
        tothframe.title('TOT/TTH evaluation')
        tothframe.configure(bg=frc)

#        print(tkFont.families())
        toth_font = tkFont.Font(family='Century', size=15)
        tothframe.option_add("*Font", toth_font)

# Widgets style

        ypad = 2
        xpad = 1

        sty = ttk.Style()
        sty.theme_use('alt')
        sty.configure('TNotebook.Tab', background=frc, width=12)
        sty.configure('TFrame', background=frc)
        sty.configure('TCheckbutton', background=frc)
        sty.configure('TButton', background=tbc, width=10)
        sty.configure('TLabel', background=frc, relief='groove', width=18)
        sty.configure('TEntry', background=frc, width=8)

# Menubar

        menubar = tk.Menu(tothframe)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Run"     , command=self.exec_toth)
        filemenu.add_command(label="Plot"    , command=self.plot_toth)
        filemenu.add_command(label="Write SF", command=self.write_sf)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=sys.exit)
        helpmenu = tk.Menu(menubar, tearoff=0)
        helpmenu.add_command(label="About", command=self.about)
        menubar.add_cascade(label="File", menu=filemenu)
        menubar.add_cascade(label="Help", menu=helpmenu)
        tothframe.config(menu=menubar)

#=========
# Frames
#=========

        toolframe = ttk.Frame(tothframe)
        entframe  = ttk.Frame(tothframe)
        rbframe   = ttk.Frame(tothframe)

        separator = ttk.Separator(tothframe, orient="horizontal")

        for frame in toolframe, separator, entframe, rbframe:
            frame.pack(side=tk.TOP, fill=tk.BOTH, anchor=tk.W, pady=ypad)

# Toolbar
        locdir = os.path.dirname(os.path.realpath(__file__))
        playfig = tk.PhotoImage(file='%s/play.gif' %locdir)
        plotfig = tk.PhotoImage(file='%s/plot.gif' %locdir)
        sffig   = tk.PhotoImage(file='%s/save.gif' %locdir)
        exitfig = tk.PhotoImage(file='%s/exit.gif' %locdir)

        but_wid =55
        play_button = ttk.Button(toolframe, command=self.exec_toth, image=playfig, width=but_wid)
        plot_button = ttk.Button(toolframe, command=self.plot_toth, image=plotfig, width=but_wid)
        sf_button   = ttk.Button(toolframe, command=self.write_sf , image=sffig, width=but_wid)
        exit_button = ttk.Button(toolframe, command=sys.exit, image=exitfig, width=but_wid)

        for but in play_button, sf_button, plot_button, exit_button:
            but.pack(side=tk.LEFT, padx=xpad)

# Options 

        self.toth_d = {}
        user = os.getenv('USER')
        dic_init = {'shot': '28053', 't_fringe': 0.,\
             'ne_exp' : 'AUGD', 'ne_ed' : 0, \
             'equ_exp': 'AUGD', 'equ_ed': 0}
        for key in ('shot', 'ne_exp', 'ne_ed', 'equ_exp', 'equ_ed', 't_fringe'):
            locframe = ttk.Frame(entframe)
            locframe.pack(side=tk.TOP, anchor=tk.W, pady=ypad)
            ttk.Label(locframe, text=key).pack(side=tk.LEFT, anchor=tk.W, padx=xpad)
            self.toth_d[key] = ttk.Entry(locframe)
            self.toth_d[key].insert(0, dic_init[key])
            self.toth_d[key].pack(side=tk.LEFT, anchor=tk.W)

#---
        locframe = ttk.Frame(rbframe)
        locframe.pack(side=tk.TOP, anchor=tk.W, pady=ypad)
        ttk.Label(locframe, text='SF exp').pack(side=tk.LEFT, anchor=tk.W, padx=xpad)

        key = 'out_exp'
        self.toth_d[key] = tk.StringVar()
        self.toth_d[key].set(user)
        for val in (user, 'AUGD'):
            ttk.Radiobutton(locframe, variable=self.toth_d[key], \
                value=val, text=val).pack(side=tk.LEFT, anchor=tk.W)
#---
        locframe = ttk.Frame(rbframe)
        locframe.pack(side=tk.TOP, anchor=tk.W, pady=ypad)
        ttk.Label(locframe, text='NBI par').pack(side=tk.LEFT, anchor=tk.W, padx=xpad)

        key = 'NBIpar'
        self.toth_d[key] = tk.StringVar()
        self.toth_d[key].set('TRANSP 2012')
        for val in ('TRANSP 2012', 'FAFNER 1990'):
            ttk.Radiobutton(locframe, variable=self.toth_d[key], \
                value=val, text=val).pack(side=tk.LEFT, anchor=tk.W)
#---
        locframe = ttk.Frame(rbframe)
        locframe.pack(side=tk.TOP, anchor=tk.W, pady=ypad)
        ttk.Label(locframe, text='Equ diag').pack(side=tk.LEFT, anchor=tk.W, padx=xpad)

        key = 'equ_dia'
        self.toth_d[key] = tk.StringVar()
        self.toth_d[key].set('GQH')
        for val in ('GQH', 'FPG', 'IDG'):
           ttk.Radiobutton(locframe, variable=self.toth_d[key], \
                value=val, text=val).pack(side=tk.LEFT, anchor=tk.W)
#---
        locframe = ttk.Frame(rbframe)
        locframe.pack(side=tk.TOP, anchor=tk.W, pady=ypad)
        ttk.Label(locframe, text='ne diag').pack(side=tk.LEFT, anchor=tk.W, padx=xpad)

        key = 'ne_diag'
        self.toth_d[key] = tk.StringVar()
        self.toth_d[key].set('DCK')
        for val in ('DCK', 'DCN', 'DCR', 'DCP', 'DCS'):
           ttk.Radiobutton(locframe, variable=self.toth_d[key], \
                value=val, text=val).pack(side=tk.LEFT, anchor=tk.W)
#---
        locframe = ttk.Frame(rbframe)
        locframe.pack(side=tk.TOP, anchor=tk.W, pady=ypad)
        ttk.Label(locframe, text='ne sig').pack(side=tk.LEFT, anchor=tk.W, padx=xpad)

        key = 'ne_sig'
        self.toth_d[key] = tk.StringVar()
        self.toth_d[key].set('H-0')
        for val in ('H-0', 'H-1'):
           ttk.Radiobutton(locframe, variable=self.toth_d[key], \
                value=val, text=val).pack(side=tk.LEFT, anchor=tk.W)
#---
        locframe = ttk.Frame(rbframe)
        locframe.pack(side=tk.TOP, anchor=tk.W, pady=ypad)
        ttk.Label(locframe, text='Force TOT').pack(side=tk.LEFT, anchor=tk.W, padx=xpad)

        key = 'force_tot'
        self.toth_d[key] = tk.BooleanVar()
        self.toth_d[key].set(False)
        ttk.Checkbutton(locframe, variable=self.toth_d[key]).pack(side=tk.LEFT, anchor=tk.W)
#---
        locframe = ttk.Frame(rbframe)
        locframe.pack(side=tk.TOP, anchor=tk.W, pady=ypad)
        ttk.Label(locframe, text='Force TTH').pack(side=tk.LEFT, anchor=tk.W, padx=xpad)

        key = 'force_tth'
        self.toth_d[key] = tk.BooleanVar()
        self.toth_d[key].set(False)
        ttk.Checkbutton(locframe, variable=self.toth_d[key]).pack(side=tk.LEFT, anchor=tk.W)

#---
        locframe = ttk.Frame(rbframe)
        locframe.pack(side=tk.TOP, anchor=tk.W, pady=ypad)
        ttk.Label(locframe, text='Auto SF-write').pack(side=tk.LEFT, anchor=tk.W, padx=xpad)

        key = 'sf_write'
        self.toth_d[key] = tk.BooleanVar()
        self.toth_d[key].set(False)
        ttk.Checkbutton(locframe, variable=self.toth_d[key]).pack(side=tk.LEFT, anchor=tk.W)

        tothframe.mainloop()


    def about(self):

        mytext = 'Documentation at <a href="http://www.aug.ipp.mpg.de/~git/tot/index.html">TOT/TTH diagnostic homepage</a>'
        h = tkhyper.HyperlinkMessageBox("Help", mytext, "500x60")


    def exec_toth(self):

        toth_in = {}
        for key, val in self.toth_d.items():
            toth_in[key] = val.get()

        n_shots = np.atleast_1d(eval(toth_in['shot']))
        for nshot in n_shots:
            self.nshot = nshot
            self.toth = exec_toth.ex_toth(nshot, toth_in)
            if toth_in['sf_write']:
                self.write_sf()


    def plot_toth(self):

        if hasattr(self, 'toth'):
            plot_toth.plot_toth(self.nshot, self.toth)
        else:
            print('Run code before plotting')


    def write_sf(self):

        if hasattr(self, 'toth'):
            force_tot = self.toth_d['force_tot'].get()
            force_tth = self.toth_d['force_tth'].get()
            sfhdir = os.path.dirname(os.path.realpath(__file__))
            if 'time' in self.toth.tot.keys():
                exp_write = self.toth_d['out_exp'].get().strip()
# TOT
                if force_tot or sfdiff.sfdiff(self.nshot, exp_write, 'TOT', self.toth.tot):
                    status = ww_20180130.write_sf(self.nshot, self.toth.tot, sfhdir, \
                             'TOT', exp=exp_write)

                self.toth.tth['TOT_file']['expr']    = exp_write
                if sf.Open('TOT', self.nshot, experiment=exp_write):
                    print('\nTOT edition used for TTH: %d\n' %sf.edition)
                    self.toth.tth['TOT_file']['edition'] = sf.edition
                    sf.Close()
# TTH
                    if force_tth or sfdiff.sfdiff(self.nshot, exp_write, 'TTH', self.toth.tth):
                        status = ww_20180130.write_sf(self.nshot, self.toth.tth, sfhdir, \
                                 'TTH', exp=exp_write)
            else:
                print('No data to store')
        else:
            print('Run code before Storing shotfiles')


TOTH()
