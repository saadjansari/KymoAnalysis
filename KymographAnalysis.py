#!/usr/bin/env python

import os, pdb, sys
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import seaborn as sns
import math, random
import glob, yaml, copy, shutil

from src.Strain import *

"""
Name: KymographAnalysis.py
Description: Parses and combines tracks from multiple kymographs for mass analysis 
"""


class KymographAnalysis:
    def __init__(self):

        self.cwd = os.getcwd()
        # Read config file
        with open("config.yaml") as f:
            self.config = yaml.load(f, Loader=yaml.CLoader)

        self.InitStrains()
        self.Analyze()

    # InitStrains {{{
    def InitStrains(self):

        # Initialize strains with track files
        self.strains = []
        # Get filenames for each strain
        for strain in self.config["strains"]:
            trackpaths = []
            for fpath in strain["path"]:
                trackpaths += glob.glob(fpath)

            # Initialize
            cstrain = Strain(trackpaths, label=strain["type"])
            cstrain.color = tuple(np.array(strain["color"]) / 255)
            self.strains += [cstrain]

        # Use Kmeans classfication if required
        if "useBipolarKmeansLabel" in self.config.keys():
            if self.config["useBipolarKmeansLabel"]:
                for strain, strain_c in zip(self.strains, self.config["strains"]):
                    strain.TrimUsingKmeansLabel(kmean_label=strain_c["kmean_label"])

    # }}}

    # Analyze {{{
    def Analyze(self):

        # Initialize graphing directory
        if "saveName" in self.config.keys():
            gdir = os.path.join(self.cwd, self.config["saveName"])
        else:
            gdir = os.path.join(self.cwd, "result")
        if os.path.exists(gdir):
            shutil.rmtree(gdir, ignore_errors=True)
        os.mkdir(gdir)
        os.chdir(gdir)

        # Analyze by groups
        if self.config["analyzeByLength"] is True:
            for group in self.config["analyzeGroups"]:

                # CD to directory
                ggdir = os.path.join(gdir, group["type"])
                os.mkdir(ggdir)
                os.chdir(ggdir)

                print("Analyzing {0}".format(group["type"]))
                strains = copy.deepcopy(self.strains)

                # Get tracks that match this spindle length
                for strain in strains:
                    strain.GetTracks(spindle_length=group["length"])

                # pdb.set_trace()
                if self.config["analyzeSpindleIntensity"] is True:
                    self.GraphSpindleIntensity(
                        strains, lrange=group["length"], gname=group["type"]
                    )
                self.Graph(strains, gname=group["type"])

        else:

            if self.config["analyzeSPBAssociatedTracks"] == 1:
                for strain in self.strains:
                    strain.TossFarTracks(self.config["SPBRegion"])

            if self.config["analyzeSPBAssociatedTracks"] == 2:
                for strain in self.strains:
                    strain.TossCloseTracks(self.config["SPBRegion"])

            # Get all tracks
            for strain in self.strains:
                strain.GetTracks()

            if self.config["analyzeSpindleIntensity"] is True:
                self.GraphSpindleIntensity(self.strains)
            self.Graph(self.strains)

        os.chdir(self.cwd)

    # }}}

    # Graph {{{
    def Graph(self, strains, gname=None):
        # Graph useful properties

        plt.rcParams.update({"font.size": 14})
        plt.rc("legend", fontsize=12)

        self.PlotTracksByState(k=1000)
        self.PlotAllTracks()
        # self.GraphStrain_EventsPerMinutePerCellViolin()
        # self.GraphStrain_SwitchFrequencyPerCellViolin()
        # self.GraphStrain_SwitchFrequencyPerCell()
        self.GraphStrain_FractionMovement()
        self.GraphStrain_SwitchFrequency2()
        self.GraphStrain_EventsPerMinute2()
        # self.GraphStrain_EventsPerMinutePerCell()
        self.GraphHistComparison()
        # self.GraphStrainMedianValues()
        # self.GraphStrain_EventsPerMinutePerCellRaw()
        # self.GraphStrain_SwitchFrequencyPerCellRaw()

        # Scatter Plots
        graphscatvars = (
            ['Run displacement','Intensity','nm','AU','scatter_intensity_runlength.pdf'],
            # ['Velocity','Intensity',r'$\mu$m/min','AU','scatter_intensity_velocity.pdf'],
            # ['Lifetime','Intensity','min','AU','scatter_intensity_lifetime.pdf'],
            # [
                # "Run length",
                # "Velocity",
                # r"$\mu$" + "m",
                # "nm/s",
                # "scatter_velocity_runlength.pdf",
            # ],
            # [
                # "Run length",
                # "Lifetime",
                # r"$\mu$" + "m",
                # "min",
                # "scatter_lifetime_runlength.pdf",
            # ],
            ['Velocity','Intensity','nm/s','AU','scatter_intensity_velocity.pdf'],
            # ["Velocity", "Lifetime", "nm/s", "min", "scatter_lifetime_velocity.pdf"],
            # [
                # "Run length",
                # "Average distance from SPB",
                # r"$\mu$" + "m",
                # r"$\mu$m",
                # "scatter_avgSPBdistance_runlength.pdf",
            # ],
            # [
                # "Velocity",
                # "Average distance from SPB",
                # "nm/s",
                # r"$\mu$m",
                # "scatter_avgSPBdistance_velocity.pdf",
            # ],
            # [
                # "Lifetime",
                # "Average distance from SPB",
                # "min",
                # r"$\mu$m",
                # "scatter_avgSPBdistance_lifetime.pdf",
            # ],
            ['Lifetime','Intensity','s','AU','scatter_intensity_lifetime.pdf'],
        )
        for x,y,xunit,yunit,figname in graphscatvars:
            self.GraphStrainScatter( strains,x,y,xlab=x,ylab=y,xunit=xunit,yunit=yunit,figname=figname )

        # self.GraphStrain_EventsPerMinute()
        # self.GraphStrain_AvgStartEnd()
        # self.GraphStrain_StateTimes()
        # self.GraphStrain_SwitchCounts()
        # self.GraphStrain_StateSwitchMatrix()

    # }}}

    # GraphHistComparison {{{
    def GraphHistComparison(self):
        def plot_median_special(ax, xloc, rel_height, col):

            (ybottom, ytop) = ax.get_ylim()
            ax.plot(
                [xloc, xloc],
                [ybottom, rel_height * ytop],
                color=col,
                linewidth=1.5,
                alpha=0.3,
                solid_capstyle="round",
            )
            ax.plot(
                [xloc],
                [rel_height * ytop],
                marker="d",
                color=col,
                alpha=0.6,
                markersize=6,
            )
            return ax

        graphhistvars = (
            [
                "GetRunLengths",
                "Run displacement",
                "Count",
                "nm",
                "strain_runlength.pdf",
            ],
            ["GetVelocities", "Velocity", "Count", "nm/s", "strain_velocity.pdf"],
            ["GetLifetimes", "Lifetime", "Count", "s", "strain_lifetime.pdf"],
            [
                "GetAverageDistances",
                "Average distance from SPB",
                "Count",
                "nm",
                "strain_avg_pos.pdf",
            ],
        )

        # Special x limits
        # xmaxes = {
        # 'Run displacement': 1.6,
        # 'Velocity': 100.0,
        # 'Lifetime': 1.6,
        # 'Average distance from SPB': 8.0,
        # }
        # if self.config['paperFigure'] == 5:
        # xmaxes['Velocity'] = 60.0
        xmaxes = {
            "Run displacement": 1600,
            "Velocity": 100.0,
            "Lifetime": 100,
            "Average distance from SPB": 8000,
        }
        ymax_scaling = {
            "Run displacement": 2000,
            "Velocity": 50,
            "Lifetime": 100,
            "Average distance from SPB": 10000,
        }
        if self.config["paperFigure"] == 5:
            xmaxes["Velocity"] = 60.0

        nStrain = len(self.strains)

        for fcn, xlab, ylab, unit, figname in graphhistvars:

            # Make a figure. Two axes (one for poleward, one for antipoleward)
            fig, ax = plt.subplots(figsize=(6, 3))
            cols1 = [cstrain.color for cstrain in self.strains]
            cols2 = [cstrain.color for cstrain in self.strains]
            # cols1 = [[68, 111, 200],[220, 95, 60]]
            # cols1 = [tuple(np.array(x)/255) for x in cols1]
            # cols2 = [[68, 111, 200],[220, 95, 60]]
            # cols2 = [tuple(np.array(x)/255) for x in cols2]

            # list for medians
            medians = {}

            original_stdout = (
                sys.stdout
            )  # Save a reference to the original standard output

            # Save to stats
            print(os.getcwd())
            with open("stats.txt", "a") as f:
                sys.stdout = f  # Change the standard output to the file we created.
                print("-" * 30)
                print("\nParameter = {0}".format(xlab))
            sys.stdout = original_stdout  # Reset the standard output to its original value            # Display

            # Make histograms for each strain
            for strain, col1, col2 in zip(self.strains, cols1, cols2):
                # Get data
                funcData = getattr(strain, fcn)
                dataPAP = funcData()
                dataAll = np.hstack((-1 * np.array(dataPAP[0]), dataPAP[1]))

                # bins and histogram
                nbins = 16
                bins = np.linspace(-1 * xmaxes[xlab], xmaxes[xlab], nbins + 1)
                # ax.hist( dataAll,  bins, density=True, edgecolor='k', alpha=0.6, color = col, label='{0} (N={1})'.format(strain.label, len(dataAll)))
                if self.config["paperFigure"] == 7:
                    print("Skip WT histograms for TD cell")
                else:
                    _, _, patches = ax.hist(
                        dataAll,
                        bins,
                        density=True,
                        edgecolor="white",
                        linewidth=1.0,
                        alpha=0.6,
                        color=col2,
                    )
                    for i in range(0, int(nbins / 2)):
                        patches[i].set_facecolor(col1)
                        patches[i].set_hatch("////")

                # Draw y-axis in middle
                ax.axvline(x=0, c="black", lw=1.5)

                # Add medians info
                medians[strain.label] = dict(
                    zip(
                        ["P", "AP"],
                        [fac * np.median(db) for db, fac in zip(dataPAP, [-1, 1])],
                    )
                )

                # ax.hist( [],  bins, edgecolor='white', linewidth=1.0, alpha=0.6, color = col2, label='{0} (N={1})'.format(strain.label, len(dataAll)))
                ax.hist(
                    [],
                    bins,
                    edgecolor="white",
                    linewidth=1.0,
                    alpha=0.6,
                    color=col2,
                    label="{0}".format(strain.label),
                )

                # Print Info
                # Save to stats
                with open("stats.txt", "a") as f:
                    sys.stdout = f  # Change the standard output to the file we created.
                    print("Strain: {0}".format(strain.label))
                    print("Poleward:")
                    print("\tN = {0}".format(len(dataPAP[0])))
                    print("\tMedian = {0:.3f} {1}".format(np.median(dataPAP[0]), unit))
                    print("\tMean = {0:.3f} {1}".format(np.mean(dataPAP[0]), unit))
                    print(
                        "\tStandard Dev = {0:.3f} {1}".format(np.std(dataPAP[0]), unit)
                    )
                    print(
                        "\tStandard Error = {0:.3f} {1}".format(
                            np.std(dataPAP[0]) / np.sqrt(len(dataPAP[0])), unit
                        )
                    )
                    print("Antipoleward:")
                    print("\tN = {0}".format(len(dataPAP[1])))
                    print("\tMedian = {0:.3f} {1}".format(np.median(dataPAP[1]), unit))
                    print("\tMean = {0:.3f} {1}".format(np.mean(dataPAP[1]), unit))
                    print(
                        "\tStandard Dev = {0:.3f} {1}".format(np.std(dataPAP[1]), unit)
                    )
                    print(
                        "\tStandard Error = {0:.3f} {1}".format(
                            np.std(dataPAP[1]) / np.sqrt(len(dataPAP[1])), unit
                        )
                    )
                sys.stdout = original_stdout  # Reset the standard output to its original value            # Display

            print("-" * 30)
            # Get max y value (ceiled to the nearest .01)
            ytop = 1.25 * max([pp.get_height() for pp in ax.patches])
            ymax = math.ceil(ytop * ymax_scaling[xlab]) / (ymax_scaling[xlab])
            ax.set_yticks([0, ymax / 2, ymax])
            ax.set(xlabel="{0} ({1})".format(xlab, unit))
            # Limits and ticks
            ax.set_xlim(left=-1 * xmaxes[xlab], right=xmaxes[xlab])
            ax.set_ylim(bottom=0, top=1.0 * ymax)

            # Plot medians
            for strain_name, col in zip(medians.keys(), cols1):
                meds = medians[strain_name].values()
                for med in meds:
                    ax = plot_median_special(ax, med, 0.9, col)

            # Legend
            if nStrain > 1 and xlab == "Velocity":
                ax.legend(frameon=False, loc="upper left")
                # ax.legend(frameon=False)
            # Set ylabel
            ax.set(ylabel="Probability density")

            # XLABELS
            if xlab == "Lifetime":
                ax.set_xticks([-100, -50, 0, 50, 100])
                ax.set_xticklabels(np.abs(ax.get_xticks()))
            elif xlab == "Velocity":
                if self.config["paperFigure"] == 5:
                    ax.set_xticks([-60, -30, 0, 30, 60])
                else:
                    ax.set_xticks([-100, -50, 0, 50, 100])
            elif xlab == "Average distance from SPB":
                ax.set_xticks([-8000, -4000, 0, 4000, 8000])
                ax.ticklabel_format(
                    style="sci", axis="y", scilimits=(0, 0), useMathText=True
                )
                ax.set_xticklabels(np.abs(ax.get_xticks()))

            elif xlab == "Run displacement":
                ax.set_xticks([-1600, -800, 0, 800, 1600])
                ax.ticklabel_format(
                    style="sci", axis="y", scilimits=(0, 0), useMathText=True
                )

            plt.tight_layout()
            plt.savefig(figname)
            plt.close()

    # }}}

    # GraphStrain_SwitchFrequency2 {{{
    def GraphStrain_SwitchFrequency2(self, figname="graph_switch_frequency.pdf"):
        # Graph comparison bar plot for events per minute

        strains = [strain.label for strain in self.strains]

        # Data
        n_events = np.zeros((len(self.strains), 2))
        dt = 0 * n_events

        for idx, strain in enumerate(self.strains):
            events, times = strain.GetSwitchFrequencyPerMinutePerCell()
            for je, jt in zip(events["P"], times["P"]):
                n_events[idx, 0] += je
                dt[idx, 0] += jt * 60
            for je, jt in zip(events["AP"], times["AP"]):
                n_events[idx, 1] += je
                dt[idx, 1] += jt * 60
        events_per_min = n_events / dt
        events_per_min_err = np.sqrt(n_events) / dt

        df = pd.DataFrame(
            {"Poleward": events_per_min[:, 0], "AntiPoleward": events_per_min[:, 1]},
            index=strains,
        )

        # Plot
        fig, ax = plt.subplots(figsize=(4, 3))
        ax = df.plot(
            kind="bar",
            ax=ax,
            color=["Green", "Red"],
            rot=0,
            yerr=events_per_min_err,
            error_kw=dict(ecolor="k"),
            legend=False,
        )
        ax.set_ylabel("Switching frequency\n(events/sec)")
        ax.set_xlabel("")

        num_cells = [0, 0]
        for idx in range(len(strains)):
            for kymo in self.strains[idx].kymographs:
                if kymo.poles != []:
                    tt = kymo.poles[0].time[-1] - kymo.poles[0].time[0]
                    if tt > 10:
                        num_cells[idx] += 1

        original_stdout = sys.stdout  # Save a reference to the original standard output
        with open("stats.txt", "a") as f:
            sys.stdout = f  # Change the standard output to the file we created.
            print("------------------------------")
            print("\nSwitching Frequency\n")
            for idx, strain in enumerate(strains):
                print("Strain: {0}".format(strain))
                print("   Num Cells: {0}".format(num_cells[idx]))
                print("   Poleward Exit")
                print("      N Events: {0}".format(n_events[idx, 0]))
                print("      Total Time: {0:.3f}".format(dt[idx, 0]))
                print(
                    "      Switching Freq: {0:.5f} sec^-1".format(
                        events_per_min[idx, 0]
                    )
                )
                print(
                    "      Error in switching freq: {0:.5f} sec^-1".format(
                        events_per_min_err[idx, 0]
                    )
                )
                print("   AntiPoleward")
                print("      N Events: {0}".format(n_events[idx, 1]))
                print("      Total Time: {0:.3f}".format(dt[idx, 1]))
                print(
                    "      Switching Freq: {0:.5f} sec^-1".format(
                        events_per_min[idx, 1]
                    )
                )
                print(
                    "      Error in switching freq: {0:.5f} sec^-1".format(
                        events_per_min_err[idx, 1]
                    )
                )

            print("------------------------------")
            sys.stdout = original_stdout  # Reset the standard output to its original value            # Display

        # Set y axis limit and ticks (ceil to nearest 0.02)
        try:
            if self.config["paperFigure"] == 2:
                # ax.set_ylim(top=2.4)
                ax.set_ylim(top=0.024)
            elif self.config["paperFigure"] == 3:
                # ax.set_ylim(top=5.0)
                ax.set_ylim(top=0.06)
            elif self.config["paperFigure"] == 4:
                # ax.set_ylim(top=5.0)
                ax.set_ylim(top=0.06)
            elif self.config["paperFigure"] == 5:
                # ax.set_ylim(top=2.4)
                ax.set_ylim(top=0.024)
            else:
                raise exception("unkown value for paperfigure parameter")
            ymax = ax.get_ylim()[1]
        except:
            ymax = np.max((data[:, :2] + data[:, 2:]).flatten())
            # ymax = math.ceil(ax.get_ylim()[1]*50)/50
            ymax = math.ceil(ymax * 50) / 50
            ax.set_ylim(top=1.5 * ymax)
        ax.set_yticks([0, ymax / 2, ymax])
        ax.set_ylim(bottom=0.0)
        ax.set(xlabel=None)
        # Scientific notation
        ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useMathText=True)

        # Set custom patch colors (Poleward_strain1, Poleward_strain2, AntiP_streain1, AntiP_strain2)
        if len(self.strains) == 1:
            c1 = self.strains[0].color
            # c1 = [68, 111, 200]
            cols = [c1, c1]
            # cols = [tuple(np.array(x)/255) for x in cols]
            for idx, (pp, col) in enumerate(zip(ax.patches, cols)):
                pp.set_facecolor(col)
                pp.set_alpha(0.6)
                pp.set_edgecolor("white")
                if idx < len(strains):
                    pp.set_hatch("////")

            cols = [c1, c1]
            cols = [tuple(np.array(x) / 255) for x in cols]
            labels = ["Poleward", "Antipoleward"]
            hatching = ["////", ""]
            handles = [
                matplotlib.patches.Rectangle(
                    (0, 0),
                    1,
                    1,
                    facecolor=cols[idx],
                    alpha=0.6,
                    label=labels[idx],
                    hatch=hatching[idx],
                    edgecolor="white",
                )
                for idx in range(len(labels))
            ]

        elif len(self.strains) == 2:

            c1 = self.strains[0].color
            c2 = self.strains[1].color
            # c1 = [68, 111, 200]
            # c2 = [220, 95, 60]
            cols = [c1, c2, c1, c2]
            # cols = [tuple(np.array(x)/255) for x in cols]
            for idx, (pp, col) in enumerate(zip(ax.patches, cols)):
                pp.set_facecolor(col)
                pp.set_alpha(0.6)
                pp.set_edgecolor("white")
                if idx < len(strains):
                    pp.set_hatch("////")

            cols = [c1, c1, c2, c2]
            # cols = [tuple(np.array(x)/255) for x in cols]
            labels = ["Poleward", "Antipoleward", "Poleward", "Antipoleward"]
            hatching = ["////", "", "////", ""]
            handles = [
                matplotlib.patches.Rectangle(
                    (0, 0),
                    1,
                    1,
                    facecolor=cols[idx],
                    alpha=0.6,
                    label=labels[idx],
                    hatch=hatching[idx],
                    edgecolor="white",
                )
                for idx in range(len(labels))
            ]

        else:
            raise Exception("only coded for 1 or 2 strains")

        # ax.legend(handles, labels, loc='upper left', frameon=False)
        ax.legend("", frameon=False)
        plt.tight_layout()
        fig.savefig(figname)
        plt.close()

    # }}}

    # GraphStrain_EventsPerMinute2 {{{
    def GraphStrain_EventsPerMinute2(self, figname="graph_events_per_second.pdf"):
        # Graph comparison bar plot for events per minute

        strains = [strain.label for strain in self.strains]

        # Data
        n_events = np.zeros((len(self.strains), 2))
        dt = 0 * n_events

        for idx, strain in enumerate(self.strains):
            events, times = strain.GetDirectionalEventsPerMinutePerCell()
            for je, jt in zip(events["P"], times["P"]):
                n_events[idx, 0] += je
                # convert times to seconds
                dt[idx, 0] += jt * 60
            for je, jt in zip(events["AP"], times["AP"]):
                n_events[idx, 1] += je
                dt[idx, 1] += jt * 60
        events_per_min = n_events / dt
        events_per_min_err = np.sqrt(n_events) / dt

        df = pd.DataFrame(
            {"Poleward": events_per_min[:, 0], "AntiPoleward": events_per_min[:, 1]},
            index=strains,
        )

        # Plot
        fig, ax = plt.subplots(figsize=(4, 3))
        ax = df.plot(
            kind="bar",
            ax=ax,
            color=["Green", "Red"],
            rot=0,
            yerr=events_per_min_err,
            error_kw=dict(ecolor="k"),
            legend=False,
        )
        ax.set_ylabel("Directional events \n per second")
        ax.set_xlabel("")

        num_cells = [0, 0]
        for idx in range(len(strains)):
            for kymo in self.strains[idx].kymographs:
                if kymo.poles != []:
                    tt = kymo.poles[0].time[-1] - kymo.poles[0].time[0]
                    if tt > 10:
                        num_cells[idx] += 1

        original_stdout = sys.stdout  # Save a reference to the original standard output
        with open("stats.txt", "a") as f:
            sys.stdout = f  # Change the standard output to the file we created.
            print("------------------------------")
            print("\nEvents per second\n")
            for idx, strain in enumerate(strains):
                print("Strain: {0}".format(strain))
                print("   Num Cells: {0}".format(num_cells[idx]))
                print("   Poleward")
                print("      N Events: {0}".format(n_events[idx, 0]))
                print("      Total Time: {0:.3f} sec".format(dt[idx, 0]))
                print(
                    "      Events per sec: {0:.5f} sec^-1".format(
                        events_per_min[idx, 0]
                    )
                )
                print(
                    "      Error in events per sec: {0:.5f} sec^-1".format(
                        events_per_min_err[idx, 0]
                    )
                )
                print("   AntiPoleward")
                print("      N Events: {0}".format(n_events[idx, 1]))
                print("      Total Time: {0:.3f} sec".format(dt[idx, 1]))
                print(
                    "      Events per sec: {0:.5f} sec^-1".format(
                        events_per_min[idx, 1]
                    )
                )
                print(
                    "      Error in events per sec: {0:.5f} sec^-1".format(
                        events_per_min_err[idx, 1]
                    )
                )
            print("------------------------------")
            sys.stdout = original_stdout  # Reset the standard output to its original value            # Display

        # Set y axis limit and ticks (ceil to nearest 0.02)
        try:
            if self.config["paperFigure"] == 2:
                # ax.set_ylim(top=0.8)
                ax.set_ylim(top=0.014)
            elif self.config["paperFigure"] == 3:
                # ax.set_ylim(top=1.8)
                ax.set_ylim(top=0.03)
            elif self.config["paperFigure"] == 4:
                # ax.set_ylim(top=1.8)
                ax.set_ylim(top=0.03)
            elif self.config["paperFigure"] == 5:
                # ax.set_ylim(top=0.8)
                ax.set_ylim(top=0.014)
            else:
                raise exception("unkown value for paperfigure parameter")
            ymax = ax.get_ylim()[1]
        except:
            ymax = np.max((data[:, :2] + data[:, 2:]).flatten())
            # ymax = math.ceil(ax.get_ylim()[1]*50)/50
            ymax = math.ceil(ymax * 50) / 50
            ax.set_ylim(top=1.5 * ymax)

        ax.set_yticks([0, ymax / 2, ymax])
        ax.set_ylim(bottom=0.0)
        ax.set(xlabel=None)
        # Scientific notation
        ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useMathText=True)

        if len(self.strains) == 1:
            c1 = self.strains[0].color
            # c1 = [68, 111, 200]
            cols = [c1, c1]
            # cols = [tuple(np.array(x)/255) for x in cols]
            for idx, (pp, col) in enumerate(zip(ax.patches, cols)):
                pp.set_facecolor(col)
                pp.set_alpha(0.6)
                pp.set_edgecolor("white")
                if idx < len(strains):
                    pp.set_hatch("////")

            cols = [c1, c1]
            cols = [tuple(np.array(x) / 255) for x in cols]
            labels = ["Poleward", "Antipoleward"]
            hatching = ["////", ""]
            handles = [
                matplotlib.patches.Rectangle(
                    (0, 0),
                    1,
                    1,
                    facecolor=cols[idx],
                    alpha=0.6,
                    label=labels[idx],
                    hatch=hatching[idx],
                    edgecolor="white",
                )
                for idx in range(len(labels))
            ]

        elif len(self.strains) == 2:

            c1 = self.strains[0].color
            c2 = self.strains[1].color
            # c1 = [68, 111, 200]
            # c2 = [220, 95, 60]
            cols = [c1, c2, c1, c2]
            # cols = [tuple(np.array(x)/255) for x in cols]
            for idx, (pp, col) in enumerate(zip(ax.patches, cols)):
                pp.set_facecolor(col)
                pp.set_alpha(0.6)
                pp.set_edgecolor("white")
                if idx < len(strains):
                    pp.set_hatch("////")

            cols = [c1, c1, c2, c2]
            # cols = [tuple(np.array(x)/255) for x in cols]
            labels = ["Poleward", "Antipoleward", "Poleward", "Antipoleward"]
            hatching = ["////", "", "////", ""]
            handles = [
                matplotlib.patches.Rectangle(
                    (0, 0),
                    1,
                    1,
                    facecolor=cols[idx],
                    alpha=0.6,
                    label=labels[idx],
                    hatch=hatching[idx],
                    edgecolor="white",
                )
                for idx in range(len(labels))
            ]

        else:
            raise Exception("only coded for 1 or 2 strains")

        ax.legend(handles, labels, loc="upper left", frameon=False)
        plt.tight_layout()
        fig.savefig(figname)
        plt.close()

    # }}}

    # GraphStrain_FractionMovement {{{
    def GraphStrain_FractionMovement(self, figname="graph_fraction_kymo_movement.pdf"):
        # Graph comparison bar plot for events per minute

        fracMove = [
            strain.GetFractionKymographsWithMovement() for strain in self.strains
        ]
        n_total = [len(strain.kymographs) for strain in self.strains]
        n_move = [int(jp * np) for jp, np in zip(fracMove, n_total)]
        strains = [strain.label for strain in self.strains]

        # Colors
        cols1 = [cstrain.color for cstrain in self.strains]

        num_cells = [0, 0]
        for idx in range(len(strains)):
            for kymo in self.strains[idx].kymographs:
                if kymo.poles != []:
                    tt = kymo.poles[0].time[-1] - kymo.poles[0].time[0]
                    if tt > 10:
                        num_cells[idx] += 1

        original_stdout = sys.stdout  # Save a reference to the original standard output
        with open("stats.txt", "a") as f:
            sys.stdout = f  # Change the standard output to the file we created.
            print("------------------------------")
            print("\nFraction kymograph movement\n\n")
            for idx, strain in enumerate(strains):
                print("  Strain: {0}".format(strain))
                print("    Percentage: {0:.3f}".format(fracMove[idx]))
                print("    N: {0}\n".format(num_cells[idx]))
            print("------------------------------")
            sys.stdout = original_stdout  # Reset the standard output to its original value            # Display

        #  Bar plot
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.bar(strains, fracMove, color=cols1, width=0.5, alpha=0.6)
        # ax.set_xlabel("Strain")
        ax.set_ylabel("Fraction of cells\nwith movement")
        ax.set_ylim(top=1.0)
        ax.set_yticks([0.0, 0.5, 1.0])
        ax.set_xlim(left=-0.75, right=len(strains) - 1 + 0.75)

        handles = [
            plt.Rectangle((0, 0), 1, 1, color=cols1[idx], alpha=0.6)
            for idx in range(len(strains))
        ]
        # plt.legend(handles, strains, loc='upper left', frameon=False)
        plt.tight_layout()
        fig.savefig(figname)
        plt.close()

    # }}}

    # GraphStrainMedianValues{{{
    def GraphStrainMedianValues(self, figname="graph_median_lifetime.pdf"):
        # Graph comparison bar plot for median lifetime

        graphhistvars = (
            [
                "GetRunLengths",
                "Run displacement",
                r"$\mu$" + "m",
                "strain_median_runlength.pdf",
            ],
            [
                "GetVelocities_nm_per_sec",
                "Velocity",
                "nm/s",
                "strain_median_velocity.pdf",
            ],
            ["GetLifetimes_min", "Lifetime", r"min", "strain_median_lifetime.pdf"],
            [
                "GetAverageDistances",
                "Average distance from SPB",
                r"$\mu$" + "m",
                "strain_median_avg_pos.pdf",
            ],
        )

        for fcn, ylab, unit, figname in graphhistvars:
            # Data
            # Row: strains
            # Col: Poleward Mean, AntiPoleward Mean, Poleward STD, Antipoleward STD
            data = np.zeros((len(self.strains), 4))
            count = np.zeros(len(self.strains))
            for idx, strain in enumerate(self.strains):
                funcData = getattr(strain, fcn)
                events = funcData()
                count[idx] = len(strain.kymographs)
                data[idx, 0] = np.mean(events[0])
                data[idx, 2] = np.std(events[0]) / np.sqrt(count[idx])
                data[idx, 1] = np.mean(events[1])
                data[idx, 3] = np.std(events[1]) / np.sqrt(count[idx])

            strains = [strain.label for strain in self.strains]

            # Create pd Dataframe for plotting
            df = pd.DataFrame(
                data,
                columns=["Poleward", "Antipoleward", "std_P", "std_AP"],
                index=strains,
            )

            # Plot
            fig, ax = plt.subplots(figsize=(4, 3))
            # convert the std columns to an array
            yerr = df[["std_P", "std_AP"]].to_numpy().T
            ax = df[["Poleward", "Antipoleward"]].plot(
                kind="bar",
                ax=ax,
                color=["Green", "Red"],
                rot=0,
                # yerr=yerr, error_kw=dict(ecolor='k'),legend=False, xlabel=None)
                legend=False,
                xlabel=None,
            )
            # ax.set_xlabel("Strain")
            ax.set_ylabel("Median\n{0}\n({1})".format(ylab, unit))

            # Set y axis limit and ticks (ceil to nearest 0.02)
            if ylab == "Velocity":  # nearest 4
                ymax = np.max((data[:, :2] + data[:, 2:]).flatten())
                ymax = math.ceil(ymax / 4) * 4
            else:
                ymax = np.max((data[:, :2] + data[:, 2:]).flatten())
                # ymax = math.ceil(ax.get_ylim()[1]*50)/50
                ymax = math.ceil(ymax * 50) / 50
            ax.set_ylim(top=1.4 * ymax)
            ax.set_yticks([0, ymax / 2, ymax])
            # for jj in range(2):
            # ax.text(jj, ymax, 'N cells = {0}'.format(count[jj]),
            # ha='center', color='black', fontsize=8)

            # Set custom patch colors (Poleward_strain1, Poleward_strain2, AntiP_streain1, AntiP_strain2)
            if len(self.strains) == 1:
                c1 = self.strains[0].color
                # c1 = [68, 111, 200]
                cols = [c1, c1]
                # cols = [tuple(np.array(x)/255) for x in cols]
                for idx, (pp, col) in enumerate(zip(ax.patches, cols)):
                    pp.set_facecolor(col)
                    pp.set_alpha(0.6)
                    pp.set_edgecolor("white")
                    if idx < len(strains):
                        pp.set_hatch("////")

                cols = [c1, c1]
                # cols = [tuple(np.array(x)/255) for x in cols]
                labels = ["Poleward", "Antipoleward"]
                hatching = ["////", ""]
                handles = [
                    matplotlib.patches.Rectangle(
                        (0, 0),
                        1,
                        1,
                        facecolor=cols[idx],
                        alpha=0.6,
                        label=labels[idx],
                        hatch=hatching[idx],
                        edgecolor="white",
                    )
                    for idx in range(len(labels))
                ]

            elif len(self.strains) == 2:

                c1 = self.strains[0].color
                c2 = self.strains[1].color
                # c1 = [68, 111, 200]
                # c2 = [220, 95, 60]
                cols = [c1, c2, c1, c2]
                # cols = [tuple(np.array(x)/255) for x in cols]
                for idx, (pp, col) in enumerate(zip(ax.patches, cols)):
                    pp.set_facecolor(col)
                    pp.set_alpha(0.6)
                    pp.set_edgecolor("white")
                    if idx < len(strains):
                        pp.set_hatch("////")

                cols = [c1, c1, c2, c2]
                # cols = [tuple(np.array(x)/255) for x in cols]
                labels = [
                    "Poleward, {0}".format(strains[0]),
                    "Antipoleward, {0}".format(strains[0]),
                    "Poleward, {0}".format(strains[1]),
                    "Antipoleward, {0}".format(strains[1]),
                ]
                hatching = ["////", "", "////", ""]
                handles = [
                    matplotlib.patches.Rectangle(
                        (0, 0),
                        1,
                        1,
                        facecolor=cols[idx],
                        alpha=0.6,
                        label=labels[idx],
                        hatch=hatching[idx],
                        edgecolor="white",
                    )
                    for idx in range(len(labels))
                ]

            else:
                raise Exception("only coded for 1 or 2 strains")

            ax.legend(handles, labels, loc="upper left", frameon=False)
            plt.tight_layout()
            fig.savefig(figname)
            plt.close()

    # }}}

    # GraphStrain_AvgStartEnd {{{
    def GraphStrain_AvgStartEnd(self, figname="graph_fraction_kymo_movement.pdf"):

        freqP = [strain.GetFractionKymographsWithMovement() for strain in self.strains]

        strains = [strain.label for strain in self.strains]
        # Create pd Dataframe for plotting
        seriesP = pd.Series(freqP, index=strains)

        # Plot
        fig, ax = plt.subplots(figsize=(6, 4))

        df = pd.DataFrame({"Movements": seriesP})
        df.plot.bar(ax=ax, color=["RebeccaPurple"], rot=0)
        ax.set_xlabel("Strain")
        ax.set_ylabel("Fraction of Cells\nwith Movement")

        plt.tight_layout()
        fig.savefig(figname)
        plt.close()

    # }}}

    # GraphStrainScatter {{{
    def GraphStrainScatter(
        self,
        strains,
        x,
        y,
        xlab=None,
        ylab=None,
        xunit="",
        yunit="",
        figname="scatter.pdf",
    ):

        # Special x limits
        xmaxes = {
            "Run displacement": 1600,
            "Velocity": 60.0,
            "Lifetime": 80,
            "Intensity": 1000,
        }

        # 2 axes. Poleward and antipoleward
        fig, axs = plt.subplots(1, 2, figsize=(6, 3), sharey=True)
        cols = sns.color_palette("husl", len(strains))
        directions = ["Poleward", "Antipoleward"]

        for strain, c in zip(strains, cols):
            if x == "Intensity":
                xx = strain.GetIntensities()
            elif x == "Run displacement":
                xx = strain.GetRunLengths()
            elif x == "Velocity":
                xx = strain.GetVelocities()
            elif x == "Lifetime":
                xx = strain.GetLifetimes()
            # elif x == "Average distance from SPB":
                # xx = strain.GetAverageDistances()
            if y == "Intensity":
                yy = strain.GetIntensities()
            elif y == "Run displacement":
                yy = strain.GetRunLengths()
            elif y == "Velocity":
                yy = strain.GetVelocities()
            elif y == "Lifetime":
                yy = strain.GetLifetimes()
            # elif y == "Average distance from SPB":
                # yy = strain.GetAverageDistances()

            for idx, ax in enumerate(axs):
                ax.scatter(
                    xx[idx],
                    yy[idx],
                    s=12,
                    alpha=0.8,
                    color=c,
                    edgecolors="none",
                    label=strain.label,
                )
                ax.set_title(directions[idx])
                # ax.grid(True)

        if xlab is not None:
            axs[0].set_xlabel("{0} ({1})".format(xlab, xunit))
            axs[1].set_xlabel("{0} ({1})".format(xlab, xunit))
        if ylab is not None:
            axs[0].set_ylabel("{0} ({1})".format(ylab, yunit))
        for ax in axs:
            ax.legend()
            ax.set_xlim(left=0, right=xmaxes[x])
            ax.set_ylim(bottom=0, top=xmaxes[y])

        plt.tight_layout()
        plt.savefig("scatter_{0}_{1}.pdf".format(x, y))
        plt.close()

    # }}}

    # GraphSpindleIntensity {{{
    def GraphSpindleIntensity(self, strains, lrange=[1, 10], gname=None):
        # Graph spindle intensity between poles

        intStrain = np.zeros((len(strains), 100))
        cols = sns.color_palette("husl", len(strains))
        xran = np.linspace(0, 1, 100)
        for k, strain in enumerate(strains):

            fig, ax = plt.subplots(figsize=(9, 6))
            # Find spindle intensities for all kymographs
            intensities = None
            for i, kymo in enumerate(strain.kymographs):
                intense = kymo.FindIntensityAlongSpindle(lrange=lrange)
                if intense is not None:
                    if intensities is None:
                        intensities = np.mean(intense, axis=0)
                    else:
                        intensities = np.vstack((intensities, np.mean(intense, axis=0)))

            try:
                intStrain[k, :] = np.mean(intensities, axis=0)
            except:
                pdb.set_trace()
                print("1")
            # Plot
            for row in intensities:
                ax.plot(xran, row, color="blue")
            ax.plot(xran, np.mean(intensities, axis=0), color="red", linewidth=4)
            ax.set_ylabel("Intensity (AU)")
            ax.set_xlabel("Position along spindle (normalized)")
            ax.set_title("Cut7 intensity - {0}".format(strain.label))
            fig.savefig("spindle_intensity_{0}.pdf".format(strain.label))
            plt.close()

        # Make a comparison figure
        fig, ax = plt.subplots(figsize=(9, 6))
        for strn, row, c in zip(strains, intStrain, cols):
            ax.plot(xran, row, color=c, linewidth=4, label=strn.label)
        ax.set_ylabel("Intensity (AU)")
        ax.set_xlabel("Position along spindle (normalized)")
        ax.set_title("Cut7 intensity")
        ax.legend()
        figname = "spindle_intensity_all.pdf"
        if gname is not None:
            figname = figname[:-4] + "_{0}.pdf".format(gname)
            fig.suptitle(gname)
        fig.savefig(figname)
        plt.close()

    # }}}

    # PlotTracksByState {{{
    def PlotTracksByState(self, k=5):
        # Plot individual curved tracks with poles

        # Plot all tracks overlayed without poles
        cols = {
            "Inactive": "blue",
            "Poleward": "green",
            "Antipoleward": "red",
        }

        for strain in self.strains:
            strain.PlotTrackByStates(cols, k=k)

    # }}}

    # PlotAllTracks {{{
    def PlotAllTracks(self):
        # Plot all tracks

        cols = {
            "Inactive": "blue",
            "Poleward": "green",
            "Antipoleward": "red",
        }

        for strain in self.strains:
            strain.PlotAllTracks(cols)

    # }}}

    # DisplayTracksStatistics {{{
    def DisplayTracksStatistics(self):
        # Display statistics about the tracks

        print("------------------------------------------")
        print("------------------------------------------")
        print("------------ Track Statistics ------------")
        print("------------------------------------------\n")
        print("Number of tracks:")
        for strain in self.strains:
            print("  {0} : {1}\n".format(strain.label, len(strain.tracks)))

        print("------------------------------------------")
        print("------------------------------------------")

    # }}}


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values - average) ** 2, weights=weights)
    std = np.sqrt(variance)
    serr = std / np.sqrt(len(values))
    return average, std, serr


#########################################################
if __name__ == "__main__":
    x = KymographAnalysis()
