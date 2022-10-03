#!/usr/bin/env python

import os, pdb
from .Load import Load
from .Kymograph import Kymograph
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import math, random
import pickle
from pathlib import Path
from scipy import interpolate, signal


class Strain:
    def __init__(self, trackpaths, label="xxx"):

        self.paths = trackpaths
        self.label = label
        self.LoadKymographs()
        self.tracks = []

    def LoadKymographs(self):
        # Initialize kymograph classes for each loaded file
        self.kymographs = []
        for pth in self.paths:
            print(pth)
            kname = pth.split("/")[-1]
            self.kymographs += [Kymograph(fname=pth)]

    def GetTracks(self, spindle_length=None):
        # Get tracks that lie within the spindle lengths defined. Trims the tracks
        # Combine tracks from all kymographs
        self.tracks = []
        for kymo in self.kymographs:
            for track in kymo.tracks:
                trimmed = track.Trim(lrange=spindle_length)
                if trimmed is not None:
                    self.tracks += [trimmed]

    def TrimUsingKmeansLabel(self, kmean_label):

        # load kmeans model
        kmeans_path = Path(self.paths[0]).parent / "kmeans.pickle"
        with open(kmeans_path, "rb") as f:
            model = pickle.load(f)

        for kymo in self.kymographs:
            # Only do stuff if its bipolar
            if len(kymo.poles) == 1:
                continue

            # Times
            time = np.array(
                sorted(np.hstack((kymo.poles[0].time, kymo.poles[1].time)))[1::10]
            )
            time = np.linspace(time[0], time[-1], int(np.ceil(time[-1] - time[0])))

            # Calculate spindle length, velocity, acceleration
            clen = np.absolute(kymo.poles[1].ifunc(time) - kymo.poles[0].ifunc(time))
            cvel = list((clen[1:] - clen[:-1]) / (time[1:] - time[:-1]))
            cvel.insert(0, cvel[0])
            cvel = np.array(cvel).reshape(-1, 1)
            # use velocity to predict label using fitted model
            labels_raw = model.predict(cvel)
            labels = self.ForceLabelsOneWay(
                self.SmoothClassifiedLabels(labels_raw, span=100)
            )
            if np.max(labels) == 0 or np.max(clen) < 2:
                AB_transition = -1
                if kmean_label == 0:
                    time_keep = [time[0], time[-1]]
                elif kmean_label == 1:
                    time_keep = [-1, -1]
            elif np.min(labels) == 1:
                AB_transition = -1
                if kmean_label == 0:
                    time_keep = [-1, -1]
                elif kmean_label == 1:
                    time_keep = [time[0], time[-1]]
            else:
                AB_transition = time[np.where((labels == 1) & (clen > 2))[0][0]]
                if kmean_label == 0:
                    time_keep = [time[0], AB_transition]
                elif kmean_label == 1:
                    time_keep = [AB_transition, time[-1]]
            # if kymo.label == '/Users/saadjansari/Documents/Projects/ImageAnalysis/KymoAnalysis/data/bipolar/wild type/MAX_1032_100msR_50msG_7Z_004_cell A_KYMOGRAPH':
            # pdb.set_trace()
            # print(AB_transition)
            # print(time_keep)
            print("Total time = {0:.2f} - {1:.2f}".format(time[0], time[-1]))
            print("Anaphase B = {0:.2f}".format(AB_transition))
            print("Kmeans Label = {0}".format(kmean_label))
            print("Time 2 keep = {0}".format(time_keep))
            kymo.TrimBasedOnTime(time_keep)

    # SmoothClassifiedLabels {{{
    def SmoothClassifiedLabels(self, label, span=100):

        # smooth_data {{{
        def smooth_data(arr, span):
            re = np.convolve(arr, np.ones(span * 2 + 1) / (span * 2 + 1), mode="same")

            # The "my_average" part: shrinks the averaging window on the side that
            # reaches beyond the data, keeps the other side the same size as given
            # by "span"
            re[0] = np.average(arr[:span])
            for i in range(1, span + 1):
                re[i] = np.average(arr[: i + span])
                re[-i] = np.average(arr[-i - span :])
            return re

        # }}}

        # Smoothed Labels
        label_new = np.where(
            np.array(smooth_data(label, min([span, int(len(label) / 2)]))) >= 0.5, 1, 0
        )

        # Once 1, always 1
        # label_perm = [max(label_new[:1+jj]) for jj in range(len(label_new))]
        return label_new

    # }}}

    # ForceLabelsOneWay {{{
    def ForceLabelsOneWay(self, label):
        labels = [np.max(label[: 1 + idx]) for idx in range(len(label))]
        return np.array(labels)

    # }}}

    def TossFarTracks(self, threshold):
        # Toss tracks that start above a threshold distance from the first pole
        self.tracks = []
        for kymo in self.kymographs:
            tracksKeep = []
            for track in kymo.tracks:
                if track.CalcPositionRelative()[0, 0] < threshold:
                    tracksKeep.append(track)
            kymo.tracks = tracksKeep
        self.GetTracks()

    def TossCloseTracks(self, threshold):
        # Toss tracks that start above a threshold distance from the first pole
        self.tracks = []
        for kymo in self.kymographs:
            tracksKeep = []
            for track in kymo.tracks:
                if track.CalcPositionRelative()[0, 0] > threshold:
                    tracksKeep.append(track)
            kymo.tracks = tracksKeep
        self.GetTracks()

    def GetSegmentsPAP(self):
        # Splits the tracks into segments and returns [poleward, antipoleward]
        segs = {"Poleward": [], "Antipoleward": []}
        bad_cnt = 0
        good_cnt = 0
        for track in self.tracks:
            segments, _ = track.SplitTrack()

            # Toss short-time segments
            for seg in segments:
                if seg.time[-1] - seg.time[0] < 2 * seg.time_step:
                    bad_cnt += 1
                elif seg.direction is "Poleward":
                    good_cnt += 1
                    segs["Poleward"] += [seg]
                elif seg.direction is "Antipoleward":
                    good_cnt += 1
                    segs["Antipoleward"] += [seg]

        # if good_cnt + bad_cnt > 0:
        # print('Major segments : {0} ({1:.2f}%)'.format(good_cnt, 100*good_cnt/(good_cnt+bad_cnt)))
        return [segs["Poleward"], segs["Antipoleward"]]

    def FilterSegments(self, segments):
        # Filter segments by imposing restrictions on velocity, run length and lifetimes
        # Velocities
        print("1")

    def GetRunLengths(self):
        # Get run lengths of poleward and antipoleward tracks (units nm)
        segsPAP = self.GetSegmentsPAP()
        runlens = []
        for segs in segsPAP:
            runlen = [
                1000
                * np.absolute(
                    seg.CalcPositionRelative()[0, -1] - seg.CalcPositionRelative()[0, 0]
                )
                for seg in segs
            ]
            runlens += [runlen]
        return runlens

    def GetVelocities(self):
        # Get velocities of poleward and antipoleward tracks (units nm/sec)
        segsPAP = self.GetSegmentsPAP()
        vels = []
        for segs in segsPAP:
            vel = [1000 * seg.CalcVelocity()[0] for seg in segs]
            vels += [vel]
        return vels

    def GetLifetimes(self):
        # Get velocities of poleward and antipoleward tracks (units sec)
        segsPAP = self.GetSegmentsPAP()
        lifes = []
        for segs in segsPAP:
            life = [seg.time[-1] - seg.time[0] for seg in segs]
            lifes += [life]
        return lifes

    def GetIntensities(self):
        # Get velocities of poleward and antipoleward tracks
        segsPAP = self.GetSegmentsPAP()
        ins = []
        for segs in segsPAP:
            inss = [seg.CalcIntensity() for seg in segs]
            ins += [inss]
        return ins

    def GetTotalSwitches(self):
        # Get total switches out of state A to state B (states: poleward,antipoleward,inactive)

        switches = {"P": 0, "AP": 0, "I": 0}
        labs = ["P", "AP", "I"]
        for track in self.tracks:
            segments, trans = track.SplitTrack()
            for lab in labs:
                switches[lab] += sum([a for k, a in trans[lab].items()])
        return switches

    def GetFractionKymographsWithMovement(self):
        # Get fraction of kymographs with movement events

        nMovement = 0
        for kymo in self.kymographs:
            nAdd = 0
            for track in kymo.tracks:
                trks_all, _ = track.SplitTrack()
                for mini in trks_all:
                    if mini.direction != "Inactive":
                        nAdd = 1
            nMovement += nAdd
        return nMovement / len(self.kymographs)

    def GetDirectionalEventsPerMinute(self):
        # Get total number of directional events per minute

        events = {"P": 0, "AP": 0}
        events_per_min = {"P": 0, "AP": 0}
        for track in self.tracks:
            segs, _ = track.SplitTrack()
            for seg in segs:
                if seg.direction == "Poleward":
                    events["P"] += 1
                elif seg.direction == "Antipoleward":
                    events["AP"] += 1

        # total kymograph time
        time_total = 0
        for kymo in self.kymographs:
            time_total += kymo.poles[0].time[-1] - kymo.poles[0].time[0]
        events_per_min["P"] = events["P"] / (time_total / 60)
        events_per_min["AP"] = events["AP"] / (time_total / 60)
        return events_per_min, events

    def GetDirectionalEventsPerMinutePerCell(self):
        # Get total number of directional events per minute per cell

        events = {"P": [], "AP": []}
        times = {"P": [], "AP": []}
        for kymo in self.kymographs:
            nP = 0
            nAP = 0
            for track in kymo.tracks:
                segs, _ = track.SplitTrack()
                for seg in segs:
                    if seg.direction == "Poleward":
                        nP += 1
                    elif seg.direction == "Antipoleward":
                        nAP += 1
            if nP + nAP > 0:
                time_total = kymo.poles[0].time[-1] - kymo.poles[0].time[0]
            else:
                time_total = 1
            if time_total > 1:
                events["P"].append(nP)
                events["AP"].append(nAP)
                times["P"].append(time_total / 60)
                times["AP"].append(time_total / 60)
        # events['P'].append( nP/(time_total/60))
        # events['AP'].append( nAP/(time_total/60))
        return events, times

    def GetTotalDirectionalTime(self):
        # Find total number of directed time

        times = {"P": 0.01, "AP": 0.01, "I": 0.01}
        for track in self.tracks:
            # Split the track into linear tracks
            segs, _ = track.SplitTrack()
            # Calculate lifetimes and sum it all up for each direction of the track
            for seg in segs:
                if seg.direction == "Poleward":
                    times["P"] += seg.time[-1] - seg.time[0]
                elif seg.direction == "Antipoleward":
                    times["AP"] += seg.time[-1] - seg.time[0]
                elif seg.direction == "Inactive":
                    times["I"] += seg.time[-1] - seg.time[0]
                else:
                    raise ValueError("what is this unknown line direction")
        return times

    def GetTotalDirectionalTimeMinutes(self):
        # Find total number of directed time

        times = {"P": 0.01, "AP": 0.01, "I": 0.01}
        for track in self.tracks:
            # Split the track into linear tracks
            segs, _ = track.SplitTrack()
            # Calculate lifetimes and sum it all up for each direction of the track
            for seg in segs:
                if seg.direction == "Poleward":
                    times["P"] += (seg.time[-1] - seg.time[0]) / 60
                elif seg.direction == "Antipoleward":
                    times["AP"] += (seg.time[-1] - seg.time[0]) / 60
                elif seg.direction == "Inactive":
                    times["I"] += (seg.time[-1] - seg.time[0]) / 60
                else:
                    raise ValueError("what is this unknown line direction")
        return times

    def GetSwitchFrequencyPerMinutePerCell(self):

        events = {"P": [], "AP": [], "I": []}
        times_all = {"P": [], "AP": [], "I": []}
        for kymo in self.kymographs:

            # Get total track time
            times = {"P": 10 ** -7, "AP": 10 ** -7, "I": 10 ** -7}
            for track in kymo.tracks:
                # Split the track into linear tracks
                segs, _ = track.SplitTrack()
                # Calculate lifetimes and sum it all up for each direction of the track
                for seg in segs:
                    if seg.direction == "Poleward":
                        times["P"] += (seg.time[-1] - seg.time[0]) / 60
                    elif seg.direction == "Antipoleward":
                        times["AP"] += (seg.time[-1] - seg.time[0]) / 60
                    elif seg.direction == "Inactive":
                        times["I"] += (seg.time[-1] - seg.time[0]) / 60
                    else:
                        raise ValueError("what is this unknown line direction")

            # Get total switches
            switches = {"P": 0, "AP": 0, "I": 0}
            labs = ["P", "AP", "I"]
            for track in kymo.tracks:
                segments, trans = track.SplitTrack()
                for lab in labs:
                    switches[lab] += sum([a for k, a in trans[lab].items()])

            # switch frequencies
            for lab in labs:
                events[lab].append(switches[lab])
                times_all[lab].append(times[lab])
        return events, times_all

    def GetStartDistances(self):
        # Get start distances

        dist_P = []
        dist_AP = []
        for track in self.tracks:
            segs, _ = track.SplitTrack()
            for seg in segs:
                if seg.direction == "Poleward":
                    dist_P.append(seg.CalcPositionRelative()[0, 0])
                elif seg.direction == "Antipoleward":
                    dist_AP.append(seg.CalcPositionRelative()[0, 0])

        return [dist_P, dist_AP]

    def GetEndDistances(self):
        # Get end distances

        dist_P = []
        dist_AP = []
        for track in self.tracks:
            segs, _ = track.SplitTrack()
            for seg in segs:
                if seg.direction == "Poleward":
                    dist_P.append(seg.CalcPositionRelative()[0, -1])
                elif seg.direction == "Antipoleward":
                    dist_AP.append(seg.CalcPositionRelative()[0, -1])

        return [dist_P, dist_AP]

    def GetAverageDistances(self):
        # Get average distances

        segsPAP = self.GetSegmentsPAP()
        avgdists = []
        for segs in segsPAP:
            avgdist = []
            for seg in segs:
                ifunc = interpolate.interp1d(
                    seg.time, seg.CalcPositionRelative()[0, :], kind="linear"
                )
                avgdist.extend(
                    [
                        1000 * dd
                        for dd in ifunc(
                            np.arange(seg.time[0], seg.time[-1], seg.time_step)
                        )
                    ]
                )
            # avgdist = [np.mean(seg.CalcPositionRelative()[0, :]) for seg in segs]
            avgdists += [avgdist]
        return avgdists

    def GraphPAP_RunLengths(self, axs, **kwargs):
        lens_pap = self.GetRunLengths()
        # Toss runlengths over 2 micron
        for idx, lens in enumerate(lens_pap):
            ld = [i for i in lens if i < 2]
            lens_pap[idx] = ld
        self.GraphPAP(lens_pap, axs, unit=r"$\mu$" + "m", **kwargs)

    def GraphPAP_Velocities(self, axs, **kwargs):
        vels_pap = self.GetVelocities()
        # Convert vel from micron/sec to nm/sec. Toss vel over 200 micron/sec
        for idx, v in enumerate(vels_pap):
            ld = [i * 1000 for i in v if i < 0.2]
            vels_pap[idx] = ld
        self.GraphPAP(vels_pap, axs, unit="nm/s", **kwargs)

    def GraphPAP_Lifetimes(self, axs, **kwargs):
        lifes_pap = self.GetLifetimes()
        # toss lifetimes over 100 sec
        for idx, life in enumerate(lifes_pap):
            ld = [i for i in life if i < 100]
            lifes_pap[idx] = ld
        self.GraphPAP(lifes_pap, axs, unit="s", **kwargs)

    def GraphPAP_StartPosition(self, axs, **kwargs):
        startPos = self.GetStartDistances()
        self.GraphPAP(startPos, axs, unit=r"$\mu$" + "m", **kwargs)

    def GraphPAP_EndPosition(self, axs, **kwargs):
        endPos = self.GetEndDistances()
        self.GraphPAP(endPos, axs, unit=r"$\mu$" + "m", **kwargs)

    def GraphPAP(
        self, dat, axs, col="m", lab=None, unit="", xmax=None, xlab=None, ylab=None
    ):

        # pdb.set_trace()
        for datt, ax in zip(dat, axs):

            # Get x axis max
            if xmax is None:
                # xmax = math.ceil( max(datt))
                xmax = max(datt)

            # find bin edges
            nbins = 16
            bins = np.array([float(el) for el in range(nbins + 1)])
            bins = np.dot(np.array(xmax / float(nbins)), bins)

            # Plot histogram
            aaa = ax.hist(datt, bins, edgecolor="k", color=col)

            # Add labels
            if xlab is not None:
                ax.set_xlabel(xlab)
            if ylab is not None:
                ax.set_ylabel(ylab)

            # Set x-limits
            ax.set_xlim([0, xmax])
            # Set y-limits and ticks
            ymax = int(math.ceil(ax.get_ylim()[1] / 10) * 10)
            ax.set_yticks([0, ymax / 2, ymax])
            ax.set_ylim([0, ymax + 2])

            # Add strain label
            if lab is not None:
                ax.text(
                    0.95,
                    0.95,
                    lab,
                    ha="right",
                    va="top",
                    transform=ax.transAxes,
                    fontsize=12,
                    weight="roman",
                )

            # Add median line
            ax.axvline(np.mean(datt), color="k", linestyle="dashed", linewidth=5)
            # Add median value label
            # mu = np.median( lens )
            mu = np.mean(datt)
            form = "%.2f"
            mu_str = np.array2string(mu, formatter={"float_kind": lambda mu: form % mu})
            std = np.std(datt)
            std = std / np.sqrt(len(datt))
            std_str = np.array2string(
                std, formatter={"float_kind": lambda std: form % std}
            )
            ax.text(
                0.95,
                0.85,
                r"{0} $\pm$ {1} {2}".format(mu_str, std_str, unit),
                ha="right",
                va="top",
                transform=ax.transAxes,
                fontsize=12,
                weight="roman",
            )
            ax.text(
                0.95,
                0.75,
                r"N = {0}".format(len(datt)),
                ha="right",
                va="top",
                transform=ax.transAxes,
                fontsize=12,
                weight="roman",
            )

    def PlotTrackByStates(self, cols, k=5):
        # Generate figure and axes and set colors

        fig, axs = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
        axsd = {"Poleward": axs[0], "Antipoleward": axs[1], "Inactive": axs[2]}

        k = min([k, len(self.tracks)])
        for idx, track in enumerate(random.sample(self.tracks, k)):

            minis, _ = track.SplitTrack()
            ax = axsd[minis[0].direction]

            # if self.label == "TD" and minis[0].direction == 'Antipoleward':
            # pdb.set_trace()
            # print('1')
            for trk in minis:
                pos_track_rel = trk.CalcPositionRelative()
                ax.plot(
                    trk.time - track.time[0],
                    np.absolute(pos_track_rel[0, :]),
                    linewidth=0.5,
                    color=cols[trk.direction],
                    alpha=0.4,
                )

        # Set x and y limits of subplots
        xl = (0, 0)
        yl = (0, 0)
        for ax in axs:
            xli = ax.get_xlim()
            yli = ax.get_ylim()
            xl = (min([xli[0], xl[0]]), max([xli[1], xl[1]]))
            yl = (min([yli[0], yl[0]]), max([yli[1], yl[1]]))

            # Force x limit
            xl = (xl[0], 400)

        # Legend
        axs[0].plot([], [], label="Poleward", color=cols["Poleward"])
        axs[0].plot([], [], label="AntiPoleward", color=cols["Antipoleward"])
        axs[0].plot([], [], label="Inactive", color=cols["Inactive"])
        axs[0].legend(frameon=False)

        axs[2].set_xlabel("Time (s)")
        # axs[0].set_ylabel(r'Distance from SPB ($\mu m$)')
        axs[1].set_ylabel(r"Distance from SPB ($\mu m$)")
        # axs[2].set_ylabel(r'Distance from SPB ($\mu m$)')
        axs[0].set_ylim(bottom=-0.01, top=yl[1])
        axs[1].set_ylim(bottom=-0.01, top=yl[1])
        axs[2].set_ylim(bottom=-0.01, top=yl[1])
        axs[0].set_xlim(left=-1, right=xl[1])
        # axs[0].set_xlim(right=300) # WT monopolar
        # axs[0].set_xlim(right=400) # KLP5D monopolar
        # axs[0].xaxis.set_ticklabels([])
        # axs[1].xaxis.set_ticklabels([])

        plt.tight_layout()
        fig.savefig("tracks_by_state_{0}.pdf".format(self.label))
        # fig.subplots_adjust(hspace = -0.2)
        plt.close()

    def PlotAllTracks(self, cols):

        fig, ax = plt.subplots(figsize=(4, 3))
        # axsd = {'Poleward': axs[0], 'Antipoleward': axs[1], 'Inactive':axs[2]}

        for idx, track in enumerate(self.tracks):

            minis, _ = track.SplitTrack()
            for trk in minis:
                pos_track_rel = trk.CalcPositionRelative()
                ax.plot(
                    trk.time - track.time[0],
                    np.absolute(pos_track_rel[0, :]),
                    linewidth=0.5,
                    color=cols[trk.direction],
                    alpha=0.6,
                )

        # Set x and y limits of subplots
        # ymax=9
        # xmax=1000
        # ax.set_xlim(left=0.0,right=xmax)
        # ax.set_ylim(bottom=-0.1,top=ymax)
        ax.set_xlim(left=0.0)
        ax.set_ylim(bottom=-0.01)
        ymax = ax.get_ylim()[1]
        xmax = ax.get_xlim()[1]
        ax.set(
            xlabel="Time (s)",
            ylabel="Distance from SPB ($\mu m$)".format(len(self.tracks)),
        )
        # Adding text inside a rectangular box by using the keyword 'bbox'
        plt.text(0.8 * xmax, 0.6 * ymax, "N = {0}".format(len(self.tracks)), fontsize=8)

        # Legend
        ax.plot([], [], label="Poleward", color=cols["Poleward"])
        ax.plot([], [], label="Antipoleward", color=cols["Antipoleward"])
        ax.plot([], [], label="Inactive", color=cols["Inactive"])
        ax.legend()

        plt.tight_layout()
        plt.savefig("tracks_{0}.pdf".format(self.label))
        # fig.subplots_adjust(hspace = -0.2)
        plt.close()


if __name__ == "__main__":
    print("no default implementation")
