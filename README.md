# KymoAnalysis

KymoAnalysis is a tool that allows users to analyze tracks extracted via KymoAnnotate (https://github.com/saadjansari/KymoAnnotate.git)

This is useful for analysis of protein movement in fluorescent microscopy.

## Prerequisites

Before you begin, ensure you have met the following requirements:
* You have `Python3` (this was tested on Python 3.9.5)
* You have installed `git`.
* You have a `MacOS` machine (not tested on `Linux/Windows`).

## Installing KymoAnalysis

To install KymoAnalysis, follow these steps:

```
git clone https://github.com/saadjansari/KymoAnalysis.git
```

## Using KymoAnalysis

To run KymoAnalysis, copy an example config.yaml file and then run KymographAnalysis.py.
```
cd KymoAnalysis
cp ExampleConfigFiles/config_bipolar_wt.yaml config.yaml
python3 KymographAnalysis.py
```

This should create a new folder named `wt_bipolar` which contains the analyzed results.

This repo has config and data files for:
* Bipolar wt
* Bipolar 989TD
* Monopolar wt
* Monopolar klp5D

## How does it work?

Kymograph analysis involves processing the saved results from KymoAnnotate. 
For a given strain (e.g. bipolar wt), analysis proceeded as follows:
1. Read results file to yield SPB and cut7 tracks (config file specifies data path).
2. Track merging. We started by merging tracks that were deemed to belong to the same global track. For example, if track 2 that began within some time interval of track 1 ending (in our case, twice the timestep), and the start position of track 2 was close to the end position of track 1 (in our case, 150 nm), then we merged the tracks (Fig. S2).
3. Track splitting. After merging, we have a number of tracks that exhibit both poleward and antipoleward motion. Here, a track is split into poleward, antipoleward, and paused segments (Fig. S3) based on the average local velocity along the track. A cutoff value of 3 nm/sec is used. Points with velocity less than the cutoff are assigned as poleward, while points with velocity greater than the cutoff are assigned as antipoleward. The remained points are assigned as paused. With this threshold, we can split a track into poleward, antipoleward and paused segments.
4. Tracks are analyzed with reference to the position of the SPB. To do that, we assigned each track to it’s nearest SPB. For monopolar spindles, there was only one pole and assignment was straightforward. For bipolar spindles with two SPBs, we picked the SPB that was closest to the track start position.
5. Track position and velocity was measured relative to the closest SPB.
6. For each kymograph, the frequency of directional events was calculated by counting the number of poleward and antipoleward tracks and dividing by the total observation time in the kymograph.
7. Switching frequency was calculated by counting the number of switches out of a state (poleward and antipoleward) divided by the lifetime of tracks in that state. For example, switching frequency for exiting the poleward state was calculated by dividing the number of transitions out of the poleward state by the total lifetime of all poleward tracks.


## Contributors

Thanks to the following people who have contributed to this project:

* [@saadjansari](https://github.com/saadjansari) 📖


## Contact

If you want to contact me you can reach me at saadjansari@gmail.com.

## License

This project uses the following license: [MIT](https://opensource.org/licenses/MIT).
