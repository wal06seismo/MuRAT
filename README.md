# MuRAT
A muti-resolution seismic attenuation tomography code - currently in its 2.2 release

Getting Started
-------

MuRAT is a Matlab Package for Seismic Attenuation Tomography at multiple Earth scales using Body and Coda Waves. 


*Features*
--------

* MuRAT is a complete software package for measuring seismic attenuation, scattering, and absorption from passive and active data, and map 2D and 3D variations of these parameters in space.

* MuRAT also comes in Python and C++ versions (under development) to provide a fully-integrated solution for seismic attenuation imaging. 

* MuRAT1.0 was first developed by Luca De Siena (Johannes Gutenberg University, Mainz, Germany) during his PhD at the INGV-Osservatorio Vesuviano, Italy, and published in 2014 while he was research assistant at the Westfälisches Wilhelms Universität, Münster.

* MuRAT2.0 is the result of the activity of the Volcano Imaging group, led by De Siena during his stint as permanent Lecturer at the University of Aberdeen, UK.

* The group of active users (providing questions, feedback, snippets of code) comprises De Siena's PhD students at the University of Aberdeen. It includes PhD students he co-supervises internationally. 

*History*
-------

* 2006-2010: MuRAT (at the time named "Multi-scale reasonable attenuation tomography analysis") is built using Matlab, c++, csh and Fortran codes mostly developed at INGV-Osservatorio.

* 2010-2013: MuRAT1.0 is developed as a 3D direct-wave attenuation imaging Matlab-only code with the contribution of Christine Thomas and Richard Aster.

* 2014: MuRAT1.0 is published in De Siena et al. 2014, JVGR, with two sample datasets (Mount St. Helens and Vesuvius) [JVGR article](https://www.sciencedirect.com/science/article/abs/pii/S0377027314000961)

* 2019: MuRAT2.0 is released including 2D scattering/absorption mapping and kernel-based inversion, and re-branded Multi-Resolution Seismic Attenuation Tomography. [GitHub Repository] (https://github.com/LucaDeSiena/MuRAT)

*Documentation*
-------------

The full documentation for MuRAT2.0 is under construction. This README file and the linked internet sites are to be used as a reference. 

*Installation and running*
------------

SYSTEM: The program works on: 1) Macbook pro with High Sierra, Matlab R2017a; 2) Ubuntu 16.04.6 LTS, Matlab R2017a.

Necessary Toolboxes: Signal Processing, Curve Fitting, Image Processing (compulsory) and Mapping (optional, for geolocalisation - Romania example).

Two sample datasets (Mount St. Helens and Romania) can be downloaded at https://doi.pangaea.de/10.1594/PANGAEA.893893 to test the code. The third sample input (Pollino) in the working directory works with a dataset that may be requested to Luca De Siena.

The current version works following these steps:

1. Download the package at https://github.com/LucaDeSiena/MuRAT.

2. Download the two sample datasets at https://doi.pangaea.de/10.1594/PANGAEA.893893. Unzip the MSH (Mount St. Helens) and Romania datasets and put the folders in Murat-master folder. Delete the zipped files.

3. MuRAT works with a single input file in .m format. Build the Input file starting from the samples (Input_MSH.m, Input_Romania.m, Input_Pollino and Input_OlDoinyo). Read the instructions presented in the next sections.

4. Run MuRAT2_2.m.

5. After the L-curve plots are visible, a prompt asks for the smoothing parameter. Input them based on the plot.

6. When applying MuRAT to your dataset, use one of the Input files as a template. Read the comments attentively and edit only the required parameters. Start with a “pa=1” or “pa=2” analysis. “pa=3” is time-consuming.

*Instructions - the input file*
------------

GENERAL - CHOICES

------------

**Murat.analysis**

*Available indexes: 1, 2, or 3 - Default -> 1*

Analysis = 1 maps pick-delay and Qc without kernels, as *De Siena et al. 2016 (EPSL)*

Analysis = 2 inverts for pick-delay and Qc with analytic kernels

Analysis = 3 inverts for pick-delay and Qc with Pacheco-Snieder kernels.

------------
**Murat.geometry.kernelTreshold**

Treshold to reduce computational time for the Pacheco-Snieder kernels. It divides the inversion grid cell by the treshold.

Higher number=longer computational time.

------------
**Murat.data.PorS**

*Available indexes: 2 or 3 - Default -> 2*

Index one in the time files contains the origin time of events. Index two in the time files contains the P-wave phase. Index three in the time files contains the S-wave phase.

------------
**Murat.data.centralFrequency**

This input is the frequency where the analyses will be carried on, depending on the data recordings.

------------

PATHS AND FIGURES

------------
**Murat.paths.workingdir**

Directory where MuRAT.m is contained - *Default = ./*

------------
**Murat.paths.datadir**

Directory containing the SAC files - *Default = ./sac_Label*

------------
**Murat.paths.label**

Name of the folder containing output files and figures. It will appear in the working directory

------------
**Murat.paths.originTime**

This input is the origin-time of the event. The best practice is to save it in the SAC file as *o* variable -> *SAChdr.time.o*. If this is unavailable set it as *Default -> \[]*

------------
**Murat.paths.PTime**

This input is the P-wave time of the event; this parameter is compulsory. The best practice is to save it in the SAC file as *a* variable -> *SAChdr.time.a*

------------
**Murat.paths.STime**

This input is the S-wave time of the event. The best practice is to save it in the SAC file as *t0* variable -> *SAChdr.time.t0*. If this is unavailable set it as *Default -> []*

------------
**Murat.figures.format**

Figures' output format as per Matlab functions - *Default -> jpeg*

------------
**Murat.figures.visibility**

Set visibility during computation, set on or off - *Default -> on*

------------
DATA DRIVEN CHOICES - WAVEFORMS

------------
**Murat.data.components**

*Available indexes: 1, 2 or 3 - Default -> 1*

MuRAT works with one vertical (1) or two horizontal (2) recordings, or with the three components(3) of motion.

------------
**Murat.data.smoothing**

*Default -> 8*

Parameter used to smooth envelopes - it will be the parameter divided by the central frequency

------------
**Murat.data.maximumPD**

This input is the maximum time from picking where the maximum of the envelope will be searched for in the peak delay analysis.

------------
**Murat.data.minimumPD**

This input is the minimum time from picking where the maximum of the envelope will be searched for in the peak delay analysis.

------------
**Murat.data.startLT**

This input is the starting time for the coda window, used for both the coda attenuation and the total-attenuation analyses.

------------
**Murat.data.codaWindow**

This input is the length of the coda window for the coda attenuation analysis.

------------
**Murat.data.spectralDecay**

*Available indexes: 0.5, 1 or 1.5 - Default -> 0.5*

This input is the expected envelope decay for surface or body waves in 2D and 3D.

------------
DATA DRIVEN CHOICES - VELOCITY

------------
**Murat.geometry.availableVelocity**

*Available indexes: 0, 1 or 2 - Default -> 0*

This input selects the velocity model to be used for ray-tracing. For now, only index 0 is valid

------------
**Murat.data.averageVelocityP**

This input is the average P-wave , in the case that the origin-time is absent from recordings or a velocity model is unavailable - *Default -> 6*

Also used to compute Vp/Vs

------------
**Murat.data.averageVelocityS**

This input is the average S-wave , in the case that the origin-time is absent from recordings or a velocity model is unavailable - *Default -> 6*

Also used to compute Vp/Vs

-------------
**namev**

The name of the file containing the 3D velocity model. Only available for **Murat.geometry.availableVelocity=3** - *Default -> []*

------------
GEOMETRY

------------
**Murat.geometry.import**

*Available indexes: 1 or 2 - Default -> 2*

It is necessary to define source and station locations for mapping: the user can import event origin-time and coordinates of events and station from an external .txt file (1) or directly from the SAC files (2). Index one is the original format of MuRAT 1.0 and requires even.txt and staz.txt files as per sample files. Index two is the ideal format, storing event and station information as latitude and longitude in the SAC header. The user needs to fill the following variables in the header:

*Event*: Name of event - *SAChdr.event.kevnm*; Station latitude - *SAChdr.event.evla*; Station longitude- *SAChdr.event.evlo*; Station depth - *SAChdr.event.evdp*

*Station*: Name of station - *SAChdr.station.kstnm*; Station latitude - *SAChdr.station.stla*; Station longitude- *SAChdr.station.stlo*; Station elevation - *SAChdr.station.stel*

-------------
**Murat.geometry.origin**

The origin point of the 2D grid for imaging is either in UTM or lat/long with depths in meters or kilometers.

-------------
**Murat.geometry.end**

The end point of the 2D grid for imaging is either in UTM or lat/long with depths in meters or kilometers.

-------------
**Murat.geometry.gridX** and **Murat.geometry.gridY**

The number of layers desired in the X and Y directions.

-------------
**Murat.geometry.degreesorutm** and **Murat.geometry.unity**

Can be set to either 1 (km) or 111 (degrees). The second is used to met meters or km. (1000 or 1)

------------
INVERSION

-------------
**Murat.inversion.sizeCheck**

*Available indexes: 2 or 4*

This input sets the dimension of the checkerboard anomalies to either twice or four times the step of the imaging grid.

-------------
**Murat.inversion.sizeCheck**

*Available indexes: 2 or 4*

This input sets the dimension of the checkerboard anomalies to either twice or four times the step of the imaging grid.

-------------
**Murat.inversion.highCheck** and **Murat.inversion.lowCheck**

This input sets the high-attenuation and low-attenuation values (as inverse Q and inverse Qc) for the checkerboard tests.

-------------
**Murat.inversion.nonlinear**

*Available indexes: 0 or 1 - Default -> 0*

With this input, the user decides to inverts coda attenuation with either a linearised (0) or non-linear (1) approach.

-------------
**Murat.inversion.fitT**

*In the linear case*
The minimum accepted inverse-Qc uncertainty used to weight the coda attenuation data.

*In the non-linear case*
The number of smaller time windows in which we divide the coda window - to be set together with next entry.

-------------
**Murat.inversion.fitL**

These inputs correspond to the length of the smaller windows in which we divide the coda window.

Example: for a total coda window of 15 seconds, the user can set either Murat.inversion.fitT=3 windows of Murat.inversion.fitL=5 seconds each or Murat.inversion.fitT=5 windows of Murat.inversion.fitL=3 seconds each. Only active for the non-linear inversion.

-------------
**Murat.inversion.minimum**, **Murat.inversion.maximum**, and **Murat.inversion.total**

With these inputs, the user searches for the Qc minimising the inversion from a minimum inverse Q (e.g. 0) to a maximum inverse Qc (e.g. 0.01). The total number of Qc are equally-spaced and defined between minimum and maximum. To be set to appropriate parameters for the non-linear inversion.

*Instructions - the output files*
------------

All the output files (.txt) and figures (in the format defined by the user) are stored in the **Murat.paths.label** folder, created in the **Murat.paths.workingdir**. The first three columns of the output files correspond to WE, SN, and depth. The fourth column is the mapped parameter. 

------------
*peakdelay.txt*

The variations of Log10 Peak delay from the average - the parameter used to map scattering attenuation. The depth is set to -1000 m.

The fifth and sixth columns correspond to the input and output of the checkerboard test, respectively. The seventh column is the first ray.

------------
*Qc.txt*

Variations of inverse Qc, alias coda attenuation - the parameter used to map absorption. The depth is set to -1000 m.
    
The fifth and sixth columns correspond to the input and output of the checkerboard test, respectively. The seventh column is the first 2D kernel.

------------

All the figures (in the **Murat.figures.format** defined by the user) are stored in the **Murat.paths.label** folder, created in the **Murat.paths.workingdir**. If the Mapping Toolbox is available and the coordinates are in latitude and longitude, the figures will show coastlines.

------------
*Rays.Figures format*

A figure to plot rays either in 2D (corresponding to peak-delay sensitivity) - it will show rays in the reference system of the pre-defined grid with sources and stations used.

------------
*Qc_sensitivity.Figures format* - Only available for **Murat.analysis = 2 or 3**

A figure to plot the source-station normalised kernel for the first source-station pair. It will show the sensitivity in the reference system of the pre-defined grid.


------------
*Qc_Peak_Delay.Figures format*

A figure to evaluate the appropriate peak-delay and coda inputs. The upper panel should show a constant inverse Qc with travel time. The lower panel should show an increasing peak delay with travel time. Red dots correspond to outliers. For Qc, the maximum uncertainty to define outliers is set as **Murat.inversion.fitT**, for the non-linear it is pre-defined. For the peak-delay, we define as outliers all the values over twice the standard deviation from the average.


------------
*Lc_Qc.Figures format and Lc_Lc_PeakDelay.Figures format*

These are the only plots appearing if **Murat.figures.visibility = 0** during computation. They show the L-curves corresponding to the peak-delay (ray sensitivity) and coda-attenuation inversions. After they appear, a prompt asks which damping parameter the user wants to pick.

------------
*Picard_Qc.Figures format and Picard_PeakDelay.Figures format*

These plots show the result of the Picard analysis, necessary to evaluate how many of the inversion parameters are correctively solved in the coda-attenuation and peak-delay inversions, respectively. The two figures do not appear during computation.

------------
*Peak_delay_map.Figures format* and *Qc_map.Figures format*

These plots show the result of the peak-delay and Qc mapping in the grid's reference system.

------------
*Qc_checkerboard.Figures format* and *PD_checkerboard.Figures format

These plots show the result of the checkerboard tests for the Qc and peak delay mapping in the grid's reference system.

------------
*Parameter_space_variations.Figures format*

The plot shows the separation of the scattering and absorption parameters in their parameter space. Grey dots correspond to parameters too near to the average to be interpreted as scattering or absorption variations - the threshold is pre-defined at 5% of the maximum variation of each parameter. Red = High scattering and absorption; Cyan = High scattering and low absorption; Orange = Low scattering and high absorption; Green = Low scattering and absorption.

------------
*Parameter_map.Figures format*

The plot of the parameters in space.

*Citing MuRAT*
------------

If you use MuRAT for your research and publications, please consider mentioning the GitHub internet site and citing the first release of the code, published as:

*De Siena, L., C. Thomas, and R. Aster. "Multi-scale reasonable attenuation tomography analysis (MuRAT): An imaging algorithm designed for volcanic regions." Journal of Volcanology and Geothermal Research 277 (2014): 22-35.*


*Disclaimer*
------------

Although we have cross-checked the whole code, we cannot warranty it is exempt of bugs. The package is provided as-is, we will neither be held responsible for any use you make of it nor for the results and conclusions you may derive using MuRAT.

*Licence*
------------

MuRAT is released under EUPL v1.1
