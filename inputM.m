% Hard coded data in .m spreadsheet
%
% Sample Workbooks: ./Input_file_MSH.m or Input_file_Romania.m

%EVERYTHING MARKED BY 'E' IS TO/CAN BE EDITED
% Which analysis do you want to perform?
% Pick delay and Qc without kernels: pa=1
% Pick delay and Qc with kernels: pa=2
% Pick delay, kernel-Qc and P/S wave attenuation with the CN method: pa=3

Murat.analysis = 3;%'E'
pa=Murat.analysis;

% INPUT PATHS AND FORMATS
% Creates folders and paths to store results and figures
Murat.paths.label = 'MSH';%'E'
Murat.paths.workingdir = './';%'E'

% Creating a list of files in folder
% Folder containing data
Murat.paths.datadir ='./sac_MSH/*.sac';%'E'
DFolder=Murat.paths.datadir;

% Output figure format 
Murat.figures.format = 'jpeg';%'E'

% Figure visibility - pick 'on' or 'off'
Murat.figures.visibility = 'on';%'E'

%set here the planes you want to inspect - only for pa=3
if pa==3
    Murat.figures.sections=[565000 575000; 5115000 5125000;...
        -2000 -10000];%'E'
end


% INPUT DATA
% Choose between P- (2) and S- (3) waves for peak delay
% RULE: First column in time files always contains origin of events
% RULE: Second column in time files always contains P-wave phase
% RULE: Third column in time files always contains S-wave phase
Murat.data.PorS = 2;

% Work with 1 vertical (1) or 2 horizontal (2) recordings, or with the
% three components(3)
Murat.data.components =  1;%'E'

% To see spectrogram (set>0) of the seismogram: set seisi='number of
% seismogram to see spectrogram'
Murat.data.spectrogram = 0;%'E'

% Central frequency (Hz) - set it according to your spectrograms
Murat.data.centralFrequency = 6;%'E'
cf=Murat.data.centralFrequency;

% Maximum window to pick pick-delays in seconds
Murat.data.maximumPD = 10;%'E'

% Minimum peak delay considering scattering in the area and frequency
Murat.data.minimumPD = 0.1;%'E'

% Lapse time for the start of the window used to measure and calculate the
% normalization energy
Murat.data.startLT = 30;%'E'

% Total coda window length for Qc and normalization
Murat.data.codaWindow = 15;%'E'

% Parameter for smoothing - must be > 2
Murat.data.smoothing = 8;%'E'

% The sped coefficient sets the spectral energy decay of the coda
% wavefield
Murat.data.spectralDecay = 0.5;%'E'

% Number of nodes along WE and SN
% If the origin time is unknown, you can set a theoretichal velocity for
% the whole area and evaluate it from picking. It must be the velocity of
% the phase you are mapping. Velocity is in km/s.
Murat.data.averageVelocity = 0;%'E'

% Name of the SAC variables where zero time and P/S pickings are saved
Murat.data.originTime = [];%'E'
Murat.data.PTime = 'SAChdr.times.a';%'E'
Murat.data.STime = [];%'E'

%start time to measure noise - before P-arrival in seconds - for pa=3
if pa==3
    % Length of the window used to measure noise and direct wave intensity
    Murat.data.bodyWindow = 1;%'E'

    % Start of the window used to measure noise
    Murat.data.startNoise = 5;%'E'
    
    % The minimum coda-to-noise energy ratio for the weighted inversion
    Murat.data.tresholdNoise = 10;%'E'

end

% INPUT GEOMETRY
% Import event origin time and coords of event and station from files
% (Import=1).
% Ideal is to have both in the SAC header (Import=2) and do it
% in lat/long.
evst = 1;%'E'
Murat.geometry.import=evst;

% UTM coordinates of the origin of the model - this must be put
% (1) at least one cell before the westernmost and southernmost
% earthquakes/station for positive latitudes and longitudes
%
% (2) at least one cell after the easternmost and northermost
% earthquakes/station for positive latitudes and longitudes
origin = [538311 5092338 3350];%'E'
Murat.geometry.origin=origin;

%Step of the grid and number of nodes for pd and Qc imaging
nxc = 10;%'E'
Murat.geometry.gridX = nxc;

nyc = 12;%'E'
Murat.geometry.gridY = nyc;

stepg = 4000;%'E'
Murat.geometry.gridStep=stepg;

%Set if in meters (1) or degrees (111)
degorutm = 1;
Murat.geometry.degreesorutm=degorutm;

%Name of the event file if importing from file
namee='even.txt';

%Name of the station file if importing from file
names='staz.txt';

% Inputs necessary for the direct-wave CN attenuation tomography
% EDIT ONLY IF pa=3
if Murat.analysis==3
    
    %set createrays to 1 if you want to create files for rays
    %WARNING: This will create A BIG .mat file
    Murat.geometry.createrays = 0;%'E'
    
    %treshold for the synthetic test, setting layer with top dtreshold
    dtreshold = -3500;%'E'
    Murat.geometry.depthTreshold=dtreshold;
    
    %Rays measured in meters or degrees.    
    if degorutm==1
        Murat.geometry.unity = 1000;
    elseif degorutm==111
        Murat.geometry.unity = 1;
    end
    
    %1D (1) or 3D (3) velocity model
    dimV = 3;%'E'
    
    %name of the velocity model
    namev='modv.txt';%'E'
    % This works in [x,y,z], created for UTM WGS84 and origins to zero.
    % In the new version if a 3D velocity model is unavailable a false 3D
    % is created from iasp91
    if dimV==1
        
        modv1=load(namev); %1D velocity model iasp91
        dend=-34000;%'E' %select maximum depth range
        
        li=length(modv1(:,1));
        modv=zeros(lxy*li,5);
        index=0;
        for i=1:nxc
            for j=1:nyc
                index1=index;
                index=index+1;
                modv(index1*li+1:index1*li+li,1)=XY(index,1)*1000;
                modv(index1*li+1:index1*li+li,2)=XY(index,2)*1000;
                modv(index1*li+1:index1*li+li,3)=...
                    -modv1(1:li,1)*1000;
                modv(index1*li+1:index1*li+li,4)=...
                    modv1(1:li,PorS+1);
            end
        end
        
        modv(:,1)=(modv(:,1)-modv(1,1));
        modv(:,2)=(modv(:,2)-modv(1,2));
        modv(modv(:,3)<dend,:)=[];
        
    elseif dimV==3
    
        modv=load(namev); %3D velocity model from text file
        modv(:,5)=0;
        modv(:,1)=modv(:,1)+origin(1);
        modv(:,2)=modv(:,2)+origin(2);
        
    end
    
    Murat.geometry.modv=modv;
    
    chx=find(modv(:,1)~=modv(1,1),1);
    chy=find(modv(:,2)~=modv(1,2),1);
    resol2x = abs(modv(chx,1)-modv(1,1))/2;
    resol2y = abs(modv(chy,2)-modv(1,2))/2;
    resol2z = abs(modv(2,3)-modv(1,3))/2;
    
    
    %Steps of the velocity models
    passox=max(modv(:,1))-min(modv(:,1));
    passoy=max(modv(:,2))-min(modv(:,2));
    passoz=max(modv(:,3))-min(modv(:,3));
    
    passo=[passox passoy passoz];
    resol=[resol2x resol2y resol2z];
    resol2=min(resol);
    Murat.geometry.resolutionMin=resol2;
    
    %Creates grid for direct waves and check for zeroes
    %Regular step of the gridD for interpolation - half of step of modv
    
    ixD=floor(passox/resol2x);%numer of x,given the step of gridD
    iyD=floor(passoy/resol2y);%numer of y,given the step of gridD
    izD=floor(passoz/resol2z);%numer of depths,given the step of gridD
    
    gridD=zeros(3,max(passo./resol));
    gridD(1,1:ixD)=origin(1):resol2x:origin(1)+passox-resol2x;
    gridD(2,1:iyD)=origin(2):resol2y:origin(2)+passoy-resol2y;
    gridD(3,1:izD)=-origin(3):resol2z:-origin(3)+passoz-resol2z;
    Murat.geometry.gridD=gridD;
    % gridD dimensions
    pvel = zeros(ixD,iyD,izD);
    
    %NUMBER OF X, Y, AND Z LAYERS
    for k=1:izD
        index=0;
        for i=1:ixD
            for j=1:iyD
                index=index+1;
                pvel(i,j,k) = modv(index,4);
            end
        end
    end
    Murat.geometry.pvel=pvel;
end

% INPUT INVERSION
% Size anomaly for testing: twice(2). Might be set to four(4). The input of
% the checkerboard must be always checked visually at the end of
% the process    
Murat.inversion.sizeCheck = 2;%'E'

% Values of attenuation for testing
Murat.inversion.highCheck = 0.02;%'E'
Murat.inversion.lowCheck = 0.001;%'E'

% Seconds of time-windows for non-linear inversion and corresponding
% number, set nonlinear=1 to activate;
Murat.inversion.nonlinear = 1;%'E'

% Uncertainty on Qc estimation
if Murat.inversion.nonlinear==0
    % Minimum R-squared for Qc fitting
    Murat.inversion.fitT = 0;%'E'
elseif Murat.inversion.nonlinear==1
    % Length of smaller time windows to compute compute coda intensity
    Murat.inversion.fitL = 3;%'E'
    Murat.inversion.fitT = 5;%'E'
    % Number of time windows to compute coda intensity
    Murat.inversion.fitN=Murat.data.codaWindow/Murat.inversion.fitL;
    %Grid search - set the minimum, maximum and spacing to search in the
    %parameter space
    m1min = 0;%'E'
    Murat.inversion.minimum=m1min; % minimum inverse Qc allowed
    m1max = 0.01;%'E'
    Murat.inversion.maximum=m1max;% minimum inverse Qc allowed
    total = 1001;%'E'
    Murat.inversion.total = total;% total number of Qc in the interval
    Murat.inversion.fit = (m1min + (m1max-m1min)/(total-1)*(0:total-1))';
end

%% CHECKS AND LOOPS - DO NOT EDIT
% How big are markers for stations and events
Murat.figures.sizeMarker=60;

% Check that user has the Mapping Toolbox installed.
hasMT = license('test', 'map_toolbox');
Murat.figures.hasMT=hasMT;
if ~hasMT
  % User does not have the toolbox installed.
  message =...
      sprintf('No Mapping Toolbox - Maps will not be geo-localised');
end

if evst==1
    
    % Opening the event file. Format is:
    % column (1) = twelve numbers for the origin time of the event
    % (date+time in seconds)
    % column (2) = UTM (WE) or latitude
    % column (3) = UTM (SN) or longitude
    % column (4) = Altitude above sea level in meters
    event = fopen(namee);
    namee= textscan(event,'%s %f %f %f');
    Murat.geometry.nameeven=namee{1};
    Murat.geometry.even=[namee{2} namee{3} namee{4}];
    fclose(event);
    
    % Opening the station file. Format is:
    % column (1) = Name of station (3 characters)
    % column (2) = UTM (WE) or latitude
    % column (3) = UTM (SN) or longitude
    % column (4) = Altitude above sea level in meters
    station=fopen(names);
    names= textscan(event,'%s %f %f %f');
    Murat.geometry.namestation=names{1};
    Murat.geometry.station=[names{2} names{3} names{4}];
    fclose(station);
    
end

%pre-define 2D matrix in space
XY = zeros(floor(nxc)*floor(nyc),2);

% Some figures require moving the point of plotting to half the resolution
halfres=stepg/2;
Murat.geometry.x=origin(1)+halfres:stepg:stepg*(nxc-1)+origin(1)+halfres;
Murat.geometry.y=origin(2)+halfres:stepg:stepg*(nyc-1)+origin(2)+halfres;

index=0;
for i=1:nxc
    for j=1:nyc
        index=index+1;
        XY(index,1:2)=[origin(1)+stepg*i origin(2)+stepg*j];
    end
end
Murat.geometry.map=XY;
lxy=length(XY(:,1));

% Name of the file containing the seismograms and rays
% creates a list of the waveforms and calculates their number.
% RULE: This is necessary for the option Murat.geometry.import = 1, where
% event and station info are in two separate files. In this case,
% the SAC file-name must have a specific format:
% the first 12 characters are the origin time of the event, as per
% event-location file). Characters 18:20 are the name of station as per
% station file.
%
% For Murat.geometry.import = 2, there is no requirement on SAC filenames.
% In the case of pa=3 the names of rays created will be identical to files.

% Get name of files in traces' folder
list = dir(DFolder);
filenames = {list.name}'; %create cell array of file names
filedir = {list.folder}'; %create cell array of file folder
listasac = filenames;
lls = length(listasac);
lung=zeros(lls,1);
lista = filenames; %This will be the name of the file listing rays

for i = 1:lls
    listasac{i,1} = cat(2,filedir{i},'/',filenames{i});
    %Here it takes the event/station name. It is necessary to adapt the
    %numbers to where even name and station name are
    if evst==1
        li = lista{i,1};
        li1 = cat(2,li(1:12),li(18:20));
        lista{i,1} = li1;
    end
end

Murat.paths.listasac=listasac;
Murat.paths.lista=lista;

if exist(cat(2,'./',Murat.paths.label),'dir')~=7
    mkdir(Murat.paths.label)
end

%clearvars -except Murat
save('Murat.mat','Murat');