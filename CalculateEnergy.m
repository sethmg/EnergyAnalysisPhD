%% Clean house and read in the data files
clc
clear all
close all

[inInsFile,inInsPath] = uigetfile('/media/Test_Data/Ins_*.mat','Please select an Instron data file');
if exist([inInsPath,'DT_',inInsFile(5:10),'_Processed_filtfilt.mat'],'file')         % check if the expected drop tower file is with the instron file, if not ask user to locate
    inDTForceFile = ['DT_',inInsFile(5:10),'_Processed_filtfilt.mat'];
    inDTForcePath = inInsPath;
else
    [inDTForceFile,inDTForcePath] = uigetfile([inInsPath,'DT_*.mat'],'Please select a DT force file');
end
if exist([inInsPath,'DT_',inInsFile(5:10),'_TEMA_Displacement_Processed_filtfilt.mat'],'file') % check if the expected drop tower file is with the instron file, if not ask user to locate
    inDTDispFile = ['DT_',inInsFile(5:10),'_TEMA_Displacement_Processed_filtfilt.mat'];
    inDTDispPath = inInsPath;
else
    [inDTDispFile,inDTDispPath] = uigetfile([inInsPath,'DT_*TEMA*.mat'],'Please select a DT displacement file');
end

sTime = inputdlg('Please enter the time of the first image in ms.','Time sync',1,{'0'});
sTime = str2num(sTime{1});

% read the input files
load([inInsPath,inInsFile]);
insTime = time;                                                             % rename the instron time vector to prevent variable name clash with other input file variables
clear time;
load([inDTForcePath,inDTForceFile]);
dTTime = time;                                                              % rename the drop tower signal data time vector to prevent variable name clash
clear time;
load([inDTDispPath,inDTDispFile]);

%% Calculate Instron energy to the max instron force
[Ins_maxF,Ins_maxFI] = max(-force);
Ins_energy = 0;
for i = 1:Ins_maxFI-1
    if displacement(i+1) > displacement(i) % if the crosshead is moving up, don't count the energy
        continue;
    end
    forceA = -force(i);
    forceB = -force(i+1);
    dispA = -displacement(i)/1000; % mm to m
    dispB = -displacement(i+1)/1000;
    Ins_energy = Ins_energy + (dispB-dispA)*mean([forceA forceB]);
end
% Ins_energy = trapz(-displacement(1:Ins_maxFI)/1000,-force(1:Ins_maxFI)); % negative values are because the displacement and force are compressive

%% Interpolated DT displacement and force into a common time spacing
DTTimeInterp = linspace(-200,500,10000);                                    % time for the interpolated drop tower
if find(isnan(TrackedImpacFilt(:,1))==0,1,'last') > find(isnan(TrackedTrochFilt(:,1))==0,1,'last')
    indexes = 1:find(isnan(TrackedTrochFilt(:,1))==0,1,'last');
else
    indexes = 1:find(isnan(TrackedImpacFilt(:,1))==0,1,'last');
end

DTImpactorDispInterp = interp1(timeDisp(indexes)+sTime,TrackedImpacFilt(indexes,1)/1000-TrackedImpacFilt(1,1)/1000,DTTimeInterp);   % interpolate the impactor displacement into new time vector
DTTrochDispInterp = interp1(timeDisp(indexes)+sTime,TrackedTrochFilt(indexes,1)/1000-TrackedTrochFilt(1,1)/1000,DTTimeInterp);      % interpolate the trochanter displacement into new time vector
DTForceInterp = interp1(dTTime(1:length(oneAxis)),oneAxis,DTTimeInterp);                       % interpolate the single axis load cell into the new time vector
DTSixAInterp = interp1(dTTime(1:length(sixAxis(:,3))),sixAxis(:,3),DTTimeInterp);                   % interoplate the six axis load cell into the new time vector
dispDefinedRange = [find(isnan(DTTrochDispInterp)==0,1,'first') find(isnan(DTTrochDispInterp)==0,1,'last')];    % find the region of the time vector where displacement is defined

%% Calculate the DT energy to max intron force
% find the index of the max instron force
DT_insMaxFI = find(DTSixAInterp(dispDefinedRange(1):dispDefinedRange(2)) > Ins_maxF,1,'first')-1;
indexes = dispDefinedRange(1):dispDefinedRange(1)+DT_insMaxFI-1;
DT_insEnergy = 0;
for i = indexes
    if DTTrochDispInterp(i) > DTTrochDispInterp(i+1) % do not count rebound ernergy
        continue;
    end
    forceA = DTSixAInterp(i);
    forceB = DTSixAInterp(i+1);
    dispA =  DTTrochDispInterp(i);
    dispB =  DTTrochDispInterp(i+1);
    DT_insEnergy = DT_insEnergy + (dispB-dispA)*mean([forceA forceB]);
end

%% Calcualte the energy to failure in the drop tower
figure
plot(DTSixAInterp(dispDefinedRange(1):dispDefinedRange(2)));
[x,y] = ginput(1);
[DT_maxF,DT_maxFI] = max(DTSixAInterp(floor(x-10):ceil(x+10)));
DT_maxFI = floor(x-10)+DT_maxFI;
% [DT_maxF,DT_maxFI] = max(DTForceInterp(dispDefinedRange(1):dispDefinedRange(2)));
indexes = dispDefinedRange(1):dispDefinedRange(1)+DT_maxFI-1;
DT_energy = 0;
for i = indexes
    if isnan(DTTrochDispInterp(i+1))
        continue;
    end
    if DTTrochDispInterp(i) > DTTrochDispInterp(i+1) % do not count rebound ernergy
        continue;
    end
    forceA = DTSixAInterp(i);
    forceB = DTSixAInterp(i+1);
    dispA = DTTrochDispInterp(i);
    dispB = DTTrochDispInterp(i+1);
    DT_energy = DT_energy + (dispB-dispA)*mean([forceA forceB]);
end
%% Write the data to an output file
if ~exist('EnergyComparisons.txt','file')   % if the data file does not exist, create and initalize it with the header line.
    outFileID = fopen('EnergyComparisons.txt','w');
    fprintf(outFileID,'Specimen\tDroptower to Fracture(J)\tInstron_to_max_force_(J)\tDrop_tower_to_max_instron_force_(J)\n');
else
    outFileID = fopen('EnergyComparisons.txt','a+');
end
specimenName = inDTDispPath(strfind(inDTDispPath,'H1'):strfind(inDTDispPath,'H1')+5);       % get the specimen name from the input path
fprintf(outFileID,'%s\t%24.5f\t%24.5f\t%35.5f\r\n',specimenName,DT_energy,Ins_energy,DT_insEnergy);
fclose(outFileID);


