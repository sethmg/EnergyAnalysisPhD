%% Clean house and read in the data files
clc
clear all
close all

[inDTForceFile,inDTForcePath] = uigetfile('/media/Test_Data/DT_*.mat','Please select an drop tower data file');
if exist([inDTForcePath,'DT_',inDTForceFile(4:9),'_TEMA_Displacement_Processed_filtfilt.mat'],'file') % check if the expected drop tower file is with the instron file, if not ask user to locate
    inDTDispFile = ['DT_',inDTForceFile(4:9),'_TEMA_Displacement_Processed_filtfilt.mat'];
    inDTDispPath = inDTForcePath;
else
    [inDTDispFile,inDTDispPath] = uigetfile([inInsPath,'DT_*TEMA*.mat'],'Please select a DT displacement file');
end

sTime = inputdlg('Please enter the time of the first image in ms.','Time sync',1,{'0'});
sTime = str2num(sTime{1});

% read the input files
load([inDTForcePath,inDTForceFile]);
dTTime = time;                                                              % rename the drop tower signal data time vector to prevent variable name clash
clear time;
load([inDTDispPath,inDTDispFile]);

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

%% Calcualte the energy to failure in the drop tower
figure
plot(DTSixAInterp(dispDefinedRange(1):dispDefinedRange(2)));
[x,y] = ginput(1);
[DT_maxF,DT_maxFI] = max(DTSixAInterp(floor(x-10):ceil(x+10)));
DT_maxFI = floor(x-10)+DT_maxFI;

indexes = dispDefinedRange(1):dispDefinedRange(1)+DT_maxFI-1;
DT_energy = 0;
for i = indexes
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
    fprintf(outFileID,'Specimen\tDroptower to Fracture(J)\tInstron to max force (J)\tDrop tower to max instron force (J)\n');
else
    outFileID = fopen('EnergyComparisons.txt','a+');
end
specimenName = inDTDispPath(strfind(inDTDispPath,'H1'):strfind(inDTDispPath,'H1')+5);       % get the specimen name from the input path
fprintf(outFileID,'%s\t%24.5f\t%24.5f\t%35.5f\r\n',specimenName,DT_energy,0,0);
fclose(outFileID);


