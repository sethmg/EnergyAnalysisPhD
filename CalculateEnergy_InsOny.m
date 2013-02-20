%% Clean house and read in the data files
clc
clear all
close all

[inInsFile,inInsPath] = uigetfile('/media/Test_Data/Ins_*.mat','Please select an Instron data file');

% read the input files
load([inInsPath,inInsFile]);
insTime = time;                                                             % rename the instron time vector to prevent variable name clash with other input file variables


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

%% Write the data to an output file
if ~exist('EnergyComparisons.txt','file')   % if the data file does not exist, create and initalize it with the header line.
    outFileID = fopen('EnergyComparisons.txt','w');
    fprintf(outFileID,'Specimen\tDroptower to Fracture(J)\tInstron to max force (J)\tDrop tower to max instron force (J)\n');
else
    outFileID = fopen('EnergyComparisons.txt','a+');
end
specimenName = inInsPath(strfind(inInsPath,'H1'):strfind(inInsPath,'H1')+5);       % get the specimen name from the input path
fprintf(outFileID,'%s\t%24.5f\t%24.5f\t%35.5f\r\n',specimenName,0,Ins_energy,0);
fclose(outFileID);


