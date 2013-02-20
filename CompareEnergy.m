%%
clc
close all
clear all

if ~exist('EnergyComparisons.txt','file')
    [inEnergyFile,inEnergyPath] = uigetfile('EnergyComparisons.txt','Please select the energy comparisons file');
else
    inEnergyFile = 'EnergyComparisons.txt';
    inEnergyPath = [pwd,filesep];
end

inEnergyFileID = fopen([inEnergyPath,inEnergyFile],'r');
C = textscan(inEnergyFileID,'%s\t%f\t%f\t%f\t%f\t%s','headerlines',1);
Specimens = C{1};
DT_energy = C{2};
Ins_energy = C{3};
DT_insEnergy = C{4};
DXA = C{5};
OP_status = C{6};

insValid = logical(Ins_energy);                     % indexes for specimens with instron data

dtValid = logical(DT_energy);                       % indexes for specimens with dt data

bothValid = logical(Ins_energy.*DT_energy);   % indexes for specimens with data for DT and Ins;

delta_energy = (DT_insEnergy - Ins_energy).*bothValid;

%%  Create some plots
 plotPosition = get(0,'screenSize');                                         % get the screen size
%plotPosition = [1684           4        1674         968];

% Instron Energy to 50% predicted failure load
f1H = figure(1);
a1H = axes;
set(f1H,'position',plotPosition,'paperpositionmode','auto');                                           % set the figure to fill the screen
plot(a1H,DXA(insValid),Ins_energy(insValid),'ko','markersize',15,'linewidth',2,'markerfacecolor','r')
grid
set(a1H,'fontname','times','fontsize',40);
xlabel('Total aBMD by DXA (g/cm^2)','fontname','times','fontsize',45);
ylabel('Energy (J)','fontname','times','fontsize',45);
title('Energy to 50% Predicted Failure in the Instron','fontname','times','fontsize',50);
% Is there a linear correlation for Energy to 50% predected Fx in the Ins and DXA?
[insFitObj,insGOF] = fit(DXA(insValid),Ins_energy(insValid),'poly1');
% if there is a correlation (I'm using r^2 > 0.2), plot it.
if insGOF.rsquare > .2
    hold
    plot(sort(DXA(insValid)),feval(insFitObj,sort(DXA(insValid))),'linewidth',2)
end


% DT energy to 50% perdicted failure load (compare to instron energy)
f2H = figure(2);
a2H = axes;
set(f2H,'position',plotPosition,'paperpositionmode','auto');
plot(a2H,DXA(bothValid),DT_insEnergy(bothValid),'ko','markersize',15,'linewidth',2,'markerfacecolor','r')
grid
set(a2H,'fontname','times','fontsize',40);
xlabel('Total aBMD by DXA (g/cm^2)','fontname','times','fontsize',45);
ylabel('Energy (J)','fontname','times','fontsize',45);
title('Energy to 50% Predicted Failure in the Drop Tower','fontname','times','fontsize',50);
% Is there a linear correlation for energy to 50% predicted Fx in the DT and DXA?
[dtInsFitObj,dtInsGOF] = fit(DXA(bothValid),DT_insEnergy(bothValid),'poly1');
% if there is a correlation (I'm using r^2 > 0.2), plot it.
if dtInsGOF.rsquare > .2
    hold
    plot(DXA(bothValid),feval(dtInsFitObj,DXA(bothValid)),'linewidth',2)
end

% difference between DT and instron energy
f3H = figure(3);
a3H = axes;
set(f3H,'position',plotPosition,'paperpositionmode','auto')
plot(a3H,DXA(bothValid),delta_energy(bothValid),'ko','markersize',15,'linewidth',2,'markerfacecolor','r')
grid
set(a3H,'fontname','times','fontsize',40);
xlabel('Total aBMD by DXA (g/cm^2)','fontname','times','fontsize',45);
ylabel('Energy (J)','fontname','times','fontsize',45);
title('\Delta Energy to 50% Failure in Instron and Drop Tower','fontname','times','fontsize',50);
% Is there a linear correlation for the difference in energy to 50% predicted Fx in the instron and DT and DXA?
[diffFitObj,diffGOF] = fit(DXA,delta_energy,'poly1');
if diffGOF.rsquare > 0.2
    hold
    plot(DXA(bothValid),feval(diffFitObj,DXA(bothValid)),'linewidth',2)
end


% a different way to look at the change in energy
f4H = figure(4);
a4H = axes;
set(f4H,'position',plotPosition,'paperpositionmode','auto');
hold on
plot(a4H,DXA(bothValid),Ins_energy(bothValid),'ko','linewidth',2,'markersize',15);
for i = 1:length(Specimens)
    if bothValid(i) == 0
        continue;
    end   
    

    
    if delta_energy(i) > 0
        plot(a4H,DXA(i),DT_insEnergy(i),'^','linewidth',2,'markersize',15,'markeredgecolor','g','markerfacecolor','g');
    else
        plot(a4H,DXA(i),DT_insEnergy(i),'v','linewidth',2,'markersize',15,'markeredgecolor','r','markerfacecolor','r');
    end
    
    if delta_energy(i) > 0
        plot(a4H,[DXA(i),DXA(i)],[Ins_energy(i),DT_insEnergy(i)],'g-','linewidth',4);
    else
        plot(a4H,[DXA(i),DXA(i)],[Ins_energy(i),DT_insEnergy(i)],'r-','linewidth',4);
    end
end
grid
set(a4H,'fontname','times','fontsize',40);
xlabel('Total aBMD by DXA (g/cm^2)','fontname','times','fontsize',45);
ylabel('Energy (J)','fontname','times','fontsize',45);
title('Energy to 50% Failure in the Instron and Drop Tower','fontname','times','fontsize',50);


% Energy to fracture in the DT
f5H = figure(5);
a5H = axes;
set(f5H,'position',plotPosition,'paperpositionmode','auto');
plot(a5H,DXA(dtValid),DT_energy(dtValid),'ko','markersize',15,'linewidth',2,'markerfacecolor','r')
grid
set(a5H,'fontname','times','fontsize',40);
xlabel('Total aBMD by DXA (g/cm^2)','fontname','times','fontsize',45);
ylabel('Energy (J)','fontname','times','fontsize',45);
title('Energy to Max Force in the Drop Tower','fontname','times','fontsize',50);
% Is there a linear correlation for the energy to fracture and DXA in the DT?
[dtFitObj,dtGOF] = fit(DXA,DT_energy,'poly1');
if dtGOF.rsquare > 0.2
    hold
    plot(DXA(dtValid),feval(dtFitObj,DXA(dtValid)),'linewidth',2)
end

% Bland Altman plot of energy difference
meanEnergy = (DT_insEnergy+Ins_energy)/2;
[delta_mean,delta_sigma,delta_meanci,delta_sigmaci] = normfit(delta_energy(bothValid));
f6H = figure(6);
a6H = axes;
set(f6H,'position',plotPosition,'paperpositionmode','auto');
plot(a6H,meanEnergy(bothValid),delta_energy(bothValid),'ko','markersize',15,'linewidth',2,'markerfacecolor','r')
grid
hold on
xlimits = xlim;
plot(a6H,xlimits,[delta_mean delta_mean],'k','linewidth',3);
plot(a6H,xlimits,[delta_meanci(1) delta_meanci(1)],'color',[.7 .7 .7],'linewidth',3)
plot(a6H,xlimits,[delta_meanci(2) delta_meanci(2)],'color',[.7 .7 .7],'linewidth',3)
xlim(xlimits);
set(a6H,'fontname','times','fontsize',40);
xlabel('Average Energy to 50% Predicted Fracture','fontname','times','fontsize',45);
ylabel('Energy Difference (J)','fontname','times','fontsize',45);
title('Bland-Altman Plot of Energy Difference','fontname','times','fontsize',50);


% %% Some simple stats to see if there are correlations
% Is the energy to 50% perdicted Fx different between DT and Ins?
[real_delta_energy,p_delata_energy,ci_delta_energy,stats_delata_energy] = ttest2(DT_insEnergy(bothValid),Ins_energy(bothValid));
% anova1([DT_insEnergy(bothValid),Ins_energy(bothValid)]);
% 
% % Is the difference between DT and Ins to 50% perdicted FX different based
% % on osteoporosis status?
[p_dxa,table_dxa,stats_dxa] = anova1(delta_energy(bothValid),OP_status(bothValid),'off');

%% Save the energy plots
% outPlotInsEng = 'Ins_energy';
% outPlotDT50Eng = 'DT_50%_energy';
% outPlotDeltaEng = 'Delta_energy';
% outPlotDTInsEng = 'DT_Ins_energy';
% outPlotDTEng = 'DT_MaxF_energy';
% outPlotBlandAlt = 'Delta_BlandAltman';
% 
% 
% 
% print(f1H,'-dpng','-r100',[outPlotInsEng,'_lowRes.png'])
% print(f2H,'-dpng','-r100',[outPlotDT50Eng,'_lowRes.png'])
% print(f3H,'-dpng','-r100',[outPlotDeltaEng,'_lowRes.png'])
% print(f4H,'-dpng','-r100',[outPlotDTInsEng,'_lowRes.png'])
% print(f5H,'-dpng','-r100',[outPlotDTEng,'_lowRes.png'])
% print(f6H,'-dpng','-r100',[outPlotBlandAlt,'_lowRes.png'])
% 
% print(f1H,'-dpng','-r300',[outPlotInsEng,'_highRes.png'])
% print(f2H,'-dpng','-r300',[outPlotDT50Eng,'_highRes.png'])
% print(f3H,'-dpng','-r300',[outPlotDeltaEng,'_highRes.png'])
% print(f4H,'-dpng','-r300',[outPlotDTInsEng,'_highRes.png'])
% print(f5H,'-dpng','-r300',[outPlotDTEng,'_highRes.png'])
% print(f6H,'-dpng','-r300',[outPlotBlandAlt,'_highRes.png'])
% 
% saveas(f1H,outPlotInsEng,'fig');
% saveas(f2H,outPlotDT50Eng,'fig');
% saveas(f3H,outPlotDeltaEng,'fig');
% saveas(f4H,outPlotDTInsEng,'fig');
% saveas(f5H,outPlotDTEng,'fig');
% saveas(f6H,outPlotBlandAlt,'fig');

