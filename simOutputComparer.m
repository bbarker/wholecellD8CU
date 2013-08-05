function [log1 states1]=simOutputComparer(file1,file2)
%
%  Authors:
%  Yiping Wang (August 2013)
%
%  extract flux information from the logger output, 
%  which we can use to compare differences between the full and surrogate submodels,
%  or between different simulations for the same model type. 
%
%
log1=load(file1);%'output/runSimulation/2013_08_04_17_25_01/1/summary.mat');
log2=load(file2);

%log1
outputstograph1={log1.mass*1e15,log1.ploidy,log1.rnas(1,:),log1.proteins(1,:)};
stringstograph1={'mass','ploidy','rna','protein'};
for i=1:length(outputstograph1)
    figure('Position',[300,300,1200,500],'Visible','off');
    a=plot(log1.time,outputstograph1{i});
    saveas(a,['../DREAMfigures/' stringstograph1{i} '.png']);
end
stateNames={'Time' 'values'
	    'Mass' 'cell'
	    'MetabolicReaction' 'fluxs'};
states1=edu.stanford.covert.cell.sim.util.SimulationEnsemble.load(file2,stateNames,[],[],1,'extract',1);
states2=edu.stanford.covert.cell.sim.util.SimulationEnsemble.load(file1,stateNames,[],[],1,'extract',1);
meanFluxdiff=mean(mean(abs(states1.MetabolicReaction.fluxs()-states2.MetabolicReaction.fluxs())));
fluxCorr=corr(states1.MetabolicReaction.fluxs(),states2.MetabollicReaction.fluxs());
end