function parameters = runPertWC(islocal, replicates, genesTusRxns, parameterTypes, parameterVals, userData)
%
% See section 2.4.2 of competition page for how perturbations were generated
% based on gold (mutant) data. Unlike there, this script applies all pertubation
% arguments at once, in case it is desired to combine multiple perturbations, 
% which may be useful for later evaluating our approximation to the gold parameters.
%
% For a "WT" simulation, just pass {} for genesTusRxns, parameterTypes, 
% and parameterVals
%
% islocal: if not islocal, it will run on bitmill (TODO)
% genesTusRxns??:  cell array? of parameters to change
% parameterTypes:  e.g. 'PromAffinity', 'HalfLife', 'RxnKcat'
% parameterVals:  multiplicative perturbations to parameters
%
% TODO:
%
% Test that parfor works
%
%
% if using gurobi, try to specify solver algorithm to use
% 
% Add switch cases for other parameterTypes? Probably not necessary in 
% this *whole cell* script.
%


% Adpated from the user's guide, 7/29/2013
import edu.stanford.covert.cell.sim.util.*;
%import edu.stanford.covert.cell.sim.util.SimulationDiskUtil.*;
%import edu.stanford.covert.cell.sim.util.SummaryLogger.*;

wholecell_root = '/home/brandon/DREAM8/WholeCell';
output_root = '/home/brandon/DREAM8/WholeCell/output';
rundir = pwd();

%try % try running whole cell simulation

cd(wholecell_root)
% Select
%(1) simulation batch output directory and
%(2) simulation output directory
simBatch = datestr (now , 'yyyy_mm_dd_HH_MM_SS');
simBatchDir = [ output_root filesep simBatch ]
% create simulation batch directories
if ~isdir(simBatchDir)
    mkdir(simBatchDir); % create simulation batch output directory
end

[sim, kbWID] = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();


opts = sim.getOptions();
%check that gurobi is a mex function
%gurobi currently seems broken in the whole cell model, but not in the COBRA/metabolic model
%if exist('gurobi') == 3 
%  opts.processes.Metabolism.linearProgrammingOptions.solver = 'gurobi';
%end
sim.applyOptions(opts);
sim.applyOptions('lengthSec', 65000)

rnaPolTuBindingProbs = sim.getRnaPolTuBindingProbs();
rnaHalfLives = sim.getRnaHalfLives();
rxnKinetics = sim.getMetabolicReactionKinetics();

%Assume the lengths of genesTusRxns, parameterTypes, and parameterVals are equal
for i = 1:size(genesTusRxns, 1)
  switch parameterTypes{i}
    case 'PromAffinity'
      tuId = genesTusRxns{i, 2};
      sim.applyRnaPolTuBindingProbs(struct(tuId, parameterVals{i, 2} * rnaPolTuBindingProbs.(tuId)));
    case 'HalfLife'
      tuId = genesTusRxns{i, 2};
      sim.applyRnaHalfLives(struct(tuId, parameterVals{i, 2} * rnaHalfLives.(tuId)));
    case 'RxnKcat'
      rxnId = genesTusRxns{i, 3};
      sim.applyMetabolicReactionKinetics(struct(rxnId, struct('for', parameterVals{i, 2} * rxnKinetics.(rxnId).for)));
  end
end

parameters = sim.getAllParameters();

%perform simulations
parfor i = 1:replicates

simDir = [ output_root filesep simBatch filesep ...
    num2str(i)];
if ~isdir (simDir)
    mkdir (simDir); % create simulation output directory
end
if size(genesTusRxns,1) > 0
  paramFileName = fullfile(simDir, sprintf('parameters_%s__%s__%s.mat', ...
    cellstrjoin(genesTusRxns,'_'), cellstrjoin(parameterTypes,'_'), ... 
    cellstrjoin(parameterVals(:,1), '_')));
else
  paramFileName = fullfile(simDir, 'parameters_WT.mat');
end

parSaveFcn(paramFileName, 'parameters');

% setup loggers
summaryLogger = SummaryLogger (1, 1); % print to command line
summaryLogger.setOptions ( struct ('outputDirectory', simDir )); % save to disk

diskLogger = DiskLogger (simDir , 10); % save complete dynamics to disk
diskLogger.addMetadata(...
    'shortDescription', 'WT whole cell simulation', ...
    'longDescription', ['logging of WT whole cell simulation: ' paramFileName], ...
    'email', userData.email, ...
    'firstName', userData.firstName, ...
    'lastName', userData.lastName, ...
    'affiliation', userData.affiliation, ...
    'knowledgeBaseWID', kbWID, ...
    'revision', 1, ...
    'differencesFromRevision', [], ...
    'userName', userData.userName, ...
    'hostName', userData.hostName, ...
    'ipAddress', userData.ip);

loggers = {summaryLogger; diskLogger};
sim.run(loggers);
end % end parfor

simBatchDir

%catch err %
%end       % end of global try 

cd(rundir);
end % end of runPertWC

function parSaveFcn(fname, parameters)
  save( fname, 'parameters' );
end