%
% Run the metabolic submodel
% Save growth rate and fluxes for inspection.
%

rundir = pwd();
cd(WCDIR) 

%import classes
import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
    
%load simulation object
sim = CachedSimulationObjectUtil.load();

opts = sim.getOptions();
%check that gurobi is a mex function
%gurobi currently seems broken in the whole cell model, but not in the COBRA/metabolic model
%if exist('gurobi') == 3 
%  opts.processes.Metabolism.linearProgrammingOptions.solver = 'gurobi';
%end
sim.applyOptions(opts);
    
%optionally, set simulation options
sim.applyOptions('seed', 1);
    
%optionally, set simulation parameter
sim.applyMetabolicReactionKinetics(struct(...
    'AtpA', struct(...
        'for', 1, ...
        'rev', -1 ...
        ) ...
    ));
    
%get handle to metabolism sub-model
met = sim.process('Metabolism');
mr = sim.state('MetabolicReaction');
    
%optionally, sample initial conditions
sim.initializeState();
    
%simulate dynamics for 100s
lengthSec = 100;
growth = zeros(lengthSec, 1);
fluxs = zeros(lengthSec, length(mr.fluxs));
for i = 1:lengthSec
    met.evolveState();
    growth(i) = mr.growth;
    fluxs(i,:)  = mr.fluxs;
end

cd(rundir);