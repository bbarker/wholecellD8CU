%set warnings and MATLAB path
%setWarnings();
%setPath();
 
rundir = pwd();
cd /home/brandon/DREAM8/WholeCell   
%import classes
import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
    
%load simulation object
sim = CachedSimulationObjectUtil.load();
    
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
    
%optionally, sample initial conditions
sim.initializeState();
    
%simulate dynamics for 100s
lengthSec = 100;
for i = 1:lengthSec
    met.evolveState();
end

cd rundir