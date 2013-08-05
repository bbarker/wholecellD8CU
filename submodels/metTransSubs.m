% 
% Run a (presumably largely unlinked, other than nucleotide demands from
% transcription) metabolic + transcriptiona submodel. Growth and fluxes
% are saved for inspection.
%
%
rundir = pwd();
cd(WCDIR) 

%import classes
import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
    
%load simulation object
sim = CachedSimulationObjectUtil.load();
    
% get state and sub-model handles
time = sim.state('Time');
met = sim.process('Metabolism');
transcription = sim.process('Transcription');

mr = sim.state('MetabolicReaction');
        
%simulate
lengthSec = 100;
growth = zeros(lengthSec, 1);
growth = zeros(lengthSec, length(mr.fluxs));

for i = 1:lengthSec
    time.values = i;
    
    met.copyFromState();
    met.evolveState();
    met.copyToState();

    growth(i) = mr.growth;
    fluxs(i,:)  = mr.fluxs;

    transcription.copyFromState();
    transcription.evolveState();
    transcription.copyToState();
end

cd(rundir);