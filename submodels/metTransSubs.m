%import classes
import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
    
%load simulation object
sim = CachedSimulationObjectUtil.load();
    
% get state and sub-model handles
time = sim.state('Time');
met = sim.process('Metabolism');
transcription = sim.process('Transcription');
        
%simulate
lengthSec = 100;
for i = 1:lengthSec
    time.values = i;
    
    met.copyFromState();
    met.evolveState();
    met.copyToState();
    
    transcription.copyFromState();
    transcription.evolveState();
    transcription.copyToState();
end