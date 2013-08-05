% import classes
import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil ;
import edu.stanford.covert.cell.sim.util.SimulationEnsemble ;
% load simulaton object
sim = CachedSimulationObjectUtil.load();
comp = sim.compartment;

met = sim.process('Metabolism');
pc = sim.state('ProteinComplex');
pm = sim.state('ProteinMonomer');
rna = sim.state('Rna');

% load data
stateNames = {
    'Time'              'values'
    'Mass'              'cell'
    'MetabolicReaction' 'fluxs'
    'ProteinComplex'    'counts'
    'ProteinMonomer'    'counts'
    'Rna'               'counts'
};
states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, ...
    'extract', simIdx );

% specify target
% reactionTitle = 'AckA';
% fluxIdx = met.reactionIndexs (reactionTitle);
% cpxIdx = pc.matureIndexs ( pc.getIndexs ('MG_357_DIMER'));
% monIdx = pm.matureIndexs ( pm.getIndexs ('MG_357_MONOMER'));
% rnaIdx = rna.matureIndexs ( rna.getIndexs ('TU_260'));

% Isolate genes of interest
genesTusRxnsCpxs = {
    'MG_006' 'TU_003' 'Tmk' 'MG_006_DIMER'
    'MG_023' 'TU_011' 'Fba' 'MG_023_DIMER'
    'MG_047' 'TU_027' 'MetK' 'MG_047_TETRAMER'
    'MG_111' 'TU_069' 'Pgi' 'MG_111_DIMER'
    'MG_272' 'TU_180' 'AceE' 'MG_271_272_273_274_192MER'
    'MG_299' 'TU_203' 'Pta' 'MG_299_DIMER'
   'MG_330' 'TU_233' 'CmkA2' 'MG_330_MONOMER'
    'MG_357' 'TU_260' 'AckA' 'MG_357_DIMER'
    'MG_407' 'TU_294' 'Eno' 'MG_407_DIMER'
    'MG_431' 'TU_307' 'TpiA' 'MG_431_DIMER'
    };

for geneIdx = 1:10;
    reactionTitle = genesTusRxnsCpxs{geneIdx, 3};
    cpxTitle = genesTusRxnsCpxs{geneIdx, 4};
    monTitle = [genesTusRxnsCpxs{geneIdx, 1} '_MONOMER'];
    tuTitle = genesTusRxnsCpxs{geneIdx, 2};

    % % specify target
    fluxIdx = met.reactionIndexs (reactionTitle);
    if pc.getIndexs(cpxTitle)
        cpxIdx = pc.matureIndexs ( pc.getIndexs (cpxTitle));
    else 
        cpxIdx = -1;
    end
    monIdx = pm.matureIndexs ( pm.getIndexs (monTitle));
    rnaIdx = rna.matureIndexs ( rna.getIndexs (tuTitle));

    % plot
    handle=figure();
    subplot(5, 1, 1);
    plot(permute(states.Time.values, [1 3 2]), permute(sum(states.Mass.cell, 2), ...
            [1 3 2]) * 1e15 );
    title_object = title([reactionTitle ':' monTitle ':' cpxTitle ':' tuTitle]);
    set(title_object, 'interpreter', 'none')
    ylabel('Mass (fg)');
    subplot(5, 1, 2);
    plot(permute(states.Time.values, [1 3 2]) , ...
    permute(states.MetabolicReaction.fluxs(fluxIdx, :, :) , [1 3 2]) * 1e-3);
    ylabel ({ 'Flux ' '(10^3 rxn s^{ -1}) '});

    subplot (5, 1, 3);
    if ~(cpxIdx==-1)
        plot(permute(states.Time.values, [1 3 2]) , ...
        permute(states.ProteinComplex.counts(cpxIdx , comp.cytosolIndexs , :) , [1 ...
            3 2]));
        ylabel ('Complex ');
    end 
    subplot (5, 1, 4);
    plot(permute(states.Time.values, [1 3 2]) , ...
    permute(states.ProteinMonomer.counts(monIdx , comp.cytosolIndexs , :) , [1 ...
        3 2]));
    ylabel('Monomer ');
    subplot(5, 1, 5);
    plot(permute(states.Time.values, [1 3 2]) , permute(states.Rna.counts(rnaIdx , ...
        comp.cytosolIndexs , :) , [1 3 2]));
    ylabel ('RNA ');
    xlabel ('Time (s)');
    
    saveas(handle, [simBatchDir filesep sprintf('roi_%d.png',geneIdx)])
end