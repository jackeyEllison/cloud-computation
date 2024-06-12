clear ; close all; clc
cfg = config_for_algo(); 
if cfg.fileName == "LSAG"
    [results,objList,energyList,penaltyList,...
        stopIteration,timeDuration] = LSAG(cfg);
end
 objPlot(objList,energyList,penaltyList,cfg.totalIterations,...
     cfg.M,cfg.K,cfg.J,cfg.fileName);





