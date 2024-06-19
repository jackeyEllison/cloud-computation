clear ; close all; clc
cfg = config_for_algo(); 
if cfg.fileName == "LSAG"
    [results,objList,energyList,penaltyList,...
        stopIteration,timeDuration] = LSAG(cfg);
 elseif cfg.fileName == "SA"
     [results,objList,energyList,penaltyList,...
         stopIteration,timeDuration] = SA(cfg);
 elseif cfg.fileName == "GA3"
    [results,objList,energyList,penaltyList,...
         stopIteration,timeDuration] = GA3(cfg);
 elseif cfg.fileName == "GWO"
      [results,objList,energyList,penaltyList,...
         stopIteration,timeDuration] = GWO(cfg);
end
   objPlot(objList,energyList,penaltyList,cfg.totalIterations,...
       cfg.M,cfg.K,cfg.J,cfg.fileName);
   




