function [bestEverSolution, objList, energyList, penaltyList, stopIteration, timeDuration] = SA(cfg)
    % Initialize outputs
    bestEverSolution = [];
    objList = [];
    energyList = [];
    penaltyList = [];
    stopIteration = -1;
    timeDuration = 0;

    obj = [];
    pen = [];
    M = cfg.M;
    J = cfg.J;
    K = cfg.K;

    FLMax = cfg.upbound.FLMax;
    FEMax = cfg.upbound.FEMax;
    PmMax = cfg.upbound.PmMax;

    numOfDimension = M*1 + M*J + M*J + M*K + M*K*J + M*K*3; 
    numOfParticles = cfg.numOfParticles;         
    totalIterations = cfg.totalIterations;       

    %粒子群的位置信息，数目， 属性+（1， 2，3）
    Population = zeros(numOfParticles, numOfDimension+3);

    fid = fopen('result.txt', 'wt');
   
    CPU = cfg.user.alpha;
    Memory = cfg.user.g;
    Csum = sum(CPU, 2);
    Msum = sum(Memory, 2);
    TotalR = [Csum Msum];

    % 均值-方差标准化
    meanTotalR = mean(TotalR);
    stdTotalR = std(TotalR);
    nor = (TotalR - meanTotalR) ./ stdTotalR;
    [SortTotal, index] = sort(nor);

    CPUSBSmax = [3*10^(10) 1*10^(10) 1.5*10^(10)];
    MemSBSmax = [2*1024*1024*1024*8 2.5*1024*1024*1024*8 1.5*1024*1024*1024*8];
    CPUSBSmax = CPUSBSmax';
    MemSBSmax = MemSBSmax';
    TotalRSBS = [CPUSBSmax MemSBSmax];
    [row1, column1] = size(TotalRSBS);
    norSBS = [];
    for i = 1:row1
        norSBS = [norSBS; norm(TotalRSBS(i,:))];
    end
    [SortTotalSBS, indexSBS] = sort(norSBS);

     % 初始化矩阵X1
      X1 = [0, 0, 1; 1, 0, 0; 0, 1, 0; 0, 1, 0; 0, 0, 1; 0, 0, 1];

      % 初始化X7, X8, X9, X10并进行赋值
      X_temp = zeros(4, 3);
      for i = 1:3
          randomNumber = rand();
          if  randomNumber < 0.5 
               X_temp(i, indexSBS(1)) = 1;
          elseif  0.5 <= randomNumber && randomNumber < 0.75
               X_temp(i, indexSBS(2)) = 1;
          else
              X_temp(i, indexSBS(3)) = 1;
          end         
      end
      X_temp(4, :) = [1, 0, 0];

      % 将X_temp中的行添加到X1
      X1 = [X1; X_temp];

    % Initialize population
    for i = 1:numOfParticles
        [Pm, X, mu, FL, FE, lambda] = ...
            vector2matrix(M, K, J, Population(i, 1:numOfDimension));
        [Pm, X, mu, FL, FE, lambda] = ...
            randPosition(M, K, J, Pm, X, mu, FL, FE, lambda, PmMax, FLMax, FEMax);
        X = X1;
        Population(i, 1:numOfDimension) = ...
            matrix2vextor(M, K, J, Pm, X, mu, FL, FE, lambda);
    end

    % Evaluate initial population
    for i = 1:numOfParticles
        [Pm, X, mu, FL, FE, lambda] = ...
            vector2matrix(M, K, J, Population(i, 1:numOfDimension));
        [oldObj, oldEnergy, oldPenalty] = ...
            objFunction(cfg, Pm, X, mu, FL, FE, lambda);
        Population(i, numOfDimension + 1) = oldObj;
        Population(i, numOfDimension + 2) = oldEnergy;
        Population(i, numOfDimension + 3) = oldPenalty;
    end

    best = zeros(totalIterations, 3); %objList= best(:,1), energy 

    bestSolution = Population(1, :);
    bestScore = Population(1, numOfDimension + 1); % Assume the first individual is the best initially

    iterationIndex = 1;
    startTime = tic;



    % Initialize current solution
    currentSolution = Population(1, :);
    currentScore = currentSolution(numOfDimension + 1);
    % Initialize temperature
    initialTemperature = 100;
    finalTemperature = 1e-3;
    coolingRate = 0.95;
    currentTemperature = initialTemperature;
    
    while iterationIndex <= totalIterations
    % Simulated Annealing main loop
    % Generate a new solution by perturbing the current solution
        for i = 1:numOfParticles
            newSolution = Population(i, :);
            perturbationIndex = randi(numOfDimension); % Randomly choosing a dimension to perturb
            newSolution(perturbationIndex) = rand; % Randomly mutate the chosen dimension

            % Ensure the new solution is within bounds
            [Pm, X, mu, FL, FE, lambda] = vector2matrix(M, K, J, newSolution(1:numOfDimension));
            [Pm, X, mu, FL, FE, lambda] = checkBound(M, K, J, Pm, X, mu, FL, FE, lambda, PmMax, FLMax, FEMax);
            newSolution(1:numOfDimension) = matrix2vextor(M, K, J, Pm, X, mu, FL, FE, lambda);

            % Evaluate the new solution
            [newObj, newEnergy, newPenalty] = objFunction(cfg, Pm, X, mu, FL, FE, lambda);
            newSolution(numOfDimension + 1) = newObj;
            newSolution(numOfDimension + 2) = newEnergy;
            newSolution(numOfDimension + 3) = newPenalty;

            % Acceptance criterion
            if newObj < currentScore
                Population(i, :) = newSolution;
                currentScore = newObj;
            else
                delta = newObj - currentScore;
                acceptanceProbability = exp(-delta / currentTemperature);
                if rand < acceptanceProbability
                    Population(i, :) = newSolution;
                    currentScore = newObj;
                end
            end

            % Update the best solution if necessary
            if newObj < bestScore
                bestScore = newObj;
                bestSolution = newSolution;
            end
        end
        % Decrease temperature
        currentTemperature = currentTemperature * coolingRate;
        best(iterationIndex, :) = [bestScore, bestSolution(numOfDimension + 2), bestSolution(numOfDimension + 3)];
       % Record the best score
        iterationIndex = iterationIndex + 1;
    end
    bestEverSolution = bestSolution(1:numOfDimension);
    objList = best(:, 1);
    energyList = best(:, 2);
    penaltyList = best(:, 3);
    
    stopIteration = iterationIndex;
    timeDuration = toc(startTime);
    save("SAobj", "energyList")
    fclose(fid);
end

