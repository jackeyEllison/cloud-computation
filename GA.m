function [bestEverSolution, objList, energyList, penaltyList, stopIteration, timeDuration] = GA(cfg)
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
    
    %objList= best(:,1), energy 
    best = zeros(totalIterations, 3); 
    %Assume the first individual is the best initially
    bestSolution = Population(1, :);
    bestScore = Population(1, numOfDimension + 1); 
    iterationIndex = 1;
    startTime = tic;

    while iterationIndex <= totalIterations
        % Selection (tournament selection)
        selectedPopulation = zeros(size(Population));
        for i = 1:numOfParticles
            idx1 = randi(numOfParticles);
            idx2 = randi(numOfParticles);
            if Population(idx1, numOfDimension + 1) < Population(idx2, numOfDimension + 1)
                selectedPopulation(i, :) = Population(idx1, :);
            else
                selectedPopulation(i, :) = Population(idx2, :);
            end
        end
        % Crossover (single-point crossover)
        newPopulation = zeros(size(Population));
        for i = 1:2:numOfParticles
            parent1 = selectedPopulation(i, 1:numOfDimension);
            parent2 = selectedPopulation(i+1, 1:numOfDimension);
            crossoverPoint = randi([1, numOfDimension-1]);
            newIndividual1 = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
            newIndividual2 = [parent2(1:crossoverPoint), parent1(crossoverPoint+1:end)];
            newPopulation(i, 1:numOfDimension) = newIndividual1;
            newPopulation(i+1, 1:numOfDimension) = newIndividual2;
        end
        % Mutation (random mutation)
        mutationRate = 0.01;
        for i = 1:numOfParticles
            for d = 1:numOfDimension
                if rand < mutationRate
                    newPopulation(i, d) = rand; % Randomly mutate within the allowable range
                end
            end
        end
        % Evaluate new population
        for i = 1:numOfParticles
            [Pm, X, mu, FL, FE, lambda] = ...
                vector2matrix(M, K, J, newPopulation(i, 1:numOfDimension));
            [Pm, X, mu, FL, FE, lambda] = ...
                checkBound(M, K, J, Pm, X, mu, FL, FE, lambda, PmMax, FLMax, FEMax);
            newPopulation(i, 1:numOfDimension) = ...
                matrix2vextor(M, K, J, Pm, X, mu, FL, FE, lambda);
            [newObj, newEnergy, newPenalty] = ...
                objFunction(cfg, Pm, X, mu, FL, FE, lambda);
            newPopulation(i, numOfDimension + 1) = newObj;
            newPopulation(i, numOfDimension + 2) = newEnergy;
            newPopulation(i, numOfDimension + 3) = newPenalty;
            % Update the best solution if necessary
            if newObj < bestScore
                bestScore = newObj;
                bestSolution = newPopulation(i, :);
            end
        end
        % Update population
        Population = newPopulation;
        [newObj, newEnergy, newPenalty] = ...
                objFunction(cfg, Pm, X, mu, FL, FE, lambda);
        % Record best score
        best(iterationIndex, :) = [newObj, newEnergy, newPenalty];
        iterationIndex = iterationIndex + 1; 
    end

    bestEverSolution = bestSolution(1:numOfDimension);
    objList = best(:, 1);
    energyList = best(:, 2);
    penaltyList = best(:, 3);

    stopIteration = iterationIndex;
    timeDuration = toc(startTime);
    save("GAobj", "EnergyList")
    fclose(fid);
end