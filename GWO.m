function [bestEverSolution, objList, energyList, penaltyList, stopIteration, timeDuration] = GWO(cfg)
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
    Position = zeros(numOfParticles, numOfDimension+3);
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
%   或最大最小标准化
%   minTotalR = min(TotalR);
%   maxTotalR = max(TotalR);
%   nor = (TotalR - minTotalR) ./ (maxTotalR - minTotalR);
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
    X1 = [1, 0, 0; 1, 0, 0; 0, 1, 0; 0, 1, 0; 0, 0, 1; 0, 0, 1];
    % 初始化X7, X8, X9, X10并进行赋值
    X_temp = zeros(4, 3);
        for i = 1:3
            X_temp(i, indexSBS(i)) = 1;
        end
    X_temp(4, :) = [0, 0, 1];
    % 将X_temp中的行添加到X1
    X1 = [X1; X_temp];
    
    
    %首先对粒子群的所有粒子的属性向量转成矩阵便于计算，然后进行随机位置， 最后将属性矩阵变成属性向量
    for i = 1:numOfParticles
        [Pm, X, mu, FL, FE, lambda] = ...
            vector2matrix(M, K, J, Position(i, 1:numOfDimension));
        [Pm, X, mu, FL, FE, lambda] = ...
            randPosition(M, K, J, Pm, X, mu, FL, FE, lambda, PmMax, FLMax, FEMax);
        X = X1;
        Position(i, 1:numOfDimension) = ...
            matrix2vextor(M, K, J, Pm, X, mu, FL, FE, lambda);
    end
    
    
    %首先将属性向量变成属性矩阵，然后通过objFunction计算SDS, SBS， MBS的cpu, memory的时间占用
    %通过时间占用计算SDS, SBS, MBS的energy返回oldObj, oldEnergy, oldPenalty分别存在每个粒子的倒数
    %第三， 第二， 第一属性中。
    for i = 1:numOfParticles
        [Pm, X, mu, FL, FE, lambda] = ...
            vector2matrix(M, K, J, Position(i, 1:numOfDimension));
        [oldObj, oldEnergy, oldPenalty] = ...
            objFunction(cfg, Pm, X, mu, FL, FE, lambda);
        Position(i, numOfDimension + 1) = oldObj;
        Position(i, numOfDimension + 2) = oldEnergy;
        Position(i, numOfDimension + 3) = oldPenalty;
    end
    

    %创建了三个行向量，行向量的维度为numofDimension+3
    %Alpha_score, Beta_score, Delta_score三个变量的初始值都为正无穷大
    Alpha_pos = zeros(1, numOfDimension+3); 
    Alpha_score = inf; 
    
    Beta_pos = zeros(1, numOfDimension+3);
    Beta_score = inf; 
    
    Delta_pos = zeros(1, numOfDimension+3);
    Delta_score = inf; 
    
    best = zeros(totalIterations, 3); %objList= best(:,1), energyList = best(:, 2) penaltyList = best(:, 3)

    
    iterationIndex = 1;
    startTime = tic; 
    
    while iterationIndex <= totalIterations
        for i = 1:numOfParticles
            [Pm, X, mu, FL, FE, lambda] = ...
                vector2matrix(M, K, J, Position(i, 1:numOfDimension));
            [Pm, X, mu, FL, FE, lambda] = ...
                checkBound(M, K, J, Pm, X, mu, FL, FE, lambda, PmMax, FLMax, FEMax);
            Position(i, 1:numOfDimension) = ...
                matrix2vextor(M, K, J, Pm, X, mu, FL, FE, lambda);
            [oldObj, oldEnergy, oldPenalty] = ...
                objFunction(cfg, Pm, X, mu, FL, FE, lambda);
            if oldObj < Alpha_score
                Alpha_score = oldObj;
                Alpha_pos = Position(i, :);
            end
            if oldObj > Alpha_score && oldObj < Beta_score
                Beta_score = oldObj;
                Beta_pos = Position(i, :);
            end
            if oldObj > Alpha_score && oldObj > Beta_score && oldObj < Delta_score
                Delta_score = oldObj;
                Delta_pos = Position(i, :);
            end
        end

        
        a = 2 - 2 * iterationIndex / totalIterations;
        for i = 1:numOfParticles
            for d = 1:numOfDimension
                r1 = rand();
                r2 = rand();

                A1 = 2 * a * r1 - a;
                C1 = 2 * r2;
                D_alpha = abs(C1 * Alpha_pos(d) - Position(i, d));
                X1 = Alpha_pos(d) - A1 * D_alpha;

                r1 = rand();
                r2 = rand();

                A2 = 2 * a * r1 - a;
                C2 = 2 * r2;
                D_beta = abs(C2 * Beta_pos(d) - Position(i, d));
                X2 = Beta_pos(d) - A2 * D_beta;

                r1 = rand();
                r2 = rand();

                A3 = 2 * a * r1 - a;
                C3 = 2 * r2;
                D_delta = abs(C3 * Delta_pos(d) - Position(i, d));
                X3 = Delta_pos(d) - A3 * D_delta;

                Position(i, d) = (X1 + X2 + X3) / 3;
            end
        end

        % Update positions and bounds
        for i = 1:numOfParticles
            [pm, X, mu, FL, FE, lambda] = ...
                vector2matrix(M, K, J, Position(i, 1:numOfDimension));
            [Pm, X, mu, FL, FE, lambda] = ...
                checkBound(M, K, J, Pm, X, mu, FL, FE, lambda, PmMax, FLMax, FEMax);
            Position(i, 1:numOfDimension) = ...
                matrix2vextor(M, K, J, Pm, X, mu, FL, FE, lambda);
            [newObj, newEnergy, oldPenalty] = ...
                objFunction(cfg, pm, X, mu, FL, FE, lambda);
            Position(i, numOfDimension + 1) = newObj;
            Position(i, numOfDimension + 2) = newEnergy;
            Position(i, numOfDimension + 3) = oldPenalty;
        end

        best(iterationIndex, :) = [Alpha_score, Beta_score, Delta_score];

        fprintf(fid, 'iterationIndex = %d | ', iterationIndex);
        fprintf(fid, 'Alpha_score = %f\n', Alpha_score);
        iterationIndex = iterationIndex + 1;
        obj = [obj Alpha_score];
        pen = [pen, oldPenalty];
        save("GWOobj", "obj")
    end

    
    % Assign output values
    bestEverSolution = Alpha_pos; % Assuming Alpha_pos is the best solution
    objList = best(:, 1);
    energyList = best(:, 2);
    penaltyList = best(:, 3);
    timeDuration = toc(startTime);
    fclose(fid);
end
 