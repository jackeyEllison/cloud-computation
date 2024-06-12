function [bestEverSolution, objList, energyList, penaltyList, stopIteration, timeDuration] = LSAG(cfg)

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
    [row, ~] = size(TotalR);
    nor = [];
    for i = 1:row
        nor = [nor; norm(TotalR(i,:))];
    end
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
    
    
    
    % 改进点1
    X1 = [1, 0, 0; 1, 0, 0; 0, 1, 0; 0, 1, 0; 0, 0, 1; 0, 0, 1];

    X7 = [0, 0, 0];
    X8 = [0, 0, 0];
    X9 = [0, 0, 0];
    X10 = [0, 0, 0];

    X7(indexSBS(1)) = 1;
    X8(indexSBS(2)) = 1;
    X9(indexSBS(3)) = 1;
    X10 = [0, 0, 1];

    X1 = [X1; X7];
    X1 = [X1; X8];
    X1 = [X1; X9];
    X1 = [X1; X10];

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
            [pm, X, mu, FL, FE, lambda] = ...
                vector2matrix(M, K, J, Position(i, 1:numOfDimension));
            [Pm, X, mu, FL, FE, lambda] = ...
                checkBound(M, K, J, Pm, X, mu, FL, FE, lambda, PmMax, FLMax, FEMax);
            Position(i, 1:numOfDimension) = ...
                matrix2vextor(M, K, J, Pm, X, mu, FL, FE, lambda);
            [oldObj, oldEnergy, oldPenalty] = ...
                objFunction(cfg, pm, X, mu, FL, FE, lambda);

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

        sigma = 0.6; 
        if iterationIndex/totalIterations < sigma
            a = -0.01 * iterationIndex / totalIterations + 0.001;
        else
            a = 2 - 0.01 * iterationIndex / totalIterations + 0.001;
        end
        
        

        for i = 1:numOfParticles
            for d = 1:numOfDimension
                
                
                r1 = rand();
                r2 = rand();

                %当前位置和当前最优位置之间的距离用D_alpha表示
                A1 = 2 * a * r1 - a;
                C1 = 2 * r2;
                D_alpha = abs(C1 * Alpha_pos(d) - Position(i, d));

                
                
                r1 = rand();
                r2 = rand();
                %当前位置和历史最优位置之间的距离用D_beta表示
                A2 = 2 * a * r1 - a;
                C2 = 2 * r2;
                D_beta = abs(C2 * Beta_pos(d) - Position(i, d));
        
                %将当前位置信息存起来放到Position_old中
                Position_old = Position(i, :);
                A = rand(); %随机生成一个随机数，如果小于0.5按照公式更新粒子位置， 如果大于0.5按照levy飞行更新粒子位置。
                if abs(A) < 0.5
                    Position(i, d) = 0.5 * (Alpha_pos(d) - A1 * D_alpha + Beta_pos(d) - A2 * D_beta);
                else
                    beta = 1.5;
                    sigma_u = ((gamma(1+beta)*sin(pi*beta/2))/(gamma((1+beta)/2)*beta*2^(0.5*(beta-1))))^(1/beta);
                    u = normrnd(0, sigma_u);
                    v = normrnd(0, 1);
                    alpha_levi = (u / abs(v)^(1/beta) * (Position(i, d) - Alpha_pos(d)) + u / abs(v)^(1/beta) * (Position(i, d) - Beta_pos(d))) / 2;
                    Position(i, d) = 0.5 * (Alpha_pos(d) - A1 * D_alpha + Beta_pos(d) - A2 * D_beta) + alpha_levi; 
                end
            end
        end

        rnew = rand();
        p = rand();
        [pm, X, mu, FL, FE, lambda] = ...
            vector2matrix(M, K, J, Position(i, 1:numOfDimension));
        [Pm, X, mu, FL, FE, lambda] = ...
            checkBound(M, K, J, Pm, X, mu, FL, FE, lambda, PmMax, FLMax, FEMax);
        Position(i, 1:numOfDimension) = ...
            matrix2vextor(M, K, J, Pm, X, mu, FL, FE, lambda); 
        [newObj, oldEnergy, oldPenalty] = ...
            objFunction(cfg, pm, X, mu, FL, FE, lambda);

        [pm, X, mu, FL, FE, lambda] = ...
            vector2matrix(M, K, J, Position_old(1:numOfDimension));
        [Pm, X, mu, FL, FE, lambda] = ...
            checkBound(M, K, J, Pm, X, mu, FL, FE, lambda, PmMax, FLMax, FEMax);
        Position(i, 1:numOfDimension) = ...
            matrix2vextor(M, K, J, Pm, X, mu, FL, FE, lambda); 
        [oldObj, oldEnergy, oldPenalty] = ...
            objFunction(cfg, pm, X, mu, FL, FE, lambda);

        temperature = 100;
        delta = oldObj - newObj; 
        if delta > 0
            % 现在的好
        else
            probability = exp(-delta / temperature);
            randomNumber = rand();
            if probability > randomNumber
                Position(i, :) = Position_old;
            end 
        end

        best(iterationIndex, :) = [Alpha_score, Beta_score, Delta_score]; 
        
        fprintf(fid, 'iterationIndex = %d | ', iterationIndex);
        fprintf(fid, 'Alpha_score = %f\n', Alpha_score);
        iterationIndex = iterationIndex + 1;
        obj = [obj Alpha_score];
        pen = [pen, oldPenalty];
        save("LSAGobj", "obj")
    end

    % Assign output values
    bestEverSolution = Alpha_pos; % Assuming Alpha_pos is the best solution
    objList = best(:, 1);
    energyList = best(:, 2);
    penaltyList = best(:, 3);
    timeDuration = toc(startTime);

    fclose(fid);
end
 