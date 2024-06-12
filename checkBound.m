
% 每次计算目标函数之前都需要对X mu lambda变量进行边界检查，否则目标函数会报错
% 所有决策变量的检查下界的时候不能直接为0，会导致后续变量(e.g. Rup Rdown...)或者在计算最终object的时候出bug
%   解决办法：如果越过下界，使决策变量的赋值能够让最后计算的penalty或者energy趋向于正无穷Inf或者一个比较大的正数
% checkPosBound函数检查位移边界
% nargin  : number of function input arguments
% varargin: length input arguments list ( 1*n cell )

function [Pm_new,X_new,mu_new,FL_new,FE_new,lambda_new] = ...
    checkBound(M,K,J,Pm,X,mu,FL,FE,lambda,PmMax,FLMax,FEMax)
    %% 2. 进行相关处理
    %% 2.1 check lower bound:
    % 决策变量X不需要检查下界，因为只选出最大值为1，其余均为0：
    X_new = zeros(M,J);
    for m = 1:M    
       [~,j] = max(X(m,:));
        X_new(m,j) = 1;
    end

    % 决策变量mu检查是否某一行全为0，如果是0则令成比较小的数
    % 当X_new(m,j)==1 mu(m,j)==0时，这样会导致mu(m,:)这一行都为0，
    % 在后续做稠密矩阵转换成稀疏矩阵的时候会报错(e.g. Rup Rdown...)
    % 解决办法就是当X_new(m,j)==1 mu(m,j)==0的时候，令mu(m,j)为一个很小的数，能够让最后计算的energy趋向于正无穷Inf
    % 当然如果在后面放缩mu的时候发现只有用户m与基站j连接，那么令mu(m,j)=1。
    for m = 1:M
        for j = 1:J
            if mu(m,j)==0
                if X_new(m,j)==1
                    mu(m,j) = 1*10^(-3);
                elseif X_new(m,j)==0
                    mu(m,j) = 0;
                end
            end
        end
    end

    % 决策变量FL和FE不能直接为0，因为在目标函数的计算中，他们是作为分母
    % 如果计算的时候分子也同时为0，即0/0=NaN，会有严重bug出现！！！！（查阅NaN相关资料）
    % 所以这里令其得1*10^(-10)，能够让最后计算的penalty趋向于正无穷Inf
    FL(FL<0) = 1*10^(-10);
    FE(FE<0) = 1*10^(-10);

    % 决策变量lambda为0，有可能会出现sum(lambda,3)为0的特殊情况，这样任务没有充分完成
    % 令其得1*10^(-10)，如果出现上面的特殊情况，交给后面去放缩（都为0.3333）
    lambda(lambda<0) = 1*10^(-10); 




    %% 2.2 check upper bound:
    % 找到Pm中大于上界的，令其为PmMax
    Pm_new = Pm;
    Pm_new(Pm>PmMax) = PmMax;

    % X_new: 取X中每一行（每一个用户）的最大值为1
    % update according to row m, each row is guaranteed to only has one '1' 
    % j: index of maximum value of each user m    
    % 写到了上面 lower bound里

    % mu_new: 按照X_new来重新按比例分配带宽
    mu_new = zeros(M,J);    
    % mu 2-D (M*J):
    for j = 1:J    % update according to column j
        % num: number of users asscociated with SBS j
        num = sum( X_new(:,j) );
        % 0. if there is no user accociated with SBS j 
        if (num == 0)
            % do nothing 

        % 1. if SBS j only associated with 1 user m, set mu_mj = 1
        elseif (num == 1)
            mu_new( X_new(:,j)==1, j ) = 1;

        % 2. SBS j associated with more than 1 user m, update mu_mj with
        % proportion in mu
        else
            % denominator: 在新的X_new分配策略下，连接第j个基站的全部用户带宽之和
            % m: m是基站j所连接的全部用户 可能是scalar 也可能是vector
            m = find( X_new(:,j)==1 );
            denominator = sum( mu(X_new(:,j)==1,j) );
            mu_new(m,j) = mu(m,j)./denominator;
        end        
    end

    % FL_new: 找到移动端超过CPU上界的用户，令其按比例缩小
    % m: 超过移动端用户CPU上界的用户
    % sum(FL(m,:),2): 超过移动端用户CPU上界的用户的CPU速度之和，作为分母进行等比例缩小
    FL_new = FL;
    m = sum(FL,2)>FLMax;
    FL_new(m,:) = FLMax .* FL(m,:) ./ sum(FL(m,:),2);

    % FE_new: 找到边缘端超过其CPU上界的基站，令其按比例缩小
    % j: 超过边缘端用户CPU上界的基站
    % sum(sum(FE(:,:,j))): 超过边缘端用户CPU上界的基站的CPU速度之和，作为分母进行等比例缩小
    FE_new = FE;
    j = sum(sum(FE)) > FEMax;
    FE_new(:,:,j) = FEMax .* FE(:,:,j) ./ sum(sum(FE(:,:,j)));

    % lambda_new: 重新按比例进行offloading 
    lambda_new = zeros(M,K,3);
    denominator = sum(lambda,3); % size(M,K)
    denominator(denominator==0) = 1; % 防止分母为0
    lambda_new = lambda_new + lambda./denominator;

end


