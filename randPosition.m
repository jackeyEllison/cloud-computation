
% 随机初始化Position
% nargin  : number of function input arguments
% varargin: length input arguments list ( 1*n cell )

function [Pm,X,mu,FL,FE,lambda] = randPosition(varargin)

    if nargin ~= 6 && nargin ~= 12 && nargin ~= 9
        disp("!!!!Some mistake in randPosition.m, check the number of input variables.");        
        return
    
    elseif nargin == 6  % 如果输入变量个数为3，是FOM策略
        %% 1. 从varargin这个输入变量cell中取出变量
        M = varargin{1};
        J = varargin{2};
        Pm = varargin{3};
        X = varargin{4};
        mu = varargin{5};   
        PmMax = varargin{6};   
        %% 2. 进行相关处理
        FL = 0;FE = 0;lambda = 0;
        % Pm 1-D (M*1): 
        Pm(1:M) = rand(M,1) * PmMax;

        % X  2-D (M*J):
        for m = 1:M    
            % update according to row m, each row is guaranteed to only has one '1' 
            % randi(J):integer values drawn from the discrete uniform distribution on 1:J
            X(m,randi(J)) = 1; 
        end

        % mu 2-D (M*J):
        for j = 1:J    % update according to column j
            % num: number of users asscociated with SBS j
            num = sum( X(:,j) );
            % 0. if there is no user accociated with SBS j 
            if (num == 0)
                % do nothing 

            % 1. if SBS j only associated with 1 user m, set mu_mj = 1
            elseif (num == 1)
                mu( X(:,j)==1, j ) = 1;

            % 2. SBS j associated with more than 1 user m, random mu_mj
            else
                % 2.1 random mu in former M-1 rows
                mu( X(1:end-1,j)==1, j ) = 1/num .* rand( sum( X(1:end-1,j)==1 ), 1 );
                % 2.2 update mu in last row: the sum of each column is guaranteed to be 1
                % 2.2.1 如果最后一个用户为1，正常更新
                if( X(end,j)==1 )
                    mu( end, j ) = 1 - sum( mu(1:end-1,j) );
                % 2.2.2 如果最后一个用户为0，需要将前面等于1的用户重新初始化
                else % x(end,j)==0
                    % index: a vector indicates the position of x==1(which user m associated with SBS j)
                    index = find( X(1:end,j)==1 );
                    % among users who associated with SBS j, random former end-1
                    mu( index(1:end-1), j ) = 1/num .* rand( num-1, 1 );
                    mu( index(end), j ) = 1 - sum( mu(index(1:end-1),j) );
                end    
            end        
        end
        
    elseif nargin == 9     % 如果输入变量个数为9，是FOS策略
        %% 1. 从varargin这个输入变量cell中取出变量        
        M = varargin{1};
        K = varargin{2};
        J = varargin{3};
        Pm = varargin{4};
        X = varargin{5};
        mu = varargin{6}; 
        FE = varargin{7};
        PmMax = varargin{8};
        FEMax = varargin{9};
        
        %% 2. 进行相关处理
        % Pm 1-D (M*1): 
        Pm(1:M) = rand(M,1) * PmMax;

        % X  2-D (M*J):
        for m = 1:M    
            % update according to row m, each row is guaranteed to only has one '1' 
            % randi(J):integer values drawn from the discrete uniform distribution on 1:J
            X(m,randi(J)) = 1; 
        end

        % mu 2-D (M*J):
        for j = 1:J    % update according to column j
            % num: number of users asscociated with SBS j
            num = sum( X(:,j) );
            % 0. if there is no user accociated with SBS j 
            if (num == 0)
                % do nothing 

            % 1. if SBS j only associated with 1 user m, set mu_mj = 1
            elseif (num == 1)
                mu( X(:,j)==1, j ) = 1;

            % 2. SBS j associated with more than 1 user m, random mu_mj
            else
                % 2.1 random mu in former M-1 rows
                mu( X(1:end-1,j)==1, j ) = 1/num .* rand( sum( X(1:end-1,j)==1 ), 1 );
                % 2.2 update mu in last row: the sum of each column is guaranteed to be 1
                % 2.2.1 如果最后一个用户为1，正常更新
                if( X(end,j)==1 )
                    mu( end, j ) = 1 - sum( mu(1:end-1,j) );
                % 2.2.2 如果最后一个用户为0，需要将前面等于1的用户重新初始化
                else % x(end,j)==0
                    % index: a vector indicates the position of x==1(which user m associated with SBS j)
                    index = find( X(1:end,j)==1 );
                    % among users who associated with SBS j, random former end-1
                    mu( index(1:end-1), j ) = 1/num .* rand( num-1, 1 );
                    mu( index(end), j ) = 1 - sum( mu(index(1:end-1),j) );
                end    
            end        
        end
        
        % FE 3-D (M*K*J):
        FE(:,:,:) = (FEMax/(M*K)).*rand(M,K,J);
        FL = 0;lambda = 0;
        
        
    elseif nargin == 12     % 如果输入变量个数为12，是 PO策略 || 原始策略
        %% 1. 从varargin这个输入变量cell中取出变量
        M = varargin{1};
        K = varargin{2};
        J = varargin{3};
        Pm = zeros(M,1);
        X  = zeros(M,J);
        mu = zeros(M,J); 
        FL = zeros(M,K);
        FE = zeros(M,K,J);
        lambda = zeros(M,K,3);
        PmMax = varargin{10}; 
        FLMax = varargin{11};
        FEMax = varargin{12};
        
        %% 2. 进行相关处理
        % Pm 1-D (M*1): 
        Pm(1:M) = rand(M,1) * PmMax;

        % X  2-D (M*J):
        for m = 1:M    
            % update according to row m, each row is guaranteed to only has one '1' 
            % randi(J):integer values drawn from the discrete uniform distribution on 1:J
            X(m,randi(J)) = 1; 
        end

       
            

        % mu 2-D (M*J):
        for j = 1:J    % update according to column j
            % num: number of users asscociated with SBS j
            num = sum( X(:,j) );
            % 0. if there is no user accociated with SBS j 
            if (num == 0)
                % do nothing 

            % 1. if SBS j only associated with 1 user m, set mu_mj = 1
            elseif (num == 1)
                mu( X(:,j)==1, j ) = 1;

            % 2. SBS j associated with more than 1 user m, random mu_mj
            else
                % 2.1 random mu in former M-1 rows
                mu( X(1:end-1,j)==1, j ) = 1/num .* rand( sum( X(1:end-1,j)==1 ), 1 );
                % 2.2 update mu in last row: the sum of each column is guaranteed to be 1
                % 2.2.1 如果最后一个用户为1，正常更新
                if( X(end,j)==1 )
                    mu( end, j ) = 1 - sum( mu(1:end-1,j) );
                % 2.2.2 如果最后一个用户为0，需要将前面等于1的用户重新初始化
                else % x(end,j)==0
                    % index: a vector indicates the position of x==1(which user m associated with SBS j)
                    index = find( X(1:end,j)==1 );
                    % among users who associated with SBS j, random former end-1
                    mu( index(1:end-1), j ) = 1/num .* rand( num-1, 1 );
                    mu( index(end), j ) = 1 - sum( mu(index(1:end-1),j) );
                end    
            end        
        end
        
        % FL 2-D (M*K):
        FL(:,:) = (FLMax/K).*rand(M,K);

        % FE 3-D (M*K*J):
        FE(:,:,:) = (FEMax/(M*K)).*rand(M,K,J);

        % lambda 3-D (M*K*3):
        % Equality, just compute k1 and k2, then k3=1-k1-k2
        lambda(:,:,1:2) = 1/2 .* rand(M,K,2);  
        lambda(:,:,3) = 1 - sum(lambda,3);
    end
end