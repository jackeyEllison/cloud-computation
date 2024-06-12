
% 把决策变量矩阵转换成Position中第i个粒子的维度
% return Position is a row vector size of (1*numOfDimention)
% nargin  : number of function input arguments
% varargin: length input arguments list ( 1*n cell )

function Position = matrix2vextor(varargin)
 
    if nargin ~= 5 && nargin ~= 9 && nargin ~= 7
        disp("!!!!Some mistake in matrix2vector.m, check the number of input variables.");        
        return
    elseif nargin == 5  % 如果输入变量个数为5，是FOM策略
        %% 1. 从varargin这个输入变量cell中取出变量
        M = varargin{1};
        J = varargin{2};
        Pm = varargin{3};
        X = varargin{4};
        mu = varargin{5};
        %% 2. 进行相关处理
        Position = zeros(1,M*1+M*J+M*J);
        Position(1:M) = Pm';
        Position(M+1:M+M*J) = reshape(X',1,M*J); % 需要先将决策变量矩阵转置后才能对应上M*J位置
        Position(M+M*J+1:M+M*J+M*J) = reshape(mu',1,M*J);
    
        
    elseif nargin == 7  % 如果输入变量个数为7，是FOS策略
        %% 1. 从varargin这个输入变量cell中取出变量
        M = varargin{1};
        K = varargin{2};
        J = varargin{3};
        Pm = varargin{4};
        X = varargin{5};
        mu = varargin{6}; 
        FE = varargin{7};
        %% 2. 进行相关处理
        Position = zeros(1,M*1+M*J+M*J+M*K*J);
        Position(1:M) = Pm';
        Position(M+1:M+M*J) = reshape(X',1,M*J); % 需要先将决策变量矩阵转置后才能对应上M*J位置
        Position(M+M*J+1:M+M*J+M*J) = reshape(mu',1,M*J);
        Position(M+M*J+M*J+1:M+M*J+M*J+M*K*J) = ...
            reshape(permute(FE,[2,1,3]),1,M*K*J); % permute(FE,[2,1,3])将FE矩阵的前两维交换

        
    elseif nargin == 9  % 如果输入变量个数为9，是 PO策略 || 原始策略
        %% 1. 从varargin这个输入变量cell中取出变量
        M = varargin{1};
        K = varargin{2};
        J = varargin{3};
        Pm = varargin{4};
        X = varargin{5};
        mu = varargin{6}; 
        FL = varargin{7};
        FE = varargin{8};
        lambda = varargin{9};
        %% 2. 进行相关处理
        Position = zeros(1,M*1+M*J+M*J+M*K+M*K*J+M*K*3);
        Position(1:M) = Pm';
        Position(M+1:M+M*J) = reshape(X',1,M*J); % 需要先将决策变量矩阵转置后才能对应上M*J位置
        Position(M+M*J+1:M+M*J+M*J) = reshape(mu',1,M*J);
        Position(M+M*J+M*J+1:M+M*J+M*J+M*K) = reshape(FL',1,M*K);
        Position(M+M*J+M*J+M*K+1:M+M*J+M*J+M*K+M*K*J) = ...
            reshape(permute(FE,[2,1,3]),1,M*K*J); % permute(FE,[2,1,3])将FE矩阵的前两维交换
        Position(M+M*J+M*J+M*K+M*K*J+1:end) = ...
            reshape(permute(lambda,[2,1,3]),1,M*K*3); % permute(lambda,[2,1,3])将lambda矩阵的前两维交换 

    end
end