

% 把Position中第i个粒子的维度转换成决策变量矩阵，方便计算
% nargin  : number of function input arguments
% varargin: length input arguments list ( 1*n cell )

function [Pm,X,mu,FL,FE,lambda] = vector2matrix(varargin)

    if nargin ~= 3 && nargin ~= 4
        disp("!!!!Some mistake in vector2matrix.m, check the number of input variables.");        
        return
    
    elseif nargin == 3 % 如果输入变量个数为3，是FOM策略
        %% 1. 从varargin这个输入变量cell中取出变量
        M = varargin{1};
        J = varargin{2};
        Position = varargin{3};
        %% 2. 进行相关处理
        Pm = Position(1:M)'; % size(M*1)
        X  = reshape( Position(M+1:M+M*J),J,M )'; % 写两个reshape所有参数的函数吧
        mu = reshape( Position(M+M*J+1:M+M*J+M*J),J,M )'; % 在reshape成矩阵之后，把矩阵转置才能对应M*J的位置
        FL = 0;FE = 0;lambda = 0;
        
   
    else  % 如果输入变量个数为4，是 FOS || PO || 原始策略
        %% 1. 从varargin这个输入变量cell中取出变量
        M = varargin{1};
        K = varargin{2};
        J = varargin{3};
        Position = varargin{4};
        %% 2. 进行相关处理
        length = size(Position,2);
        
        % FOS策略
        if length == M*1 + M*J + M*J + M*K*J 
            Pm = Position(1:M)'; % size(M*1)
            X  = reshape( Position(M+1:M+M*J),J,M )'; % 写两个reshape所有参数的函数吧
            mu = reshape( Position(M+M*J+1:M+M*J+M*J),J,M )'; % 在reshape成矩阵之后，把矩阵转置才能对应M*J的位置
            FE = permute(... % permute():将多维矩阵的某些维度调换，这里需要将3维矩阵FE的2维1维调换，第3维保持不变
                 reshape( Position(M+M*J+M*J+1:M+M*J+M*J+M*K*J),K,M,J ),...
                 [2,1,3]); % FE size(M,K,J)
            FL = 0;lambda = 0;
            
        % PO策略 || 原始策略    
        elseif length == M*1 + M*J + M*J + M*K + M*K*J + M*K*3
            Pm = Position(1:M)'; % size(M*1)
            X  = reshape( Position(M+1:M+M*J),J,M )'; % 写两个reshape所有参数的函数吧
            mu = reshape( Position(M+M*J+1:M+M*J+M*J),J,M )'; % 在reshape成矩阵之后，把矩阵转置才能对应M*J的位置
            FL = reshape( Position(M+M*J+M*J+1:M+M*J+M*J+M*K),K,M )';
            FE = permute(... % permute():将多维矩阵的某些维度调换，这里需要将3维矩阵FE的2维1维调换，第3维保持不变
                 reshape( Position(M+M*J+M*J+M*K+1:M+M*J+M*J+M*K+M*K*J),K,M,J ),...
                 [2,1,3]); % FE size(M,K,J) 
            lambda = permute(... % permute():将多维矩阵的某些维度调换，这里需要将3维矩阵lambda的2维1维调换，第3维保持不变
                 reshape( Position(M+M*J+M*J+M*K+M*K*J+1:M+M*J+M*J+M*K+M*K*J+M*K*3),K,M,3 ),...
                 [2,1,3]); % lambda size(M,K,3)
        end
    end
end