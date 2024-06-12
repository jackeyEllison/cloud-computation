%% OBJFUNCTION  1. Establish the system(model) of 
%               Multi user Multi base_station(MM) 
%               2. Compute augment/sum_energy/penalty 
% 
%  Input:       
%               cfg: 
%                   configuration of system
%               optimized results:
%                   (1-D)Pm,...
%                   (2-D)X,mu,FL,...
%                   (3-D)FE,lambda 
%  Output:      augmentedObj = sum_energy + penalty -> resultObj
%               sum_energy                          -> returnProfit
%               penalty                             -> returnCurrentPenalty

% 注意：只有在主函数中调用的时候，才有可能用到dataSMD及后面的返回值，否则返回值只有前三个
function [augmentedObj,sum_energy,penalty,...
    dataSMD,dataSBS,dataMBS,energySMD,energySBS,E0] = ...
    objFunction(varargin)
    if nargin ~= 4 && nargin ~= 7 && nargin ~= 5
        disp("!!!!Some mistake in objFunction.m, check the number of input variables.");        
        return
    elseif nargin == 4  % 如果输入变量个数为4，是FOM策略 
        cfg = varargin{1};
        Pm = varargin{2};
        X = varargin{3};
        mu = varargin{4};
        energySMD = 0;energySBS = 0;dataSMD = 0;dataSBS = 0;dataMBS = 0;
        %% ASSIGN
        % lambda for local/SBS/MBS is k1/k2/k3 
        k2 = 0;     % size(M*K) all tasks are offloaded to MBS 
        k3 = 1;     % size(M*K) 

        % Basic params
        J = cfg.J;

        % user  2-D M*K
        I       = cfg.user.I;       % size of task k of user m 
        alpha   = cfg.user.alpha;   % CPU cycles/bit
        g       = cfg.user.g;       % memory/bit

        % user  1-D M*1
        d       = cfg.user.d;       % distance between user m and SBS it connected
        Pm0     = cfg.user.Pm0;     % idle power of user m  
        Pmr     = cfg.user.Pmr;     % download data power of user m 
        rho     = cfg.user.rho;     % amplifier coefficient of user m

        % SBS   1-D 1*J
        Pf      = cfg.SBS.Pf;       % Pf: transmitted power of SBS
        BWup    = cfg.SBS.BWup;     % BWup: uplink channel bandwidth of SBS
        BWdown  = cfg.SBS.BWdown;   % BWdown: downlink channel bandwidth of SBS

        % SBS   2-D M*J
        Rup     = cfg.SBS.Rup;      % data_rate of uplink/downlink for m,j
        Rdown   = cfg.SBS.Rdown;    % data_rate of uplink/downlink for m,j

        % MBS   scalar
        f0       = cfg.MBS.f0;       % f0: MBS computing speed
        r0       = cfg.MBS.r0;       % r0: data_rate between MBS and SBS
        e0       = cfg.MBS.e0;       % e0: energy consumption per CPU cycle 

        % constant
        v        = cfg.cons.v;        % exponent of path loss 
        h1       = cfg.cons.h1;       % fading coefficient of uplink 
        h2       = cfg.cons.h2;       % fading coefficient of downlink 
        N0       = cfg.cons.N0;       % White noise 
        Psigma1  = cfg.cons.Psigma1;  % from SBS to MBS transmission power 
        Psigma2  = cfg.cons.Psigma2;  % from MBS to SBS transmission power 
        beta1    = cfg.cons.beta1;    % uplink overhead coefficient between user and SBS 
        beta2    = cfg.cons.beta2;    % downlink overhead coefficient between user and SBS 
        beta3    = cfg.cons.beta3;    % uplink overhead coefficient between SBS and MBS 
        beta4    = cfg.cons.beta4;    % downlink overhead coefficient between SBS and MBS 

        %% Upper bound
        PmMax   = cfg.upbound.PmMax;   % PmMax: max transmission power of user m                   0.1W
        LMax    = cfg.upbound.LMax;    % LMax: max latency of user m      ~U[5,10]s
        E0Max   = cfg.upbound.E0Max;   % E0Max: max energy of MBS                      10J
        AMBSMax = cfg.upbound.AMBSMax; % AMBSMax: max CPU cycles of MBS    1*10^(10)
        GMBSMax = cfg.upbound.GMBSMax; % GMBSMax: max memory of MBS        4GB

        inequality = 0;
%         equality=0;
        sum_energy = 0;
        augmentedObj = 0;
        penalty = 0;
        amplify = cfg.amplify;

        %% COMPUTE

        %% local latency
        each_user_local_time = 0;            % size(M*1)

        %% SBS latency
        % data rate of uplink/downlink from user to SBS
        Rup   = Rup + mu .* BWup .*   log2( 1 + (Pm.*d.^(-v).*abs(h1)^2)./N0 );   % size(M*J) 
        Rdown = Rdown+mu .* BWdown .* log2( 1 + (Pf.*d.^(-v).*abs(h2)^2)./N0 );   % size(M*J) 

        % 将稀疏矩阵Rup（每一行只有一个元素不为0）
        % 变成稠密矩阵Rup_dense（把每行中不为0的元素提取出来成为一个M*1的向量）
        Rup_dense   = sparse2dense(Rup);        % size(M*1)
        Rdown_dense = sparse2dense(Rdown);      % size(M*1)
        mu_dense    = sparse2dense(mu);         % size(M*1)

        % Tup/Tdown: 每一个用户m上行/下行链路花费的时间 这里k2=0
        Tup   = sum( beta1.*(k2+k3).*I, 2 ) ./ Rup_dense;            % size(M*1)
        Tdown = sum( beta2.*(k2+k3).*I, 2 ) ./ Rdown_dense;          % size(M*1)
%         Tup
%         Tdown
        % SBS的时间分为三项：处理时间、上行链路时间、下行链路时间
        % 这里要注意！SBS_time是M*1维的向量不是J*1，只是因为时间约束关于用户m的时间约束
        % 这里time_need_in_SBS = 0
        SBS_time = Tup + Tdown;                   % size(M*1)


        %% MBS latency
        % each_user_need_MBS_time : 每一个用户分配到MBS的任务需要处理的时间
        % each_MBS_uplink_time    : 每一个MBS上行链路的时间
        % each_MBS_downlink_time  : 每一个MBS下行链路的时间
        cpu_need_in_MBS = I .* k3 .* alpha;                         % size(M*K) 
        AMBS = sum(sum( cpu_need_in_MBS ));                         % size(1*1) 
        GMBS = sum(sum(   I .* k3 .* g  ));                         % size(1*1)    

        each_user_need_MBS_time = sum( cpu_need_in_MBS./f0 , 2 );   % size(M*1)
        MBS_uplink_time         = sum( beta3.*k3.*I, 2 ) ./ r0 ;    % size(M*1)
        MBS_downlink_time       = sum( beta4.*k3.*I, 2 ) ./ r0 ;    % size(M*1)
%         each_user_need_MBS_time
%         MBS_uplink_time
%         MBS_downlink_time
        % MBS的时间分为三项：处理时间、上行链路时间、下行链路时间
        MBS_time = each_user_need_MBS_time + ...                    % size(M*1)
            MBS_uplink_time + MBS_downlink_time;


        %% local energy
        % Em: each user local energy
        Em = 0;                 % size(M*1)

        %% SBS energy
        each_user_uplink_energy   = Tup .* ( Pm0+rho.*mu_dense.*Pm );    % size(M*1)
        each_user_downlink_energy = Tdown .* Pmr                    ;    % size(M*1)

        each_SBS_uplink_energy    = zeros(J,1);                          % size(J*1)
        each_SBS_downlink_energy  = zeros(J,1);                          % size(J*1)
        
        % 需要将用户(M*1)表示的能量转换为基站(J*1)表示的能量
        for j=1:J
            % m是基站j所连接的全部用户(scalar or vector) num是基站j所连接的全部用户的个数
            num = sum( X(:,j) );
            m = find( X(:,j)==1 );
            if num == 0 % 没有用户和基站j连接
                each_SBS_uplink_energy(j)  = 0;
                each_SBS_downlink_energy(j)= 0;
            else
                each_SBS_uplink_energy(j) = sum( each_user_uplink_energy(m) );
                each_SBS_downlink_energy(j) = sum( each_user_downlink_energy(m) );            
            end
        end
%         each_SBS_uplink_energy
%         each_SBS_downlink_energy
        % Ej: each SBS j total energy: 3 parts : uplink/downlink/compute
        Ej = each_SBS_uplink_energy + each_SBS_downlink_energy;               % size(J*1)
%         Ej

        %% MBS energy
        MBS_uplink_energy   = sum(sum( MBS_uplink_time   .* Psigma1 ));     % size(1*1)
        MBS_downlink_energy = sum(sum( MBS_downlink_time .* Psigma2 ));     % size(1*1)
        MBS_compute_energy  = sum(sum( cpu_need_in_MBS   .* e0      ));     % size(1*1)
%         MBS_uplink_energy
%         MBS_downlink_energy
%         MBS_compute_energy
        % 注意E0是一个scalar！
        E0 = MBS_compute_energy + ...
            MBS_uplink_energy + MBS_downlink_energy;                        % size(1*1)

        %%
        % Tm: user m total time 
        Tm = max( each_user_local_time, (SBS_time+MBS_time) );       % size(M*1)
        sum_energy = sum_energy + sum(Em(:)) + sum(Ej(:)) + E0;                  % total energy size(1*1)

        %% penalty 
        inequality = inequality + sum( max(Tm - LMax, 0).^2 ) + ... % latency penalty
                 max(AMBS-AMBSMax, 0).^2   + ... % latency penalty
                 max(GMBS-GMBSMax, 0).^2   + ... % latency penalty
                 max(E0-E0Max, 0).^2   + ... % latency penalty
            sum( max(Pm-PmMax, 0).^2 ) + ...            % transmission power penalty
            sum(sum(     max(mu    -1, 0).^2 ));

        % 注意这里不需要计算等式的约束(equality=0)，因为在等式变量lambda和mu更新的过程中
        % 已经满足了相应的等式条件
        penalty = penalty + inequality; 
        augmentedObj = augmentedObj + sum_energy + amplify*penalty;
     
        
     elseif nargin == 5  % 如果输入变量个数为5，是FOS策略
        cfg = varargin{1};
        Pm = varargin{2};
        X  = varargin{3};
        mu = varargin{4};
        FE = varargin{5};
        energySMD = 0;energySBS = 0;dataSMD = 0;dataSBS = 0;dataMBS = 0;
        %% ASSIGN
        % lambda for local/SBS/MBS is k1/k2/k3 
        k2 = 1;     % size(M*K) all tasks are offloaded to SBS : k2 = 1
        k3 = 0;     % size(M*K) 

        % Basic params
        M = cfg.M;
        J = cfg.J;
        K = cfg.K;

        % user  2-D M*K
        I       = cfg.user.I;       % size of task k of user m 
        alpha   = cfg.user.alpha;   % CPU cycles/bit
        g       = cfg.user.g;       % memory/bit

        % user  1-D M*1
        d       = cfg.user.d;       % distance between user m and SBS it connected
        Pm0     = cfg.user.Pm0;     % idle power of user m  
        Pmr     = cfg.user.Pmr;     % download data power of user m 
        rho     = cfg.user.rho;     % amplifier coefficient of user m

        % SBS   1-D 1*J
        Pf      = cfg.SBS.Pf;       % Pf: transmitted power of SBS
        BWup    = cfg.SBS.BWup;     % BWup: uplink channel bandwidth of SBS
        BWdown  = cfg.SBS.BWdown;   % BWdown: downlink channel bandwidth of SBS
        sgmEg   = cfg.SBS.sigma;    % sgmEg: depend on chip architecture of j

        % SBS   2-D M*J
        Rup     = cfg.SBS.Rup;      % data_rate of uplink/downlink for m,j
        Rdown   = cfg.SBS.Rdown;    % data_rate of uplink/downlink for m,j

        % constant
        v        = cfg.cons.v;        % exponent of path loss 
        h1       = cfg.cons.h1;       % fading coefficient of uplink 
        h2       = cfg.cons.h2;       % fading coefficient of downlink 
        N0       = cfg.cons.N0;       % White noise 
        beta1    = cfg.cons.beta1;    % uplink overhead coefficient between user and SBS 
        beta2    = cfg.cons.beta2;    % downlink overhead coefficient between user and SBS 

        %% Upper bound
        PmMax   = cfg.upbound.PmMax;   % PmMax: max transmission power of user m                   0.1W
        FEMax   = cfg.upbound.FEMax;   % FEMax: max edge cpu processing speed of SBS j             5GHz(cycels/sec)
        LMax    = cfg.upbound.LMax;    % LMax: max latency of user m      ~U[5,10]s
        EjMax   = cfg.upbound.EjMax;   % EjMax: max energy of SBS j                    10J
        ASBSMax = cfg.upbound.ASBSMax; % ASBSMax: max CPU cycles of SBS j         5*10^(9)
        GSBSMax = cfg.upbound.GSBSMax; % GSBSMax: max memory of SBS j      2GB

        inequality = 0;
%         equality=0;
        sum_energy = 0;
        augmentedObj = 0;
        penalty = 0;
        amplify = cfg.amplify;

        %% COMPUTE
        %% local latency
        each_user_local_time = 0;            % size(M*1)
        
        %% SBS latency
        % time_need_in_SBS : 每一个用户分配到SBS的任务需要处理的时间 size(M*1)
        FE_new = zeros(M,K);
        for m = 1:M
            % X(m,:)==1 是第m个用户选择的第几个SBSj在矩阵中的位置
            % FE_new是FE压缩之后的矩阵，因为在运算的时候FE中的大部分维数是没有用到的
            FE_new(m,:) = FE(m,:, X(m,:)==1 );                              % size(M*K)
        end
        % cpu_need_in_SBS : the CPU cycles needed in SBS for each user m
        cpu_need_in_SBS    = I .* k2 .* alpha;                              % size(M*K)    
        memory_need_in_SBS = I .* k2 .*  g;                                 % size(M*K)    
        time_need_in_SBS = sum( cpu_need_in_SBS./FE_new , 2 );              % size(M*1)

        % data rate of uplink/downlink from user to SBS
        Rup   = Rup + mu .* BWup .*   log2( 1 + (Pm.*d.^(-v).*abs(h1)^2)./N0 );   % size(M*J) 
        Rdown = Rdown+mu .* BWdown .* log2( 1 + (Pf.*d.^(-v).*abs(h2)^2)./N0 );   % size(M*J) 

        % 将稀疏矩阵Rup（每一行只有一个元素不为0）
        % 变成稠密矩阵Rup_dense（把每行中不为0的元素提取出来成为一个M*1的向量）
        Rup_dense   = sparse2dense(Rup);        % size(M*1)
        Rdown_dense = sparse2dense(Rdown);      % size(M*1)
        mu_dense    = sparse2dense(mu);         % size(M*1)

        % Tup/Tdown: 每一个用户m上行/下行链路花费的时间
        Tup   = sum( beta1.*(k2+k3).*I, 2 ) ./ Rup_dense;            % size(M*1)
        Tdown = sum( beta2.*(k2+k3).*I, 2 ) ./ Rdown_dense;          % size(M*1)

        % SBS的时间分为三项：处理时间、上行链路时间、下行链路时间
        % 这里要注意！SBS_time是M*1维的向量不是J*1，只是因为时间约束关于用户m的时间约束
        SBS_time = time_need_in_SBS + Tup + Tdown;                   % size(M*1)


        %% MBS latency
        MBS_time = 0;


        %% local energy
        % Em: each user local energy
        Em = 0;                 % size(M*1)

        %% SBS energy
        each_user_uplink_energy   = Tup .* ( Pm0+rho.*mu_dense.*Pm );    % size(M*1)
        each_user_downlink_energy = Tdown .* Pmr                    ;    % size(M*1)

        each_SBS_compute_energy   = zeros(J,1);                          % size(J*1)
        each_SBS_uplink_energy    = zeros(J,1);                          % size(J*1)
        each_SBS_downlink_energy  = zeros(J,1);                          % size(J*1)
        ASBS = zeros(J,1);                                               % size(J*1)
        GSBS = zeros(J,1);                                               % size(J*1)
%         X
%         cpu_need_in_SBS
        
        % 需要将用户(M*1)表示的能量转换为基站(J*1)表示的能量
        for j=1:J
            % m是基站j所连接的全部用户(scalar or vector) num是基站j所连接的全部用户的个数
            num = sum( X(:,j) );
            m = find( X(:,j)==1 );
            if num == 0 % 没有用户和基站j连接
                each_SBS_compute_energy(j) = 0;
                each_SBS_uplink_energy(j)  = 0;
                each_SBS_downlink_energy(j)= 0;
                ASBS(j) = 0;
                GSBS(j) = 0;
            else
                each_SBS_compute_energy(j) = ... 
                    sum(sum( cpu_need_in_SBS(m,:) .* FE(m,:,j) .* FE(m,:,j) .* sgmEg(j) ));
                each_SBS_uplink_energy(j) = sum( each_user_uplink_energy(m) );
                each_SBS_downlink_energy(j) = sum( each_user_downlink_energy(m) );            
                ASBS(j) = sum(sum( cpu_need_in_SBS(m,:)    ));
                GSBS(j) = sum(sum( memory_need_in_SBS(m,:) ));
            end
        end

        % Ej: each SBS j total energy: 3 parts : uplink/downlink/compute
        Ej = each_SBS_compute_energy + ...
            each_SBS_uplink_energy + each_SBS_downlink_energy;               % size(J*1)


        %% MBS energy
        E0 = 0;                        % size(1*1)


        %%
        % Tm: user m total time 
        Tm = max( each_user_local_time, (SBS_time+MBS_time) );       % size(M*1)
        sum_energy = sum_energy + sum(Em(:)) + sum(Ej(:)) + E0;                  % total energy size(1*1)


        %% penalty 
        inequality = inequality + sum( max(Tm - LMax, 0).^2 ) + ... % latency penalty
            sum( max(ASBS-ASBSMax, 0).^2 ) + ... % latency penalty
            sum( max(GSBS-GSBSMax, 0).^2 ) + ... % latency penalty
            sum( max(Ej-EjMax, 0).^2 ) + ... % latency penalty
            sum( max(Pm-PmMax, 0).^2 ) + ...            % transmission power penalty
            sum(sum(     max(mu    -1, 0).^2 ));

        if K ~=1 % 只有一个K的时候求和的时候有bug，用条件语句解决
            inequality = inequality + sum( max(sum(sum(FE))-FEMax, 0).^2 );      % transmission power penalty
        else
            inequality = inequality + sum( max(sum(FE)-FEMax, 0).^2 );      % transmission power penalty        
        end

        % 注意这里不需要计算等式的约束(equality=0)，因为在等式变量lambda和mu更新的过程中
        % 已经满足了相应的等式条件
        penalty = penalty+inequality; 
        augmentedObj = augmentedObj+sum_energy + amplify*penalty;
        
        
        
    
    elseif nargin == 7  % 如果输入变量个数为7，是 PO策略 || 原始策略
        cfg = varargin{1};
        Pm = varargin{2};
        X  = varargin{3};
        mu = varargin{4};
        FL = varargin{5};
        FE = varargin{6};
        lambda = varargin{7};
        
        %% ASSIGN
        % lambda for local/SBS/MBS is k1/k2/k3 
        k1 = lambda(:,:,1);     % size(M*K) 
        k2 = lambda(:,:,2);     % size(M*K) 
        k3 = lambda(:,:,3);     % size(M*K) 

        % Basic params
        M = cfg.M;
        J = cfg.J;
        K = cfg.K;

        % user  2-D M*K
        I       = cfg.user.I;       % size of task k of user m 
        alpha   = cfg.user.alpha;   % CPU cycles/bit
        g       = cfg.user.g;       % memory/bit

        % user  1-D M*1
        sgmLc   = cfg.user.sigma;   % chip architecture coefficient
        d       = cfg.user.d;       % distance between user m and SBS it connected
        Pm0     = cfg.user.Pm0;     % idle power of user m  
        Pmr     = cfg.user.Pmr;     % download data power of user m 
        rho     = cfg.user.rho;     % amplifier coefficient of user m

        % SBS   1-D 1*J
        Pf      = cfg.SBS.Pf;       % Pf: transmitted power of SBS
        BWup    = cfg.SBS.BWup;     % BWup: uplink channel bandwidth of SBS
        BWdown  = cfg.SBS.BWdown;   % BWdown: downlink channel bandwidth of SBS
        sgmEg   = cfg.SBS.sigma;    % sgmEg: depend on chip architecture of j

        % SBS   2-D M*J
        Rup     = cfg.SBS.Rup;      % data_rate of uplink/downlink for m,j
        Rdown   = cfg.SBS.Rdown;    % data_rate of uplink/downlink for m,j

        % MBS   scalar
        f0       = cfg.MBS.f0;       % f0: MBS computing speed
        r0       = cfg.MBS.r0;       % r0: data_rate between MBS and SBS
        e0       = cfg.MBS.e0;       % e0: energy consumption per CPU cycle 

        % constant
        v        = cfg.cons.v;        % exponent of path loss 
        h1       = cfg.cons.h1;       % fading coefficient of uplink 
        h2       = cfg.cons.h2;       % fading coefficient of downlink 
        N0       = cfg.cons.N0;       % White noise 
        Psigma1  = cfg.cons.Psigma1;  % from SBS to MBS transmission power 
        Psigma2  = cfg.cons.Psigma2;  % from MBS to SBS transmission power 
        beta1    = cfg.cons.beta1;    % uplink overhead coefficient between user and SBS 
        beta2    = cfg.cons.beta2;    % downlink overhead coefficient between user and SBS 
        beta3    = cfg.cons.beta3;    % uplink overhead coefficient between SBS and MBS 
        beta4    = cfg.cons.beta4;    % downlink overhead coefficient between SBS and MBS 

        %% Upper bound
        PmMax   = cfg.upbound.PmMax;   % PmMax: max transmission power of user m                   0.1W
        FLMax   = cfg.upbound.FLMax;   % FLMax: max local cpu processing speed user m              1GHz(cycels/sec)
        FEMax   = cfg.upbound.FEMax;   % FEMax: max edge cpu processing speed of SBS j             5GHz(cycels/sec)
        LMax    = cfg.upbound.LMax;    % LMax: max latency of user m      ~U[5,10]s
        EmMax   = cfg.upbound.EmMax;   % EmMax: max energy of user m      ~U[1,4]J
        EjMax   = cfg.upbound.EjMax;   % EjMax: max energy of SBS j                    10J
        E0Max   = cfg.upbound.E0Max;   % E0Max: max energy of MBS                      10J
        ASBSMax = cfg.upbound.ASBSMax; % ASBSMax: max CPU cycles of SBS j         5*10^(9)
        GSBSMax = cfg.upbound.GSBSMax; % GSBSMax: max memory of SBS j      2GB
        AMBSMax = cfg.upbound.AMBSMax; % AMBSMax: max CPU cycles of MBS    1*10^(10)
        GMBSMax = cfg.upbound.GMBSMax; % GMBSMax: max memory of MBS        4GB

        inequality = 0;
        equality=0;
        sum_energy = 0;
        augmentedObj = 0;
        penalty = 0;
        amplify = cfg.amplify;

        %% COMPUTE

        %% local latency
        cpu_need_in_local = I .* k1 .* alpha;                               % size(M*K) 
        each_user_local_time = sum( cpu_need_in_local./FL , 2 );            % size(M*1)

        %% SBS latency
        % time_need_in_SBS : 每一个用户分配到SBS的任务需要处理的时间 size(M*1)
        FE_new = zeros(M,K);
        for m = 1:M
            % X(m,:)==1 是第m个用户选择的第几个SBSj在矩阵中的位置
            % FE_new是FE压缩之后的矩阵，因为在运算的时候FE中的大部分维数是没有用到的
            FE_new(m,:) = FE(m,:, X(m,:)==1 );                              % size(M*K)
        end
        % cpu_need_in_SBS : the CPU cycles needed in SBS for each user m
        cpu_need_in_SBS    = I .* k2 .* alpha;                              % size(M*K)    
        memory_need_in_SBS = I .* k2 .*  g;                                 % size(M*K)    
        time_need_in_SBS = sum( cpu_need_in_SBS./FE_new , 2 );              % size(M*1)

        % data rate of uplink/downlink from user to SBS
        Rup   = Rup + mu .* BWup .*   log2( 1 + (Pm.*d.^(-v).*abs(h1)^2)./N0 );   % size(M*J) 
        Rdown = Rdown+mu .* BWdown .* log2( 1 + (Pf.*d.^(-v).*abs(h2)^2)./N0 );   % size(M*J) 

        % 将稀疏矩阵Rup（每一行只有一个元素不为0）
        % 变成稠密矩阵Rup_dense（把每行中不为0的元素提取出来成为一个M*1的向量）
        Rup_dense   = sparse2dense(Rup);        % size(M*1)
        Rdown_dense = sparse2dense(Rdown);      % size(M*1)
        mu_dense    = sparse2dense(mu);         % size(M*1)

        % Tup/Tdown: 每一个用户m上行/下行链路花费的时间
        Tup   = sum( beta1.*(k2+k3).*I, 2 ) ./ Rup_dense;            % size(M*1)
        Tdown = sum( beta2.*(k2+k3).*I, 2 ) ./ Rdown_dense;          % size(M*1)

        % SBS的时间分为三项：处理时间、上行链路时间、下行链路时间
        % 这里要注意！SBS_time是M*1维的向量不是J*1，只是因为时间约束关于用户m的时间约束
        SBS_time = time_need_in_SBS + Tup + Tdown;                   % size(M*1)


        %% MBS latency
        % each_user_need_MBS_time : 每一个用户分配到MBS的任务需要处理的时间
        % each_MBS_uplink_time    : 每一个MBS上行链路的时间
        % each_MBS_downlink_time  : 每一个MBS下行链路的时间
        cpu_need_in_MBS = I .* k3 .* alpha;                         % size(M*K) 
        AMBS = sum(sum( cpu_need_in_MBS ));                         % size(1*1) 
        GMBS = sum(sum(   I .* k3 .* g  ));                         % size(1*1)    

        each_user_need_MBS_time = sum( cpu_need_in_MBS./f0 , 2 );   % size(M*1)
        MBS_uplink_time         = sum( beta3.*k3.*I, 2 ) ./ r0 ;    % size(M*1)
        MBS_downlink_time       = sum( beta4.*k3.*I, 2 ) ./ r0 ;    % size(M*1)

        % MBS的时间分为三项：处理时间、上行链路时间、下行链路时间
        MBS_time = each_user_need_MBS_time + ...                    % size(M*1)
            MBS_uplink_time + MBS_downlink_time;


        %% local energy
        % Em: each user local energy
        Em = sum( sgmLc.*cpu_need_in_local.*FL.*FL , 2 );                 % size(M*1)


        %% SBS energy
        each_user_uplink_energy   = Tup .* ( Pm0+rho.*mu_dense.*Pm );    % size(M*1)
        each_user_downlink_energy = Tdown .* Pmr                    ;    % size(M*1)

        each_SBS_compute_energy   = zeros(J,1);                          % size(J*1)
        each_SBS_uplink_energy    = zeros(J,1);                          % size(J*1)
        each_SBS_downlink_energy  = zeros(J,1);                          % size(J*1)
        ASBS = zeros(J,1);                                               % size(J*1)
        GSBS = zeros(J,1);                                               % size(J*1)

        % 需要将用户(M*1)表示的能量转换为基站(J*1)表示的能量
        for j=1:J
            % m是基站j所连接的全部用户(scalar or vector) num是基站j所连接的全部用户的个数
            num = sum( X(:,j) );
            m = find( X(:,j)==1 );
            if num == 0 % 没有用户和基站j连接
                each_SBS_compute_energy(j) = 0;
                each_SBS_uplink_energy(j)  = 0;
                each_SBS_downlink_energy(j)= 0;
                ASBS(j) = 0;
                GSBS(j) = 0;
            else
                each_SBS_compute_energy(j) = ... 
                    sum(sum( cpu_need_in_SBS(m,:) .* FE(m,:,j) .* FE(m,:,j) .* sgmEg(j) ));
                each_SBS_uplink_energy(j) = sum( each_user_uplink_energy(m) );
                each_SBS_downlink_energy(j) = sum( each_user_downlink_energy(m) );            
                ASBS(j) = sum(sum( cpu_need_in_SBS(m,:)    ));
                GSBS(j) = sum(sum( memory_need_in_SBS(m,:) ));
            end
        end

        % Ej: each SBS j total energy: 3 parts : uplink/downlink/compute
        Ej = each_SBS_compute_energy + ...
            each_SBS_uplink_energy + each_SBS_downlink_energy;               % size(J*1)


        %% MBS energy
        MBS_uplink_energy   = sum(sum( MBS_uplink_time   .* Psigma1 ));     % size(1*1)
        MBS_downlink_energy = sum(sum( MBS_downlink_time .* Psigma2 ));     % size(1*1)
        MBS_compute_energy  = sum(sum( cpu_need_in_MBS   .* e0      ));     % size(1*1)

        % 注意E0是一个scalar！
        E0 = MBS_compute_energy + ...
            MBS_uplink_energy + MBS_downlink_energy;                        % size(1*1)


        %%
        % Tm: user m total time 
        Tm = max( each_user_local_time, (SBS_time+MBS_time) );              % size(M*1)
        dataSMD   = sum(sum(k1.*I));
        dataSBS   = sum(sum(k2.*I));
        dataMBS   = sum(sum(k3.*I));
        energySMD = sum(Em(:));
        energySBS = sum(Ej(:));
        sum_energy = sum_energy + energySMD + energySBS + E0;             % total energy size(1*1)


        %% penalty 
        inequality = inequality + sum( max(Tm - LMax, 0).^2 ) + ... % latency penalty
            sum( max(ASBS-ASBSMax, 0).^2 ) + ... % latency penalty
            sum( max(GSBS-GSBSMax, 0).^2 ) + ... % latency penalty
                 max(AMBS-AMBSMax, 0).^2   + ... % latency penalty
                 max(GMBS-GMBSMax, 0).^2   + ... % latency penalty
            sum( max(Em-EmMax, 0).^2 ) + ... % latency penalty
            sum( max(Ej-EjMax, 0).^2 ) + ... % latency penalty
                 max(E0-E0Max, 0).^2   + ... % latency penalty
            sum( max(Pm-PmMax, 0).^2 ) + ...            % transmission power penalty
            sum( max(sum(FL,2)-FLMax, 0).^2 ) + ...     % transmission power penalty
            sum(sum(     max(mu    -1, 0).^2 )) + ...
            sum(sum(sum( max(lambda-1, 0).^2 )));

        if K~=1 && M~=1 % M==1 || K==1 的时候求和的时候有bug，用条件语句解决(只需要求一次和或者不求和)
            inequality = inequality + sum( max(sum(sum(FE))-FEMax, 0).^2 );
        else
            if K==1 && M==1
                inequality = inequality + sum( max(FE-FEMax, 0).^2 );
            else
                inequality = inequality + sum( max(sum(FE)-FEMax, 0).^2 );  
            end
        end

        % 注意这里不需要计算等式的约束(equality=0)，因为在等式变量lambda和mu更新的过程中
        % 已经满足了相应的等式条件
        penalty = penalty + inequality; 
        augmentedObj = augmentedObj + sum_energy + amplify*penalty;

    end
end