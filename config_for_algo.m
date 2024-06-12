function cfg = config_for_algo(varargin)
    
    if nargin == 0   % 如果输入变量个数为0，说明是main_run_once.m调用了config()，只需要执行一遍程序
%        cfg.fileName = "GLPSO";
%        cfg.fileName = "GA"; 
%        cfg.fileName = "GWO";
         cfg.fileName = "LSAG";

        %% Basic params
        cfg.M = 10;      % user number M 5
        cfg.J = 3;      % SBS  number J 3
        cfg.K = 3;      % task number K 3 
    end
    
    cfg.execMode = "train";

    cfg.numOfParticles  = 100;   %粒子数
    cfg.totalIterations = 600; %一共的迭代次数
    cfg.amplify = 1;             %扩充
    
    if cfg.execMode == "train"
  
        
        
        cfg.user.I = reshape(...
            generateUniform(cfg.M*cfg.K, 500*1024*8, 1000*1024*8),cfg.M,cfg.K); %cfg.M*cfg.K=10*3 low = 500*1024*8 high=1000*1024*8
        % CPU cycles/bit for m,k                    
        cfg.user.alpha = reshape(...
            generateUniform(cfg.M*cfg.K, 20, 30),cfg.M,cfg.K);
        % memory/bit for m,k,                   
        cfg.user.g = reshape(...
            generateUniform(cfg.M*cfg.K, 50, 100),cfg.M,cfg.K);
        cfg.user.sigma = (generateUniform(cfg.M, 0.8*10^(-25), 1.3*10^(-25)))';
        cfg.user.d = reshape(...
            generateUniform(cfg.M*cfg.J, 20, 100),cfg.M,cfg.J);      
        cfg.user.Pm0 = (generateUniform(cfg.M, 0.2, 0.8))'; 
        cfg.user.Pmr = (generateUniform(cfg.M, 0.2, 0.8))'; 
        cfg.user.rho = (generateUniform(cfg.M, 16, 20))'; 
        
        
        
     
        cfg.SBS.Pf = generateUniform(cfg.J, 0.1, 0.4); 
        cfg.SBS.BWup = generateUniform(cfg.J, 1*10^(7), 2*10^(7));        
        cfg.SBS.BWdown = generateUniform(cfg.J, 1*10^(7), 2*10^(7));          
        cfg.SBS.sigma = generateUniform(cfg.J, 0.5*10^(-26), 1*10^(-26)); 
        cfg.SBS.Rup = zeros(cfg.M,cfg.J); 
        cfg.SBS.Rdown = zeros(cfg.M,cfg.J);
            
        
        
        
         cfg.MBS.f0 = 5*10^(9);
         cfg.MBS.r0 = 1024*1024*1024*8;           
         cfg.MBS.e0 = 1*10^(-9);
         

         
         cfg.cons.v = 4;            
         cfg.cons.h1 = 0.73;        
         cfg.cons.h2 = 0.73;        
         cfg.cons.N0 = 10^(-8)/625; 
         cfg.cons.Psigma1 = 0.001;  
         cfg.cons.Psigma2 = 0.001;  
         cfg.cons.beta1 = 1;        
         cfg.cons.beta2 = 0.3;            
         cfg.cons.beta3 = 1;       
         cfg.cons.beta4 = 0.3;      

         cfg.upbound.PmMax = 0.1;       
         cfg.upbound.FLMax = 10^(9);    
         cfg.upbound.FEMax = 5*10^(9);  
         cfg.upbound.LMax  = (generateUniform(cfg.M, 3, 7))';    
         cfg.upbound.EmMax = (generateUniform(cfg.M, 3, 6))';    
         cfg.upbound.EjMax = 13;                     
         cfg.upbound.E0Max = 2*cfg.upbound.EjMax;    
         cfg.upbound.ASBSMax = 2*10^(10);               
         cfg.upbound.GSBSMax = 2*1024*1024*1024*8;     
         cfg.upbound.AMBSMax = 3*cfg.upbound.ASBSMax;   
         cfg.upbound.GMBSMax = 3*cfg.upbound.GSBSMax;   
          
    end

