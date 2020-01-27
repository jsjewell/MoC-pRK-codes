% This FUNCTION code outputs the MoC-RK3 solution
% (see Sec. 3.1 in Lakoba-Jewell paper).
% Inputs for this function are:
%  - the solution at the current time level, and
%  - the boundary conditions at the new level. 
	
function  [ypp_n, ymm_n] = LakobaJewell_MoC_pRK4_FUNC_MoCRK3( ...
                           M, dt, NumComp, fpp, fmm, ypp, ymm, BCleft, BCright );
					
	% Compute the stage derivative k1:
    for   jNC = 1 : NumComp
        k1_pp(jNC, 1:M-1) = fpp{jNC}( ypp(:, 1:M-1), ymm(:, 1:M-1) );
        k1_mm(jNC, 2 : M) = fmm{jNC}( ypp(:, 2 : M), ymm(:, 2 : M) );
    end  
  
    %  Compute the MoC-SE (simple Euler) solution
    %  (the notation "o1" stands for "order 1").
    %  It is denoted as y{pp,mm}_(1) in the paper. 
    ypp_o1(:, 2 : M) = ypp(:, 1:M-1) + dt*k1_pp(:, 1:M-1);
    ymm_o1(:, 1:M-1) = ymm(:, 2 : M) + dt*k1_mm(:, 2 : M);
        
    ypp_o1(:, 1) = BCleft;  
    ymm_o1(:, M) = BCright;
    
    % Compute the MoC-ME (modified Euler) solution
    % (the notation "o2" stands for "order 2").
    % It is denoted as y{pp,mm}_(2) in the paper. 
    for   jNC = 1 : NumComp
        ypp_o2(jNC, 2 : M) = (1/2)*( ypp(jNC, 1:M-1) + ypp_o1(jNC, 2 : M) + ...
                                     dt*fpp{jNC}( ypp_o1(:, 2 : M), ymm_o1(:, 2 : M) ) );                    
        ymm_o2(jNC, 1:M-1) = (1/2)*( ymm(jNC, 2 : M) + ymm_o1(jNC, 1:M-1) + ...
                                     dt*fmm{jNC}( ypp_o1(:, 1:M-1), ymm_o1(:, 1:M-1) ) );
    end

    ypp_o2(:, 1) = BCleft;  
    ymm_o2(:, M) = BCright;
                                                   
    
    % Comute the solution at the half-node.
	% For ypp,  matlab_index = 2 corresponds to grid_index = 1 + 1/2;
	% For ymm,  matlab_index = 1 corresponds to grid_index = 1 + 1/2 = 2 - 1/2.
    ypp_H(:, 2 : M) = (3/4)*ypp(:, 1:M-1) + (1/4)*ypp_o2(:, 2 : M) + (dt/4)*k1_pp(:, 1:M-1);
    ymm_H(:, 1:M-1) = (3/4)*ymm(:, 2 : M) + (1/4)*ymm_o2(:, 1:M-1) + (dt/4)*k1_mm(:, 2 : M);
    
	% Boundary values ypp(:,1) and ymm(:,M) are not needed. 

    
    % The following variables are used by the stage derivative k2 (computed later):
    ypp_fork2(:, 1:M-1) = ypp(:, 1:M-1) + (dt/2)*k1_pp(:, 1:M-1);
    ymm_fork2(:, 2 : M) = ymm(:, 2 : M) + (dt/2)*k1_mm(:, 2 : M);

    
    % Compute k2:
    for   jNC = 1 : NumComp
        k2_pp(jNC, 1:M-1) = fpp{jNC}( ypp_fork2(:, 1:M-1), ymm_H(:,     1:M-1) );
        k2_mm(jNC, 2 : M) = fmm{jNC}( ypp_H(:,     2 : M), ymm_fork2(:, 2 : M) );
    end
    
    % The following variables are used in the stage derivative k3:
    ypp_fork3(:, 1:M-1) = ypp(:, 1:M-1) - dt*k1_pp(:, 1:M-1) + 2*dt*k2_pp(:, 1:M-1);
    ymm_fork3(:, 2 : M) = ymm(:, 2 : M) - dt*k1_mm(:, 2 : M) + 2*dt*k2_mm(:, 2 : M);

    % Compute k3:
    for   jNC = 1 : NumComp
        k3_pp(jNC, 1:M-1) = fpp{jNC}( ypp_fork3(:, 1:M-1), ymm_o2(:,    2 : M) );
        k3_mm(jNC, 2 : M) = fmm{jNC}( ypp_o2(:,    1:M-1), ymm_fork3(:, 2 : M) );
    end
    

    % The RK3 solution:
    ypp_n(:, 2 : M) = ypp(:, 1:M-1) + ...
                      (dt/6)*( k1_pp(:, 1:M-1) + 4*k2_pp(:, 1:M-1) + k3_pp(:, 1:M-1) );
    ymm_n(:, 1:M-1) = ymm(:, 2 : M) + ...
                      (dt/6)*( k1_mm(:, 2 : M) + 4*k2_mm(:, 2 : M) + k3_mm(:, 2 : M) );
    
    ypp_n(:, 1) = BCleft;  
    ymm_n(:, M) = BCright;

    

	


		
		
	
