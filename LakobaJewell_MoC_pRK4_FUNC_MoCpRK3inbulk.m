% This FUNCTION computes the MoC-pRK3 solution
% (see Sec. 6.1 in Lakoba-Jewell paper).
% Inputs for this function are:
%  - the solution at the current time level;
%  - stage derivatives  k1  computed at the current and previous levels in the MAIN code;
%  - stage derivative  k2  computed at the previous level in the MAIN code; and
%  - the boundary conditions at the new level. 
	
function  [ypp_n, ymm_n] = LakobaJewell_MoC_pRK4_FUNC_MoCpRK3inbulk( ...
                           M, dt, NumComp, fpp, fmm, ypp, ymm, BCleft, BCright, ...
                           k1_pp_nmin1, k1_mm_nmin1, ...
                           k2_pp_nmin1, k2_mm_nmin1, ...
                           k1_pp_nmin0, k1_mm_nmin0 );
    
    
    % a)   Compute the solution at the next time level by MoC-ME. 
  
    %      MoC-SE solution:
    ypp_o1(:, 2 : M) = ypp(:, 1:M-1) + dt*k1_pp_nmin0(:, 3:M+1);
    ymm_o1(:, 1:M-1) = ymm(:, 2 : M) + dt*k1_mm_nmin0(:, 1:M-1);
    
    ypp_o1(:, 1) = BCleft;
    ymm_o1(:, M) = BCright;
    
    %      MoC-ME solution:
    for   jNC = 1 : NumComp
        ypp_o2(jNC, 2 : M) = (1/2)*( ypp(jNC, 1:M-1) + ypp_o1(jNC, 2 : M) + ...
	                                 dt*fpp{jNC}( ypp_o1(:, 2 : M), ymm_o1(:, 2 : M) ) );
        ymm_o2(jNC, 1:M-1) = (1/2)*( ymm(jNC, 2 : M) + ymm_o1(jNC, 1:M-1) + ...
                                     dt*fmm{jNC}( ypp_o1(:, 1:M-1), ymm_o1(:, 1:M-1) ) ); 
    end
    
    ypp_o2(:, 1) = BCleft;
    ymm_o2(:, M) = BCright;

    
    % b)   Compute k2^(n-0) for the pRK3, which we denote ell2 
    %      (so as not to confuse it with the k2 of the pRK4 scheme).     
    for   jNC = 1 : NumComp
        ell2_pp(jNC, 3:M+1) = fpp{jNC}( ypp_o1(:, 2 : M), ymm_o2(:, 2 : M) );
        ell2_mm(jNC, 1:M-1) = fmm{jNC}( ypp_o2(:, 1:M-1), ymm_o1(:, 1:M-1) );
    end
    
    % NOTE:  One can use  k2  from the MAIN code (i.e., pRK4) in the pRK3 calculation;
    %        however, one *cannot* use  k2  (i.e., ell2) from the above pRK3 calculation
    %        to compute the pRK4 solution in the MAIN code.
    
    % c) Compute the MoC-pRK3 solution:    
    ypp_n(:, 2 : M) = ypp(:, 1:M-1) + (dt/12)*( 13*k1_pp_nmin0(:, 3:M+1) + 5*ell2_pp(:, 3 : M + 1) - ...
                                                 1*k1_pp_nmin1(:, 2 : M) - 5*k2_pp_nmin1(:, 2 : M) );
    ymm_n(:, 1:M-1) = ymm(:, 2 : M) + (dt/12)*( 13*k1_mm_nmin0(:, 1:M-1) + 5*ell2_mm(:, 1 : M - 1) - ...
                                                 1*k1_mm_nmin1(:, 2 : M) - 5*k2_mm_nmin1(:, 2 : M) );
                                  
    ypp_n(:, 1) = BCleft;    
    ymm_n(:, M) = BCright;   

    
    
    
   
                                 
                                 