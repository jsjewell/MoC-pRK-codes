% This FUNCTION computes the solution by the rotated MoC-pRK3 at a virtual node,
% as described in Sec. 7.2 of the Lakoba-Jewell paper.

% Inputs correspond to figure 11 in the paper.
%
% NOTATIONS:
% whichboundary  = 1 (number "one", not letter "ell", at LEFT), 
%                = 2 (at RIGHT).
% lev{X}p{y}     - stands for p(oint) Y at lev(el) X, as described below;
% lev1p{1,2}     - points 1, 2 in level 1 ( LEFT: m = 1; tn = 4, 2;  RIGHT: m = M, tn = 2, 4),
%                  where values of  tn  are referenced for the calculation at time level n = 3
%                  (they generalize straightforwardly for the general n > 3);
% lev2p{1,2}     - points 1, 2 in level 2 ( LEFT: m = 2; tn = 5, 1;  RIGHT: m = M-1, tn = 1, 5).
%                  Note: 1) lev1 is closer to the computed point than lev2. 
%                        2) The order of points within a level is chosen so that components 
%                           along the "positive" characteristic are listed first. 
%
% "v" in the output refers to the virtual node. 

function [vpp, vmm] = LakobaJewell_MoC_pRK4_FUNC_rotatedMoCpRK3( ...
                      dt, NumComp, fpp, fmm, whichboundary, ...
                      lev1p1_pp, lev1p1_mm, ...
                      lev1p2_pp, lev1p2_mm, ...    
                      lev2p1_pp, lev2p1_mm, ...
                      lev2p2_pp, lev2p2_mm  );

    if whichboundary == 1      %'LEFT' 

        rotsign = -1;     % stencil is rotated by minus 90 degrees
        
    elseif  whichboundary == 2 %'RIGHT'
        
        rotsign = 1;      % stencil is rotated by plus 90 degrees
        
    end
        
    % 1pp)   Compute k1(+) at lev2p1 (this is k1^(n-1)) and ypp_o1 at lev1p1:

    for   jNC = 1 : NumComp
        k1_pp_nmin1(jNC, 1) = fpp{jNC}(lev2p1_pp, lev2p1_mm);   
    end
    
    ypp_o1 = lev2p1_pp + rotsign*dt*k1_pp_nmin1;

    % 1mm)   Compute k1(-) at lev2p2 (this is k1^(n-1)) and ymm_o1 at lev1p2:

    for   jNC = 1 : NumComp
        k1_mm_nmin1(jNC, 1) = fmm{jNC}(lev2p2_pp, lev2p2_mm);
    end
    
    ymm_o1 = lev2p2_mm - rotsign*dt*k1_mm_nmin1;
    
    % 2pp)   Compute k2(+) at lev2p1 (this is k2^(n-1)):
    
    for   jNC = 1 : NumComp
        k2_pp_nmin1(jNC, 1) = fpp{jNC}(ypp_o1,    lev1p1_mm);
    end    
    
    % 2mm)   Compute k2(-) at lev2p2 (this is k2^(n-1)):
    
    for   jNC = 1 : NumComp
        k2_mm_nmin1(jNC, 1) = fmm{jNC}(lev1p2_pp, ymm_o1);
    end
    
    % 3pp)   Compute k1(+) at lev1p1 (this is k1^(n-0)) and ypp_o1 at "v":

    for   jNC = 1 : NumComp
        k1_pp_nmin0(jNC, 1) = fpp{jNC}(lev1p1_pp, lev1p1_mm);
    end
    
    ypp_o1 = lev1p1_pp + rotsign*dt*k1_pp_nmin0;  % we recycle the notation ypp_o1
    
    % 3mm)   Compute k1(-) at lev1p2 (this is k1^(n-0)) and ymm_o1 at "v":

    for   jNC = 1 : NumComp
        k1_mm_nmin0(jNC, 1) = fmm{jNC}(lev1p2_pp, lev1p2_mm);
    end
    
    ymm_o1 = lev1p2_mm - rotsign*dt*k1_mm_nmin0;
    
    % 4)   Compute the MoC-ME solution (order 2, i.e. local error O(h^3)) at "v":
    
    for   jNC = 1 : NumComp
        ypp_o2(jNC, 1) = (1/2)*( lev1p1_pp(jNC, 1) + ypp_o1(jNC, 1) + ...
                                 rotsign*dt*fpp{jNC}(ypp_o1, ymm_o1) );
    
        ymm_o2(jNC, 1) = (1/2)*( lev1p2_mm(jNC, 1) + ymm_o1(jNC, 1) - ...
                                 rotsign*dt*fmm{jNC}(ypp_o1, ymm_o1) );
    end
    
    % 5pp)   Compute k2(+) at lev1p1 (this is k2^(n-0)):
    
    for   jNC = 1 : NumComp
        k2_pp_nmin0(jNC, 1) = fpp{jNC}(ypp_o1, ymm_o2);
    end
    
    % 5mm)   Compute k2(-) at lev1p2 (this is k2^(n-0)):
    
    for   jNC = 1 : NumComp
        k2_mm_nmin0(jNC, 1) = fmm{jNC}(ypp_o2, ymm_o1);
    end
    
    % 6)   Compute the pRK3 solution at "v":
    
    vpp = lev1p1_pp + rotsign*(dt/12)*( ...
                      13*k1_pp_nmin0 + 5*k2_pp_nmin0 - 1*k1_pp_nmin1 - 5*k2_pp_nmin1 );
    vmm = lev1p2_mm - rotsign*(dt/12)*( ...
                      13*k1_mm_nmin0 + 5*k2_mm_nmin0 - 1*k1_mm_nmin1 - 5*k2_mm_nmin1 );
    

