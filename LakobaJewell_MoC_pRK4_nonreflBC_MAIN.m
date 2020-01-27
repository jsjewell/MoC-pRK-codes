% 
% This code implements the MoC-pRK4 scheme with *nonreflecting* boundary conditions.
% The equations are:
%       ypp_t + ymm_x = fpp(ypp, ymm)
%       ymm_t - ymm_x = fmm(ypp, ymm)
% where:
% variables ypp, ymm can have an arbitrary (but equal) number of components, NumComp.
% 
% The code plots the Fourier spectrum of the solution at times specified in the
% "Auxiliary computational parameters" block (see "The user can set:" below).
%
%     ---  This particular code uses:  --- 
% NumComp = 2;
% fpp = (P++) * ypp  +  (P+-) * ymm
% fmm = (P-+) * ypp  +  (P--) * ymm
% where (P++) etc are (2 x 2) matrices. 
% The entire (4 x 4) matrix P = [(P++) (P+-);  (P-+) (P--)] equals matrix P1 of the paper.
% Initial conditions = white noise of magnitude 10^-10 in x-space;
% Homogeneous boundary conditions.
% 
%     ---  The user can set:  ---
%   - Auxiliary computational parameters,
%     such as domain length, number of points, max time, etc.;
%   - Functions  fpp, fmm  (which include parts of matrix P);
%   - Initial Conditions;
%   - Boundary conditions (BC).
%   All of these, except BC, are defined in blocks of code that are respectively labeled.
%   The BC are set at the beginning of the calculations for time levels n = 1, 2, 3 and
%   also at the beginning of the Main Loop over n, which starts at n = 4.
%
%     ---   This MAIN code uses auxiliary FUNCtion codes:  --- 
%   - LakobaJewell_MoC_pRK4_FUNC_MoCRK3 - 
%     for computing the solution at levels n = 2 and 3 by the MoC-RK3 scheme
%     (described in Sec. 3.1 of Lakoba-Jewell paper);
%   - LakobaJewell_MoC_pRK4_FUNC_rotatedMoCpRK3 - 
%     for computing solution at virtual nodes by the rotated MoC-pRK3
%     (described in Sec. 7.2 of Lakoba-Jewell paper);
%   - LakobaJewell_MoC_pRK4_FUNC_MoCpRK3inbulk - 
%     for computing the (auxiliary) solution inside the Main Loop by the MoC-pRK3
%     (described in Sec. 6.1 of Lakoba-Jewell paper).
%
%     ---   Note about code performance:   ---
% This code is optimized for *clarity of organization*, but *not for speed*. 
% That is, it:
%    1) emphasizes the fact that the number of components (NumComp) in y{pp,mm}
%       can be arbitrary (but equal for pp and mm), and
%    2) uses the most compact form of presenting array operations.
%       For example, it would use:
%           (i)   ypp_new = ypp + ymm  
%       rather than
%           (ii)  for j = 1:NumComp;  ypp_new(j,:) = ypp(j,:) + ymm(j,:);  end
%       even though (ii) is about 30-40% (for NumComp=2) faster than (i). 
%       Incidentally, a way to further speed up (i) would be (examplified for NumComp=2):
%           (iii) yppA = ypp(1,:); yppB = ypp(2,:); ymmA = ymm(1,:); ymmB = ymm(2,:);
%                 yppA_new = yppA + ymmA; yppB_new = yppB + ymmB;
%       As far as this particular code is concerned, using way (iii) instead of way (i)
%       would speed it up by a factor around 2. However, it would also make it less clear
%       and more difficult to change should one use it later for a different NumComp. 

 
clear all;

% *************************************************
% Auxiliary computational parameters:
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%
%   x-space:
L = 200;        % length of the computational domain (actually, the length is L-h)
h = .05;        % dx
xmax = L/2;     xmin = -xmax + h;
x = [xmin : h : xmax];
M = length(x);

%   time:
tmin = 0;       tmax = 400;   % maximum simulation time
dt = h;         % by design, the mesh size in x and t is the same for the MoC
t = [tmin : dt : tmax];

tplot = 100;               % time interval at which solution's spectrum is plotted on screen
nplot = round(tplot/dt);   % number of steps after which the plotting is done

%   Fourier space:
%   The wavenumber vector, defined below, is used only to plot the Fourier spectrum of the solution. 
dkk = 2*pi/L;
kk = [0 : (M/2 -1)   -M/2 : -1]*dkk;   % wavenumber is named kk, not k, in order to avoid confusion 
                                       % with stage derivatives in the pRK scheme

spatial_window = exp(-(x/(L/3)).^8);   % Since this code uses nonperiodic BC, the Fourier spectrum of 
                                       % the solution will have high-freqnecy artifacts. To avoid this,
                                       % one multiplies the solution by  spatial_window (only before
                                       % making plots of the spectrum, not in the calculations).
%                                       
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     



% *************************************************                                       
% Conventions for naming variables and other quantities
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

%
%  1) Variables propagating along the positive or negative characteristic 
%     are denoted as  ypp or ymm, respectively. 
%     Here "pp" ("mm") stands for "plus" ("minus").
%     Double letters are used since they are more conspicuous than single ones (i.e., "pp" vs "p"). 
%     The same naming convention is applied to k1, k2, the r.h.s. functions, 
%     and all other related quantities.
%  The code allows any number of variables (denoted as NumComp) per characteristic.
%  This particular version has NumComp = 2. 
%  The only place where one needs to change anything manually is when
%  defining the rhs functions  fpp,  fmm  and the boundary conditions.  
%                                       
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                       


% *************************************************                                       
% RHSs of the equations (matrix P and fpp, fmm)
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

NumComp = 2;                 % Number of components along each chracteristic

P = [ 0   0   0  -1; ...
     1/2  0  1/2  0; ...
      0  -1   0   0; ...
     1/2  0  1/2  0];

fpp{1} = @(spp , smm) ...
          P(1, 1:2)*[spp(1,:); spp(2,:)] + P(1, 3:4)*[smm(1,:); smm(2,:)];
fpp{2} = @(spp , smm) ...
          P(2, 1:2)*[spp(1,:); spp(2,:)] + P(2, 3:4)*[smm(1,:); smm(2,:)];
fmm{1} = @(spp , smm) ...
          P(3, 1:2)*[spp(1,:); spp(2,:)] + P(3, 3:4)*[smm(1,:); smm(2,:)];
fmm{2} = @(spp , smm) ...
          P(4, 1:2)*[spp(1,:); spp(2,:)] + P(4, 3:4)*[smm(1,:); smm(2,:)];
%      
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                       


% *************************************************                                       
% Initial conditions
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

for  jNC = 1 : NumComp
    ypp_0(jNC,:) = zeros(1, M); 
    ymm_0(jNC,:) = zeros(1, M); 
end
% Thus, ypp, ymm are arrays of size  NumComp x M.

noise_ampl = 10^(-10);       % amplitude of noise added to the initial condition
randn('state',0);            % seed of the random number generator, for reproducibility of results
for  jNC = 1 : NumComp
    ypp(jNC,:) = ypp_0(jNC,:) + noise_ampl*randn(size(x));
    ymm(jNC,:) = ymm_0(jNC,:) + noise_ampl*randn(size(x));
end
%      
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                       



% *************************************************                                       
% The MoC-pRK4 computation begins here ----------------------------------- 
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

% Conventions:
% 1) Label "n" denotes the time level *at* which the solution is being computed.
% 2) k{pp,mm}1,2 are stage derivatives of the pRK4 solver, denoted as kappa^{+,-}_1,2 in the paper.
% 3) Their m-indices are defined as follows:
%    kpp(1 : M+1) in the code correspond to kpp(-1 : M-1) in the paper's notations;
%    kmm(1 : M+1) in the code correspond to kmm(2 : M+2) in the paper's notations.
%    This should be contrasted with the indices for y{pp,mm}, which are the same as in the paper
%    (i.e., from 1 to M). 


% n = 1 -------------------------------------------------------------------

% IG-a-1) Solution at this level is given by the initial condition

% IG-a-2) Preallocated an array of a specific size for the boundary conditions:
BCleft =  zeros(NumComp, 1);
BCright = zeros(NumComp, 1);


% IG-b) Compute stage derivative k1 at this level. In the MoC-pRK4 algorithm, 
%       we need stage derivatives at three levels, i.e. levels "n minus 2", "n minus 1" 
%       and "n minus 0".
%       The k1 computed at the current level (n = 1) will first be used in the Main Loop 
%       to find solution at n = 4. Thus, it corresponds to k1^(n-2) in the algorithm;
%       hence we denote it  k1_nmin2. The one at n = 2 will be k1^(n-1) (k1_nmin1), and
%       the one at n = 3 will be k1^n  (k1_nmin0); they all will be used to compute the 
%       MoC-pRK4 solution at n = 4. 
%   
%       Leave 2 entries for k1_pp on the left and for k1_mm on the right, 
%       to be supplied at n = 3. 

for  jNC = 1 : NumComp
    k1_pp_nmin2(jNC, 3:M+1) = fpp{jNC}( ypp(:, 1:M-1), ymm(:, 1:M-1) );
    k1_mm_nmin2(jNC, 1:M-1) = fmm{jNC}( ypp(:, 2 : M), ymm(:, 2 : M) );
end


% IG-c) Compute SE (simple Euler) solutions at level  n+1 ( = 2), which will be used 
%       at the next level for the computation of k2_nmin2. 
%       The notation "o1" stands for "order 1".
%       This solution is denoted as y{pp,mm}_(1) in the paper. 


ypp_o1(:, 2 : M) = ypp(:, 1:M-1) + dt*k1_pp_nmin2(:, 3:M+1);
ymm_o1(:, 1:M-1) = ymm(:, 2 : M) + dt*k1_mm_nmin2(:, 1:M-1);

% The boundary values of ypp and ymm appear *not* be be needed at the next level.
% (Moreover, it would not have been correct to assign them to equal the exact solutions
%  because it is critical that the o1-variables used in k2 be computed by the SE and not 
%  be exact! (This wouldn't the case for the o2- or o3-variables; those, when they are needed,
%  can be exact or higher-order accurate.) )



% n = 2 -------------------------------------------------------------------

% Set the BC at n = 2

% Set the nonreflecting boundary conditions:

for   jNC = 1 : NumComp
    BCleft(jNC) = 0;
    BCright(jNC) = 0;
end

% Do calculations OUTSIDE the grid first. ..........................................

% OG-a) Extrapolate the values at nodes (m,n) = (0,1) and (M+1,1). 
%       We denote the values at the "virtual" nodes by "v". 
%       They will be first used at the step n = 3.
%       They will be (NumComp x 2) arrays, 
%       where the "2" is the number of the virtual nodes added on each side of the grid. 
%       Thus, 
%       v{pp,mm}_left(:,  1 or 2) means  y{pp,mm}(:, m =  0  or -1 ),
%       v{pp,mm}_right(:, 1 or 2) means  y{pp,mm}(:, m = M+1 or M+2).
%       - NOTES:  1. Notations v{pp,mm} will be recycled at every time level.
%                 2. Note that these extrapolations obtain v-values at the 
%                    *previous* level (which right now is n = 1). 

vpp_left(:, 1) = 4*ypp(:, 1) - 6*ypp(:, 2) + 4*ypp(:, 3) - ypp(:, 4);
vmm_left(:, 1) = 4*ymm(:, 1) - 6*ymm(:, 2) + 4*ymm(:, 3) - ymm(:, 4);

vpp_right(:, 1) = 4*ypp(:, M) - 6*ypp(:, M-1) + 4*ypp(:, M-2) - ypp(:, M-3);
vmm_right(:, 1) = 4*ymm(:, M) - 6*ymm(:, M-1) + 4*ymm(:, M-2) - ymm(:, M-3);


% OG-b) Extrapolate to m = -1, M+2:

vpp_left(:, 2) = 4*vpp_left(:, 1) - 6*ypp(:, 1) + 4*ypp(:, 2) - ypp(:, 3);
vmm_left(:, 2) = 4*vmm_left(:, 1) - 6*ymm(:, 1) + 4*ymm(:, 2) - ymm(:, 3);

vpp_right(:, 2) = 4*vpp_right(:, 1) - 6*ypp(:, M) + 4*ypp(:, M-1) - ypp(:, M-2);
vmm_right(:, 2) = 4*vmm_right(:, 1) - 6*ymm(:, M) + 4*ymm(:, M-1) - ymm(:, M-2);


% OG-c) Store boundary and near-boundary values (denoted as BV). 
% 		At the left boundary, the values are stored in the usual order, 
% 		i.e., from left (farthest from the boundary on the left) to the right.
% 		At the right boundary, the values are stored in the *reverse* order, 
% 		i.e., from right (farthest from the boundary on the right) to the left.
% 		These values are stored in arrays of size  ( NumComp x 5 x 4),  
% 		where the 5 is number of time levels and 4 is the number of stored nodes per level. 
%       I.e., at the first 5 time levels, the indices mean the following:
%       BV_left(  {ind. of field's component}, {time level n = 1 : 5}, {m = -1,   0,  1, 2  } ),
%       BV_right( {ind. of field's component}, {time level n = 1 : 5}, {m = M+2, M+1, M, M-1} ).
%       At subsequent time levels, 4 levels (2nd index) of the BV values will be 
%       reassigned from old ones and one new level will be computed.

BVpp_left(:, 1, 1:2) = vpp_left(:, [2,1]); 
BVmm_left(:, 1, 1:2) = vmm_left(:, [2,1]); 

BVpp_left(:, 1, 3:4) = ypp(:, [1,2]); 
BVmm_left(:, 1, 3:4) = ymm(:, [1,2]); 

BVpp_right(:, 1, 1:2) = vpp_right(:, [2,1]); 
BVmm_right(:, 1, 1:2) = vmm_right(:, [2,1]);  

BVpp_right(:, 1, 3:4) = ypp(:, [M,M-1]); 
BVmm_right(:, 1, 3:4) = ymm(:, [M,M-1]); 
 
 
% Do calculations INSIDE the grid. ..........................................

% IG-a)    MoC-RK3 solution:
    
[ypp, ymm] = LakobaJewell_MoC_pRK4_FUNC_MoCRK3( ...
             M, dt, NumComp, fpp, fmm, ypp, ymm, BCleft, BCright );


% IG-b)   Compute stage derivative k2 at n = 1 using the MoC-SE solutions 
%         computed at the end of the n = 1 step.
%         (This is new compared to the n = 1 step, but will be repeated at the general-n step.)

for   jNC = 1 : NumComp
    k2_pp_nmin2(jNC, 3:M+1) = fpp{jNC}( ypp_o1(:, 2 : M), ymm(:,    2 : M) );
    k2_mm_nmin2(jNC, 1:M-1) = fmm{jNC}( ypp(:,    1:M-1), ymm_o1(:, 1:M-1) );
end


% IG-c)   Compute k1. It will first be used in the Main Loop to find solution at n = 4,
%         so denote it k1^(n-1) (the one at n = 3 will be k1^n). 
%     
%         Leave 2 entries for k1_pp on the left and for k1_mm on the right, 
%         to be supplied at n = 4. 

for   jNC = 1 : NumComp
    k1_pp_nmin1(jNC, 3:M+1) = fpp{jNC}( ypp(:, 1:M-1), ymm(:, 1:M-1) );
    k1_mm_nmin1(jNC, 1:M-1) = fmm{jNC}( ypp(:, 2 : M), ymm(:, 2 : M) );
end

% IG-d)   Compute the MoC-SE solutions at level  n+1 ( = 3), which will be used 
%         at the next level for the computation of k2 at level n = 2. 

ypp_o1(:, 2 : M) = ypp(:, 1:M-1) + dt*k1_pp_nmin1(:, 3:M+1);
ymm_o1(:, 1:M-1) = ymm(:, 2 : M) + dt*k1_mm_nmin1(:, 1:M-1);

% Again, as for n = 1, boundary values to the MoC-SE solution are not needed.



% n = 3 -------------------------------------------------------------------

% Set the BC at n = 3

for   jNC = 1 : NumComp
    BCleft(jNC) = 0;
    BCright(jNC) = 0;
end


% Do calculations OUTSIDE the grid first. ..........................................

% OG-a)   Extrapolate to m = 0, M+1  (same as for n = 2) 
%         Note again that these extrapolations obtain v-values at the 
%         *previous* level, n = 2. 


vpp_left(:, 1) = 4*ypp(:, 1) - 6*ypp(:, 2) + 4*ypp(:, 3) - ypp(:, 4);
vmm_left(:, 1) = 4*ymm(:, 1) - 6*ymm(:, 2) + 4*ymm(:, 3) - ymm(:, 4);

vpp_right(:, 1) = 4*ypp(:, M) - 6*ypp(:, M-1) + 4*ypp(:, M-2) - ypp(:, M-3);
vmm_right(:, 1) = 4*ymm(:, M) - 6*ymm(:, M-1) + 4*ymm(:, M-2) - ymm(:, M-3);


% OG-b)   Extrapolate to m = -1, M+2  (same as for n = 2):

vpp_left(:, 2) = 4*vpp_left(:, 1) - 6*ypp(:, 1) + 4*ypp(:, 2) - ypp(:, 3);
vmm_left(:, 2) = 4*vmm_left(:, 1) - 6*ymm(:, 1) + 4*ymm(:, 2) - ymm(:, 3);

vpp_right(:, 2) = 4*vpp_right(:, 1) - 6*ypp(:, M) + 4*ypp(:, M-1) - ypp(:, M-2);
vmm_right(:, 2) = 4*vmm_right(:, 1) - 6*ymm(:, M) + 4*ymm(:, M-1) - ymm(:, M-2);


% OG-c)   Store boundary and near-boundary values  (same as for n = 2, 
%         but at the next level) 
BVpp_left(:, 2, 1:2) = vpp_left(:, [2,1]); 
BVmm_left(:, 2, 1:2) = vmm_left(:, [2,1]); 

BVpp_left(:, 2, 3:4) = ypp(:, [1,2]); 
BVmm_left(:, 2, 3:4) = ymm(:, [1,2]); 

BVpp_right(:, 2, 1:2) = vpp_right(:, [2,1]); 
BVmm_right(:, 2, 1:2) = vmm_right(:, [2,1]);  

BVpp_right(:, 2, 3:4) = ypp(:, [M,M-1]); 
BVmm_right(:, 2, 3:4) = ymm(:, [M,M-1]); 


% OG-d)   Assign a previously computed value to the 1st (= closest to the grid) 
%         virtual node at the *previous* level (n-2) ( = 1 here):

vpp_left(:, 1) = BVpp_left(:, 1, 2);
vmm_left(:, 1) = BVmm_left(:, 1, 2);   

vpp_right(:, 1) = BVpp_right(:, 1, 2);
vmm_right(:, 1) = BVmm_right(:, 1, 2);

% OG-e)   Compute k1 at level (n-2) ( = 1 here) at the 1st virtual node:

for   jNC = 1 : NumComp
    k1_pp_nmin2(:, 2) = fpp{jNC}( vpp_left(:,  1), vmm_left(:,  1) );
    k1_mm_nmin2(:, M) = fmm{jNC}( vpp_right(:, 1), vmm_right(:, 1) );
end


% OG-f)  Compute k2 at level (n-2) at the 1st virtual node: 

%        Compute MoC-SE solution at level  n-1 ( = 2 here), which will be used 
%        to compute k2^(n-2). 

vpp_o1_left =  vpp_left(:,  1) + dt*k1_pp_nmin2(:, 2);
vmm_o1_right = vmm_right(:, 1) + dt*k1_mm_nmin2(:, M);

% Compute k2

for   jNC = 1 : NumComp
    k2_pp_nmin2(jNC, 2) = fpp{jNC}( vpp_o1_left, ymm(:, 1) );
    k2_mm_nmin2(jNC, M) = fmm{jNC}( ypp(:, M),   vmm_o1_right );
end



% Do calculations INSIDE the grid.  ........................................................

% IG-a)   MoC-RK3 solution:
    
[ypp, ymm] = LakobaJewell_MoC_pRK4_FUNC_MoCRK3( ...
             M, dt, NumComp, fpp, fmm, ypp, ymm, BCleft, BCright );


% IG-b)   Store boundary and near-boundary values  (same as for n = 2, but only INSIDE the grid): 

BVpp_left(:, 3, 3:4) = ypp(:, [1,2]); 
BVmm_left(:, 3, 3:4) = ymm(:, [1,2]); 

BVpp_right(:, 3, 3:4) = ypp(:, [M,M-1]); 
BVmm_right(:, 3, 3:4) = ymm(:, [M,M-1]); 


% IG-c)   Compute k2 at n = 2 using the MoC-SE solution computed at the end of the n = 2 step
%         (same as for the n = 2 step and for the general-n step):

for   jNC = 1 : NumComp
    k2_pp_nmin1(jNC, 3:M+1) = fpp{jNC}( ypp_o1(:, 2 : M), ymm(:,    2 : M) );
    k2_mm_nmin1(jNC, 1:M-1) = fmm{jNC}( ypp(:,    1:M-1), ymm_o1(:, 1:M-1) );
end


% IG-d)   Compute k1. It will first be used in the Main Loop to find solution at the next level, 
%         i.e., n = 4;  so denote it k1^(n-0). 
%   
%         Leave 2 entries for k1_pp on the left and for k1_mm on the right, 
%         to be supplied at n = 4 and 5. 

for   jNC = 1 : NumComp
    k1_pp_nmin0(jNC, 3:M+1) = fpp{jNC}( ypp(:, 1:M-1), ymm(:, 1:M-1) );
    k1_mm_nmin0(jNC, 1:M-1) = fmm{jNC}( ypp(:, 2 : M), ymm(:, 2 : M) );
end


% IG-e)   Compute MoC-SE solution at level  n+1 ( = 4 here), which will be used at the 
%         next level for the computation of k2^n ( n = 3 here) (same as the general-n step): 

ypp_o1(:, 2 : M) = ypp(:, 1:M-1) + dt*k1_pp_nmin0(:, 3:M+1);
ymm_o1(:, 1:M-1) = ymm(:, 2 : M) + dt*k1_mm_nmin0(:, 1:M-1);

% Again, the boundary values of the o1-variables are not needed.



% The Main Loop begins here:   -----------------------------------------------------

for n = 4 : length(t)-1    
     
   % Set the BC at level n:
   
    for   jNC = 1 : NumComp
        BCleft(jNC) = 0;
        BCright(jNC) = 0;
    end
    
    tic
    
    % Do calculations OUTSIDE the grid first. ..........................................

    % OG-a)   Compute the low-order solution at m = 2 and M-1 at level  n:
    if  n > 4 
       
        %     Compute solution only at m = 2 and M-1  by MoC-ME (local error O(h^3)):
        
        for   jNC = 1 : NumComp
            ypp_left_o2(jNC, 1)  = (1/2)*( ypp(jNC, 1) +   ypp_o1(jNC, 2) + ...
                                           dt*fpp{jNC}( ypp_o1(:, 2), ymm_o1(:, 2) ) );
            ymm_left_o2(jNC, 1)  = (1/2)*( ymm(jNC, 3) +   ymm_o1(jNC, 2) + ...
                                           dt*fmm{jNC}( ypp_o1(:, 2), ymm_o1(:, 2) ) );
            
            ypp_right_o2(jNC, 1) = (1/2)*( ypp(jNC, M-2) + ypp_o1(jNC, M-1) + ...
                                           dt*fpp{jNC}( ypp_o1(:, M-1), ymm_o1(:, M-1) ) );
            ymm_right_o2(jNC, 1) = (1/2)*( ymm(jNC, M) +   ymm_o1(jNC, M-1) + ...
                                           dt*fmm{jNC}( ypp_o1(:, M-1), ymm_o1(:, M-1) ) );
        end
        
    end
    
    % OG-b)   Obtain k1 and k2 at the 1st (= closest to the grid) virtual node 
    %         at the *previous* level (n-2).
    %         This consists of several sub-steps.
    
    %     OG-b-1)   Compute or assign (see below) a value at this virtual node:

    if  n == 4  % assign a previously stored value:
        
        vpp_left(:, 1) = BVpp_left(:, 2, 2);
        vmm_left(:, 1) = BVmm_left(:, 2, 2);

        vpp_right(:, 1) = BVpp_right(:, 2, 2);
        vmm_right(:, 1) = BVmm_right(:, 2, 2);
        
    else   % if n >= 5, compute this value by the rotated MoC-pRK3 and store it:
                
        [vpp_left(:, 1),  vmm_left(:, 1)] = ...
                      LakobaJewell_MoC_pRK4_FUNC_rotatedMoCpRK3( ...
                      dt, NumComp, fpp, fmm, 1, ... % "1" stands for LEFT boundary 
                      BVpp_left(:, end,   3), BVmm_left(:, end,   3), ...
                      ... % for n = 5, end = 4;  for n >= 6, end = 5
                      BVpp_left(:, end-2, 3), BVmm_left(:, end-2, 3), ...
                      ypp_left_o2,            ymm_left_o2,           ...
                      BVpp_left(:, end-3, 4), BVmm_left(:, end-3, 4) );
                                    
        [vpp_right(:,1),  vmm_right(:,1)] = ...
                      LakobaJewell_MoC_pRK4_FUNC_rotatedMoCpRK3( ...
                      dt, NumComp, fpp, fmm, 2, ... % "2" stands for RIGHT boundary 
                      BVpp_right(:, end-2, 3), BVpp_right(:, end-2, 3), ...
                      BVpp_right(:, end,   3), BVmm_right(:, end,   3), ...
                      BVpp_right(:, end-3, 4), BVmm_right(:, end-3, 4), ...
                      ypp_right_o2,            ymm_right_o2 );
                                    
        BVpp_left(:, end-1, 2) = vpp_left(:, 1);
        BVmm_left(:, end-1, 2) = vmm_left(:, 1);
        
        BVpp_right(:, end-1, 2) = vpp_right(:, 1);
        BVmm_right(:, end-1, 2) = vmm_right(:, 1);

    end
    

    %     OG-b-2)  Compute k1 at level (n-2) at this virtual node:

    for   jNC = 1 : NumComp
        k1_pp_nmin1(jNC, 2) = fpp{jNC}(vpp_left(:,  1), vmm_left(:,  1));
        k1_mm_nmin1(jNC, M) = fmm{jNC}(vpp_right(:, 1), vmm_right(:, 1));
    end
    
    %     OG-b-3)  Compute k2 at level (n-2) at this virtual node: 

    %              Compute MoC-SE solution at level  n-1, 
    %              which will be used to compute k2^(n-2): 

    vpp_o1_left =  vpp_left(:,  1) + dt*k1_pp_nmin1(:, 2);

    vmm_o1_right = vmm_right(:, 1) + dt*k1_mm_nmin1(:, M);

    %             Compute k2^(n-2):

    for   jNC = 1 : NumComp
        k2_pp_nmin1(jNC, 2) = fpp{jNC}(vpp_o1_left, ymm(:, 1));
        k2_mm_nmin1(jNC, M) = fmm{jNC}(ypp(:, M),   vmm_o1_right);
    end

    
    % OG-c)   Obtain k1 and k2 at the 2nd ( = 2nd-closest to the grid) virtual node 
    %         *two levels ago*, i.e., at (n-3).
    %         This consists of several sub-steps.
    
    %     OG-c-1)   Compute or assign (see below) a value at this virtual node:

    if  n == 4  ||  n == 5  % assign a previously stored value
        
        vpp_left(:, 2) = BVpp_left(:, n-3, 1);
        vmm_left(:, 2) = BVmm_left(:, n-3, 1);

        vpp_right(:, 2) = BVpp_right(:, n-3, 1);
        vmm_right(:, 2) = BVmm_right(:, n-3, 1);
        
    else   % if n > 5, compute this value by the rotated MoC-pRK3 and store it
        
        [vpp_left(:,2),  vmm_left(:,2)] = ...
                      LakobaJewell_MoC_pRK4_FUNC_rotatedMoCpRK3( ...
                      dt, NumComp, fpp, fmm, 1, ... % "1" stands for LEFT boundary 
                      BVpp_left(:, 4, 2), BVmm_left(:, 4, 2), ...
                      BVpp_left(:, 2, 2), BVmm_left(:, 2, 2), ...
                      BVpp_left(:, 5, 3), BVmm_left(:, 5, 3), ...
                      BVpp_left(:, 1, 3), BVmm_left(:, 1, 3) );
                  
        [vpp_right(:,2),  vmm_right(:,2)] = ...
                      LakobaJewell_MoC_pRK4_FUNC_rotatedMoCpRK3( ...
                      dt, NumComp, fpp, fmm, 2, ... % "2" stands for RIGHT boundary 
                      BVpp_right(:, 2, 2), BVmm_right(:, 2, 2), ...
                      BVpp_right(:, 4, 2), BVmm_right(:, 4, 2), ...
                      BVpp_right(:, 1, 3), BVmm_right(:, 1, 3), ...
                      BVpp_right(:, 5, 3), BVmm_right(:, 5, 3) );
                  
        BVpp_left(:, 3, 1) = vpp_left(:, 2);
        BVmm_left(:, 3, 1) = vmm_left(:, 2);
        
        BVpp_right(:, 3, 1) = vpp_right(:, 2);
        BVmm_right(:, 3, 1) = vmm_right(:, 2);
        
    end
    
    %     OG-c-2)  Compute k1 at level (n-3) at this virtual node:

    for   jNC = 1 : NumComp
        k1_pp_nmin2(jNC, 1) =   fpp{jNC}(vpp_left(:,  2), vmm_left(:,  2));
        k1_mm_nmin2(jNC, M+1) = fmm{jNC}(vpp_right(:, 2), vmm_right(:, 2));
    end
    
    %     OG-c-3)  Compute k2 at level (n-3) at this virtual node: 

    %              Compute MoC-SE solution at level  n-2, 
    %              which will be used to compute k2^(n-3): 

    vpp_o1_left =  vpp_left(:,2)  + dt*k1_pp_nmin2(:,1);   
	               % we recycle the notation  v{pp,mm}_o1_{left,right}
    vmm_o1_right = vmm_right(:,2) + dt*k1_mm_nmin2(:,M+1);

    %              Compute k2^(n-3):

    for   jNC = 1 : NumComp
        k2_pp_nmin2(jNC, 1) =   fpp{jNC}(vpp_o1_left,    vmm_left(:,1));
        k2_mm_nmin2(jNC, M+1) = fmm{jNC}(vpp_right(:,1), vmm_o1_right);
    end
    
    
    
    % Do calculations INSIDE the grid.  ........................................................
    % ..........................................................................................

    % IG-a)    Compute the MoC-pRK3 solution at the next time level:
    
    [ypp_o3, ymm_o3] = LakobaJewell_MoC_pRK4_FUNC_MoCpRK3inbulk( ...
                       M, dt, NumComp, fpp, fmm, ypp, ymm, BCleft, BCright, ...
                       k1_pp_nmin1, k1_mm_nmin1, ...
                       k2_pp_nmin1, k2_mm_nmin1, ...
                       k1_pp_nmin0, k1_mm_nmin0 );
    
            
    % IG-b)   Compute k2^(n-0), to be used to compute the solution at this level:
    
    for   jNC = 1 : NumComp
        k2_pp_nmin0(jNC, 3:M+1) = fpp{jNC}( ypp_o1(:, 2 : M), ymm_o3(:, 2 : M) );
        k2_mm_nmin0(jNC, 1:M-1) = fmm{jNC}( ypp_o3(:, 1:M-1), ymm_o1(:, 1:M-1) );
    end
    
    
    % IG-c)   Compute the MoC-pRK4 solution:
    
    ypp_n(:, 2 : M) = ypp(:, 1:M-1) + (dt/24)*( 37*k1_pp_nmin0(:, 3:M+1) +  9*k2_pp_nmin0(:, 3:M+1) - ...
                                                14*k1_pp_nmin1(:, 2 : M) - 18*k2_pp_nmin1(:, 2 : M) + ...
                                                 1*k1_pp_nmin2(:, 1:M-1) +  9*k2_pp_nmin2(:, 1:M-1) );
                                 
    ymm_n(:, 1:M-1) = ymm(:, 2 : M) + (dt/24)*( 37*k1_mm_nmin0(:, 1:M-1) +  9*k2_mm_nmin0(:, 1:M-1) - ...
                                                14*k1_mm_nmin1(:, 2 : M) - 18*k2_mm_nmin1(:, 2 : M) + ...
                                                 1*k1_mm_nmin2(:, 3:M+1) +  9*k2_mm_nmin2(:, 3:M+1) );
                                             
    ypp_n(:, 1) = BCleft;    
    ymm_n(:, M) = BCright;   
    
    %         Update the solution:
    
    ypp = ypp_n; 
    ymm = ymm_n; 


    % IG-d)   Store or (reassign + store) new near-boundary values:
               
    if  n == 4  ||  n == 5   %  <--  store only
        
        BVpp_left(:, end+1, 3:4) = ypp(:, [1,2]);  % for n = 4, end = 3;   for n = 5, end = 4
        BVmm_left(:, end+1, 3:4) = ymm(:, [1,2]); 

        BVpp_right(:, end+1, 3:4) = ypp(:, [M,M-1]); 
        BVmm_right(:, end+1, 3:4) = ymm(:, [M,M-1]); 

    else   % i.e. for n > 5   <--  reassign + store
        
        for jBV = 1 : 4
            
            BVpp_left(:, jBV, :) = BVpp_left(:, jBV+1, :);
            BVmm_left(:, jBV, :) = BVmm_left(:, jBV+1, :);
            
            BVpp_right(:, jBV, :) = BVpp_right(:, jBV+1, :);
            BVmm_right(:, jBV, :) = BVmm_right(:, jBV+1, :);
            
        end
        
        BVpp_left(:, 5, 3:4) = ypp(:, [1,2]);
        BVmm_left(:, 5, 3:4) = ymm(:, [1,2]);
        
        BVpp_right(:, 5, 3:4) = ypp(:, [M,M-1]); 
        BVmm_right(:, 5, 3:4) = ymm(:, [M,M-1]);
        
    end
    
    
    % IG-e)   Reassign k1, k2 from two previous levels:
    
    k1_pp_nmin2(:, 2:M+1) = k1_pp_nmin1(:, 2:M+1);
    k1_mm_nmin2(:, 1 : M) = k1_mm_nmin1(:, 1 : M);
    
    k1_pp_nmin1(:, 3:M+1) = k1_pp_nmin0(:, 3:M+1);
    k1_mm_nmin1(:, 1:M-1) = k1_mm_nmin0(:, 1:M-1);
    
    k2_pp_nmin2(:, 2:M+1) = k2_pp_nmin1(:, 2:M+1);
    k2_mm_nmin2(:, 1 : M) = k2_mm_nmin1(:, 1 : M);
    
    k2_pp_nmin1(:, 3:M+1) = k2_pp_nmin0(:, 3:M+1);
    k2_mm_nmin1(:, 1:M-1) = k2_mm_nmin0(:, 1:M-1);
    
    
    % IG-f)   Compute k1^(n-0):
    
    %         Leave 2 entries for k1_pp on the left and for k1_mm on the right, 
    %         to be supplied at the next two levels. 

    for   jNC = 1 : NumComp
        k1_pp_nmin0(jNC, 3:M+1) = fpp{jNC}( ypp(:, 1:M-1), ymm(:, 1:M-1) );
        k1_mm_nmin0(jNC, 1:M-1) = fmm{jNC}( ypp(:, 2 : M), ymm(:, 2 : M) );
    end
    
    
    % IG-g)   Compute MoC-SE solution at level  n+1, which will be used 
    %         at the next level for the computation of k2^n:

    ypp_o1(:, 2 : M) = ypp(:, 1:M-1) + dt*k1_pp_nmin0(:, 3:M+1);
    ymm_o1(:, 1:M-1) = ymm(:, 2 : M) + dt*k1_mm_nmin0(:, 1:M-1);

    % The boundary values of the o1-variables are not needed.

        
% ---------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------

    % This loop plots the Fourier spectrum of the solution
    if rem(n, nplot) == 0
        total_spectrum_sq = zeros(size(kk));
        for   jNC = 1 : NumComp
            total_spectrum_sq = total_spectrum_sq + ...
                  abs(fft(ypp(jNC, :).*spatial_window)).^2 + abs(fft(ymm(jNC, :).*spatial_window)).^2;
        end
        
         figure(41)
         plot(h*fftshift(kk), fftshift(log10( sqrt(total_spectrum_sq) )), ...
              'k', 'linewidth', 2)           
         set(gca,'fontsize',20);
         xlim([-0.0, pi+0.002]);
         ShowAxisX = [0  pi/2  pi ];
         ShowAxisXlabel = {'0',  '\pi/2',  '\pi'};
         set(gca,'xtick', ShowAxisX,'xticklabel', ShowAxisXlabel,'fontsize',20)
         xlabel('z');
         ylabel('$\log_{10}\| F[\mathbf{y}] \|$','Interpreter','Latex');
         title(['spectrum at t = ',num2str(n*dt)],'fontsize',14)
         pause(0.5)
    end
    
end    % end of the main loop over  n 

