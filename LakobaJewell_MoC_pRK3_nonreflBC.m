% 
% This code implements the MoC-pRK3 scheme with *nonreflecting* boundary conditions.
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
%   The BC are set at the beginning of the calculations for time levels n = 1, 2 and
%   also at the beginning of the Main Loop over n, which starts at n = 3.
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
L = 100;        % length of the computational domain (actually, the length is L-h)
h = .05;        % dx
xmax = L/2;     xmin = -xmax + h;
x = [xmin : h : xmax];
M = length(x);

%   time:
tmin = 0;       tmax = 1000;   % maximum simulation time
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
% The MoC-pRK3 computation begins here ----------------------------------- 
% vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

% Conventions:
% 1) Label "n" denotes the time level *at* which the solution is being computed.
% 2) k{pp,mm}1,2 are stage derivatives of the pRK3 solver, denoted as kappa^{+,-}_1,2 in the paper.
% 3) Their m-indices are defined as follows:
%    kpp(1 : M) in the code correspond to kpp(0 : M-1) in the paper's notations;
%    kmm(1 : M) in the code correspond to kmm(2 : M+1) in the paper's notations.
%    This should be contrasted with the indices for y{pp,mm}, which are the same as in the paper
%    (i.e., from 1 to M). 


% n = 1 -------------------------------------------------------------------

% IG-a-1) Solution at this level is given by the initial condition

% IG-a-2) Preallocated an array of a specific size for the boundary conditions:
BCleft =  zeros(NumComp, 1);
BCright = zeros(NumComp, 1);
    
	
% IG-b) Store boundary values (BV) for later use. These will be arrays of size
%       NumComp x 3,  where the "3" is specific to MoC-pRK3 (as opposed to -pRK4).

BVpp_left(:, 1) = ypp(:, 1);    
BVmm_left(:, 1) = ymm(:, 1);

BVpp_right(:, 1) = ypp(:, M);
BVmm_right(:, 1) = ymm(:, M);


% IG-c) Compute k1 at this level. In the MoC-pRK3 algorithm, 
%    we need stage derivatives at two levels, i.e. levels "n minus 1" and "n minus 0", 
%    hence the "_nmin1" notation now and "_nmin0" notation later. 
%   
%    Leave 1 entry for k1_pp on the left and for k1_mm on the right; see Conventions 3) above.
%    These missing entries will be supplied at n = 3. 

for  jNC = 1 : NumComp
    k1_pp_nmin1(jNC, 2 : M) = fpp{jNC}( ypp(:, 1:M-1), ymm(:, 1:M-1) );
    k1_mm_nmin1(jNC, 1:M-1) = fmm{jNC}( ypp(:, 2 : M), ymm(:, 2 : M) );
end


% n = 2 -------------------------------------------------------------------

% Set the nonreflecting boundary conditions:

for   jNC = 1 : NumComp
    BCleft(jNC) = 0;
    BCright(jNC) = 0;
end


% OG-a) Extrapolate the values at nodes (m,n) = (0,1) and (M+1,1),
%       denoting the values at the "virtual" nodes by "v". 
%       They will be used at the step n = 3.
  
vpp_left = 3*ypp(:, 1) - 3*ypp(:, 2) + ypp(:, 3);
vmm_left = 3*ymm(:, 1) - 3*ymm(:, 2) + ymm(:, 3);

vpp_right = 3*ypp(:, M) - 3*ypp(:, M-1) + ypp(:, M-2);
vmm_right = 3*ymm(:, M) - 3*ymm(:, M-1) + ymm(:, M-2);

 

% IG-a) Compute the solution at n = 2 by MoC-ME (Modified Euler)

% y{pp,mm}_o1 refer to the MoC-SE (Simple Euler) solution ("o1" stands for "order 1").
% It is denoted as y{pp,mm}_(1) in the paper. 

ypp_o1(:, 2 : M) = ypp(:, 1:M-1) + dt*k1_pp_nmin1(:, 2 : M);
ymm_o1(:, 1:M-1) = ymm(:, 2 : M) + dt*k1_mm_nmin1(:, 1:M-1);

% Assign boundary values:
ypp_o1(:, 1) = BCleft;  
ymm_o1(:, M) = BCright;

% y{pp,mm}_o2 refer to the MoC-ME solution ("o2" stands for "order 2").
% It is denoted as y{pp,mm}_(2) in the paper. 

for   jNC = 1 : NumComp
    ypp_o2(jNC, 2 : M) = (1/2)*( ypp(jNC, 1:M-1) + ypp_o1(jNC, 2 : M) + ...
                                 dt*fpp{jNC}( ypp_o1(:, 2 : M), ymm_o1(:, 2 : M) ) );                    
    ymm_o2(jNC, 1:M-1) = (1/2)*( ymm(jNC, 2 : M) + ymm_o1(jNC, 1:M-1) + ...
                                 dt*fmm{jNC}( ypp_o1(:, 1:M-1), ymm_o1(:, 1:M-1) ) );
end

ypp_o2(:, 1) = BCleft;  
ymm_o2(:, M) = BCright;


% IG-b) Reassign the solution and store boundary values:

ypp = ypp_o2;
ymm = ymm_o2;

BVpp_left(:, 2) = ypp(:, 1);    
BVmm_left(:, 2) = ymm(:, 1);

BVpp_right(:, 2) = ypp(:, M);
BVmm_right(:, 2) = ymm(:, M);


% IG-c) Compute k2 at the *previous* level:

for   jNC = 1 : NumComp
    k2_pp_nmin1(jNC, 2 : M) = fpp{jNC}( ypp_o1(:, 2 : M), ymm_o2(:, 2 : M) );
    k2_mm_nmin1(jNC, 1:M-1) = fmm{jNC}( ypp_o2(:, 1:M-1), ymm_o1(:, 1:M-1) );
end

% d) Compute the k1-stage derivatives at the current level (n=2):

for   jNC = 1 : NumComp
    k1_pp_nmin0(jNC, 2 : M) = fpp{jNC}( ypp(:, 1:M-1), ymm(:, 1:M-1) );
    k1_mm_nmin0(jNC, 1:M-1) = fmm{jNC}( ypp(:, 2 : M), ymm(:, 2 : M) );
end


% Main loop: for n >= 3 ----------------------------------------------------------
% --------------------------------------------------------------------------------

for n = 3 : length(t)-1 
    % Recall Convention 1) before the "n = 1" case above:
    % The solution is computed at level  n  using the information at levels (n-1) and below.
         
    % Nonreflecting boundary conditions:
    
    for   jNC = 1 : NumComp
        BCleft(jNC) = 0;
        BCright(jNC) = 0;
    end

    % Do calculations OUTSIDE the grid ........................................
	%
    % OG-a) Compute the solution at the virtual node two level ago (n-2):
    
    if  n == 3  % in this case only, retrieve the value found at n = 1 by extrapolation
        
        % No calculations, since v{pp,mm}_{left,right} have already been found at n = 1.
        
        % Preallocate an array of correct size for y{pp,mm}_o1_{left,right}:
        
        ypp_o1_left = zeros(NumComp, 1);
        ymm_o1_left = zeros(NumComp, 1);
        
        ypp_o1_right = zeros(NumComp, 1);
        ymm_o1_right = zeros(NumComp, 1);
        
    else  % for n >= 4, find the solution at the virtual nodes by the rotated MoC-ME 
           
        % Rotated MoC-SE for the left boundary:

        for   jNC = 1 : NumComp
            ypp_o1_left(jNC) = BVpp_left(jNC, 3) - dt*fpp{jNC}( BVpp_left(:, 3), BVmm_left(:, 3) );
            ymm_o1_left(jNC) = BVmm_left(jNC, 1) + dt*fmm{jNC}( BVpp_left(:, 1), BVmm_left(:, 1) );
        end

        % Rotated MoC-ME for the left boundary:

        for   jNC = 1 : NumComp
            vpp_left(jNC) = (1/2)*( BVpp_left(jNC, 3) + ypp_o1_left(jNC) - ...
                                    dt*fpp{jNC}( ypp_o1_left, ymm_o1_left ) );
            vmm_left(jNC) = (1/2)*( BVmm_left(jNC, 1) + ymm_o1_left(jNC) + ...
                                    dt*fmm{jNC}( ypp_o1_left, ymm_o1_left ) );
        end

        % Rotated MoC-SE for the right boundary:

        for   jNC = 1 : NumComp
            ypp_o1_right(jNC) = BVpp_right(jNC, 1) + dt*fpp{jNC}( BVpp_right(:, 1), BVmm_right(:, 1) );
            ymm_o1_right(jNC) = BVmm_right(jNC, 3) - dt*fmm{jNC}( BVpp_right(:, 3), BVmm_right(:, 3) );
        end

        % Rotated MoC-ME for the right boundary:

        for   jNC = 1 : NumComp
            vpp_right(jNC) = (1/2)*( BVpp_right(jNC, 1) + ypp_o1_right(jNC) + ...
                                     dt*fpp{jNC}( ypp_o1_right, ymm_o1_right ) );
            vmm_right(jNC) = (1/2)*( BVmm_right(jNC, 3) + ymm_o1_right(jNC) - ...
                                     dt*fmm{jNC}( ypp_o1_right, ymm_o1_right ) ); 
        end
        
        % Reassign stored boundary values:
        
        for  jBV = 1:2
            BVpp_left(:, jBV) = BVpp_left(:, jBV+1);
            BVmm_left(:, jBV) = BVmm_left(:, jBV+1);
            
            BVpp_right(:, jBV) = BVpp_right(:, jBV+1);
            BVmm_right(:, jBV) = BVmm_right(:, jBV+1);
        end
            
    end
    
    
    % OG-b) Compute k{1,2}^(n-2) at m = 1,M and append them to the existing vectors:
    
    for   jNC = 1 : NumComp
        k1_pp_nmin1(jNC, 1) = fpp{jNC}( vpp_left,  vmm_left );
        k1_mm_nmin1(jNC, M) = fmm{jNC}( vpp_right, vmm_right );
    end
    
    % For k2, we need the following quantities first:
    
    ypp_o1_left = vpp_left + dt*k1_pp_nmin1(:, 1);
                  % Note: 
                  % We can recycle the name y{pp,mm}_o1_{left,right} 
                  % since we won't need it until next n.
    ymm_o2_left = ymm(:, 1);     
                  % Note:
                  % The notation "o2" on the lhs in this and the next line is 
                  % *intentionally incorrect*:  indeed, the solution ymm(:,1) 
                  % has already been computed with accuracy O(h^3).
                  % We use "o2" in the notations to maintain uniformity with 
                  % the notations used at the level  n = 2. 

    ypp_o2_right = ypp(:, M);
    ymm_o1_right = vmm_right + dt*k1_mm_nmin1(:, M);

    % Compute k2 at the virtual nodes and append to k2 in the bulk of the grid:

    for   jNC = 1 : NumComp
        k2_pp_nmin1(jNC, 1) = fpp{jNC}( ypp_o1_left,  ymm_o2_left );
        k2_mm_nmin1(jNC, M) = fmm{jNC}( ypp_o2_right, ymm_o1_right );
    end

   
    % Do calculations INSIDE the grid ........................................
	%
    % IG-a) Compute the MoC-ME solution which will be used to compute the new k2 (i.e., k2_nmin0)
    %       at the "previous" level (which is (n-1));  see the Note right after the start
    %       of the Main loop). 
        
    %   Compute the MoC-SE solution and append BC to it:
    
    ypp_o1(:, 2 : M) = ypp(:, 1:M-1) + dt*k1_pp_nmin0(:, 2 : M);
    ymm_o1(:, 1:M-1) = ymm(:, 2 : M) + dt*k1_mm_nmin0(:, 1:M-1);

    ypp_o1(:,1) = BCleft;    
    ymm_o1(:,M) = BCright;   
    
    %   Compute the MoC-ME solution proper and append BC to it:
    
    for   jNC = 1 : NumComp
        ypp_o2(jNC, 2 : M) = (1/2)*( ypp(jNC, 1:M-1) + ypp_o1(jNC, 2 : M) + ...
                                     dt*fpp{jNC}( ypp_o1(:, 2 : M), ymm_o1(:, 2 : M) ) );
        ymm_o2(jNC, 1:M-1) = (1/2)*( ymm(jNC, 2 : M) + ymm_o1(jNC, 1:M-1) + ...
                                     dt*fmm{jNC}( ypp_o1(:, 1:M-1), ymm_o1(:, 1:M-1) ) );
    end
                         
    ypp_o2(:, 1) = BCleft;    
    ymm_o2(:, M) = BCright;  
    
    
    % IG-b) Compute k2 at the *previous* level (n-1), in the bulk only:
    
    for   jNC = 1 : NumComp
        k2_pp_nmin0(jNC, 2 : M) = fpp{jNC}( ypp_o1(:, 2 : M), ymm_o2(:, 2 : M) );
        k2_mm_nmin0(jNC, 1:M-1) = fmm{jNC}( ypp_o2(:, 1:M-1), ymm_o1(:, 1:M-1) );
    end    
    
    
    % IG-c) Compute the MoC-pRK3 solution:
    
    ypp_n(:, 2 : M) = ypp(:, 1:M-1) + (dt/12)*( 13*k1_pp_nmin0(:, 2 : M) + 5*k2_pp_nmin0(:, 2 : M) -...
                                                 1*k1_pp_nmin1(:, 1:M-1) - 5*k2_pp_nmin1(:, 1:M-1) );                                 
    ymm_n(:, 1:M-1) = ymm(:, 2 : M) + (dt/12)*( 13*k1_mm_nmin0(:, 1:M-1) + 5*k2_mm_nmin0(:, 1:M-1) -...
                                                 1*k1_mm_nmin1(:, 2 : M) - 5*k2_mm_nmin1(:, 2 : M) );
                                  
    ypp_n(:, 1) = BCleft;   
    ymm_n(:, M) = BCright;  
        
    
    % IG-d) Update the solution and store new boundary values:
    
    ypp = ypp_n;
    ymm = ymm_n;
           
    BVpp_left(:, 3) = ypp(:, 1);    
    BVmm_left(:, 3) = ymm(:, 1);  

    BVpp_right(:, 3) = ypp(:, M);
    BVmm_right(:, 3) = ymm(:, M);
    
    
    % IG-e) Reassign k1 and k2 one level up:
    
    k1_pp_nmin1(:, 2 : M) = k1_pp_nmin0(:, 2 : M);
    k1_mm_nmin1(:, 1:M-1) = k1_mm_nmin0(:, 1:M-1);
    
    k2_pp_nmin1(:, 2 : M) = k2_pp_nmin0(:, 2 : M);
    k2_mm_nmin1(:, 1:M-1) = k2_mm_nmin0(:, 1:M-1);

    
    % IG-f) Compute k1 at the *current* level:
    
    for   jNC = 1 : NumComp
        k1_pp_nmin0(jNC, 2 : M) = fpp{jNC}( ypp(:, 1:M-1), ymm(:, 1:M-1) );
        k1_mm_nmin0(jNC, 1:M-1) = fmm{jNC}( ypp(:, 2 : M), ymm(:, 2 : M) );
    end
    
        
    % This block plots the Fourier spectrum of the solution  --------------------------------
    % ---------------------------------------------------------------------------------------
    
    if rem(n, nplot) == 0
        total_spectrum_sq = zeros(size(kk));
        for   jNC = 1 : NumComp
            total_spectrum_sq = total_spectrum_sq + ...
                  abs(fft(ypp(jNC, :).*spatial_window)).^2 + abs(fft(ymm(jNC, :).*spatial_window)).^2;
        end
        %
         figure(31)
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
    
end

