% Boar & Midrigan (2022) - Efficient Redistribution
% Authors: Originally written by Robert Kirkby, later modified by
% Alessandro Di Nola

clear,clc,close all
mytoolkit = 'C:\Users\dinolaa\Documents\GitHub\VFIToolkit-matlab';
addpath(genpath(mytoolkit));

%% Grid sizes
n_d = 51;  % Grid points for labor supply
n_a = 301; % Grid points for assets
n_z = 12;  % Grid points for labor productivity
           % (11 to discretize AR(1), plus one for the super-star state)

figure_c=0; % counter for figures

%% Toolkit options
vfoptions.gridinterplayer = 1;
vfoptions.ngridinterp     = 20;
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;

%% Parameters

% Preferences
Par.beta  = 0.975; % discount factor
Par.theta = 1;     % curvature of utility of consumption (CRRA parameter)
Par.gamma = 2;     % curvature of utility of leisure (inverse Frisch elasticity)

% Production
Par.alpha = 0.333; % capital share of income
Par.delta = 0.06;  % depreciation rate

% Taxes
% Consumption tax
Par.tau_s = 0.065; 
% Income tax
Par.tau = 0.263;
Par.xi  = 0.049;
% Wealth tax
Par.tau_a = 0;
Par.xi_a  = 0;
% Lump-sum transfer 
Par.iota_target = 0.167; % relative to per-capita GDP
Par.iota        = 0.1; % guess, will be calibrated so we hit iota_target

% Government debt-to-GDP
Par.Bbar=1;

% Earnings productivity process
Par.rho_z   = 0.982;  % autocorrelation of z
Par.sigma_e = 0.2;    % std dev of innovations to z
Par.p       = 2.2e-6; % probability to enter the super-star state
Par.q       = 0.990;  % probability to stay in the super-star state
Par.zbar    = 504.3;  % ability super-star state relative to mean

% Following parameters/prices are deteremined in general eqm, these are just
% initial guesses (I had worse guesses on the first run, these are updated
% so the initial general eqm is quicker to solve by using less bad guesses)
Par.r    = 0.05;   
Par.iota = 0.3; 
Par.B    = 1.8;   
Par.G    = 0.05;

%% Generate grids and transitions
d_grid=linspace(0,1,n_d)'; % labor supply grid

Par.maxa = 40;  % max assets
Par.mina = 0.0; % min assets
a_grid = Par.mina+(Par.maxa-Par.mina)*(linspace(0,1,n_a)'.^3); % assets grid

% BM2022 use Rouwenhorst
[z_grid_pre,pi_z_pre]=discretizeAR1_Rouwenhorst(0,Par.rho_z,Par.sigma_e,n_z-1); 
% Make it so that mean of z is 1
z_grid_pre=exp(z_grid_pre);
[meanz,~,~,~]=MarkovChainMoments(z_grid_pre,pi_z_pre);
z_grid_pre=z_grid_pre./meanz;
[meanz,~,~,statdistz]=MarkovChainMoments(z_grid_pre,pi_z_pre);
% Add in the super-star state
z_grid=[z_grid_pre; Par.zbar*meanz]; % meanz is very close to one in any case
pi_z=[(1-Par.p)*pi_z_pre, Par.p*ones(n_z-1,1); (1-Par.q)*statdistz', Par.q];
% Page 82, "We assume that agents transit from the normal to the super-star
% state with a constant probability p and remain there with probability q.
% When agents return to the normal state, they draw a new ability from the
% ergodic distribution associated with the AR(1) process."

% BM write that "0.02% of households are in the super-star state at any
% point in time". Let's check if this is the case

[~,~,~,statdistz]=MarkovChainMoments(z_grid_pre,pi_z_pre);

%% Return fn and discount factor
DiscountFactorParamNames={'beta'};

ReturnFn=@(h,aprime,a,z,r,tau_s,tau,xi,tau_a,xi_a,iota,theta,gamma,delta,alpha)...
    BoarMidrigan2022_ReturnFn(h,aprime,a,z,r,tau_s,tau,xi,tau_a,xi_a,iota,theta,gamma,delta,alpha);

%% Test the value fn and policy fn
tic;
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Par,DiscountFactorParamNames,[],vfoptions);
vftime=toc

%% Test for Agent Dist
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,simoptions);

%% Setup model moments
FnsToEvaluate.K=@(h,aprime,a,z) a;
FnsToEvaluate.L=@(h,aprime,a,z) h*z;
FnsToEvaluate.TaxRevenue=@(h,aprime,a,z,r,tau,xi,iota,delta,alpha,tau_a,xi_a) BM2022_IncomeTaxRevenue(h,aprime,a,z,r,tau,xi,iota,delta,alpha) + BM2022_WealthTaxRevenue(h,aprime,a,z,tau_a,xi_a);

% Test the model moments
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist,Policy,FnsToEvaluate,Par,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

%% Solve initial stationary general eqm
GEPriceParamNames_pre={'r','iota','B','G'};

GeneralEqmEqns_pre.capitalmarket=@(r,K,L,alpha,delta) r-(alpha*(K^(alpha-1))*(L^(1-alpha))-delta); % interest rate equals marginal product of capital (net of depreciation)
GeneralEqmEqns_pre.iotacalib=@(iota,iota_target,K,L,alpha) iota_target - iota/((K^(alpha))*(L^(1-alpha))); % get iota to GDP-per-capita ratio correct
% Note: because iota is same for everyone, just use it directly here rather than needing to calculate it as an aggregate var.
GeneralEqmEqns_pre.govbudget=@(r,B,G,TaxRevenue) (1+r)*B+G-(B+TaxRevenue); % Balance the goverment budget
% Include the calibration target as a general eqm constraint
GeneralEqmEqns_pre.govdebtcalib=@(Bbar,B,K,L,alpha) Bbar - B/((K^(alpha))*(L^(1-alpha))); % Gov Debt-to-GDP ratio of Bbar

heteroagentoptions.verbose=1;
[p_eqm_init,GECondns_init]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns_pre, Par, DiscountFactorParamNames, [], [], [], GEPriceParamNames_pre,heteroagentoptions, simoptions, vfoptions);
% Update Par based on general eqm
Par.r=p_eqm_init.r;
Par.iota=p_eqm_init.iota;
Par.B=p_eqm_init.B;
Par.G=p_eqm_init.G;

[V_init,Policy_init]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Par,DiscountFactorParamNames,[],vfoptions);
StationaryDist_init=StationaryDist_Case1(Policy_init,n_d,n_a,n_z,pi_z,simoptions);

% Double-check that no-one is hitting the max of the a_grid
asset_cdf=cumsum(sum(StationaryDist_init,2));
figure_c=figure_c+1;
figure(figure_c);
plot(a_grid,asset_cdf)
title('cdf of asset in initial stationary general eqm (check not hitting the top of grid)')

% Add some further FnsToEvaluate so can plot the kinds of outputs shown in Figure 3
FnsToEvaluate2=FnsToEvaluate;
FnsToEvaluate2.Consumption=@(h,aprime,a,z,r,tau_s,tau,xi,tau_a,xi_a,iota,delta,alpha) BM2022_Consumption(h,aprime,a,z,r,tau_s,tau,xi,tau_a,xi_a,iota,delta,alpha);

AggVars_init=EvalFnOnAgentDist_AggVars_Case1(StationaryDist_init,Policy_init,FnsToEvaluate2,Par,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
AllStats_init=EvalFnOnAgentDist_AllStats_Case1(StationaryDist_init,Policy_init,FnsToEvaluate2,Par,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

wage_init=(1-Par.alpha)*((p_eqm_init.r+Par.delta)/Par.alpha)^(Par.alpha/(Par.alpha-1));
output_init=(AggVars_init.K.Mean^(Par.alpha))*(AggVars_init.L.Mean^(1-Par.alpha));

C_t=AggVars_init.Consumption.Mean;
C_tplus1=AggVars_init.Consumption.Mean^(-Par.theta); % C_t and C_t+1 are same thing in stationary general eqm
laborwedge_init=wage_init*(C_t^(-Par.theta))/(AggVars.L.Mean^Par.gamma); % Rearrange eqn at bottom of page 80 for varthetabar
savingswedge_init=Par.beta*(1+Par.r)*(C_t^(-Par.theta))/(C_tplus1^(-Par.theta)); % Rearrange eqn (4) on pg 81 to get zetabar_t, the aggregate capital wedge 

%% From now on, just keep G and B unchanged
GEPriceParamNames={'r','iota'};
GeneralEqmEqns.capitalmarket=GeneralEqmEqns_pre.capitalmarket;
GeneralEqmEqns.govbudget=GeneralEqmEqns_pre.govbudget;

%% Boar and Midrigan do optimal, which involves comparing lots of possible tax reforms.
% According to their computational appendix, they first did a rough grid on
% tax rates, then started an optimization routine from the best point on
% this grid [to get the welfare maximizing tax rates].
% Note, once we solve one we can use solution as initial guess for the next, so
% runtime would be much less than simply repeating this exercise.

%% Here we just do one tax reform, to show how it is done.

% Tax reform (this is the one they found to be optimal for 'utilitarian welfare'; I had to eyeball the tau and tau_a out of Figure 4 as the exact number does not appear to be in paper, the xi and xi_a are explicit in Figure 4)
Par.tau=0.56; % income tax
Par.xi=0.065;
Par.tau_a=-0.002; % wealth tax
Par.xi_a=0.0017;

%% Final stationary general eqm
tic;
[p_eqm_final,GECondns_final]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Par, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
GEtime2=toc
% Update Par based on general eqm
Par.r=p_eqm_final.r;
Par.iota=p_eqm_final.iota;

[V_final,Policy_final]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Par,DiscountFactorParamNames,[],vfoptions);
StationaryDist_final=StationaryDist_Case1(Policy_final,n_d,n_a,n_z,pi_z,simoptions);


save BM2022pre.mat

%% Solve the transition path
T=100; % BM2022 Fig 3 has this as x-axis, so seems appropriate

% Initial guess for price path
PricePath0.r=[linspace(p_eqm_init.r,p_eqm_final.r,ceil(T/3)),p_eqm_final.r*ones(1,T-ceil(T/3))];
% PricePath0.iota=[linspace(p_eqm_init.iota,p_eqm_final.iota,ceil(T/3)),p_eqm_final.iota*ones(1,T-ceil(T/3))];
PricePath0.iota=[linspace(p_eqm_init.iota,0.8*p_eqm_final.iota,ceil(T/3)),0.8*p_eqm_final.iota*ones(1,T-ceil(T/3)-1),p_eqm_final.iota]; % deliberately start without enough iota

% Parameter path is trivial, as they are preannounced one-off reforms
ParamPath.tau=Par.tau*ones(1,T);
ParamPath.xi=Par.xi*ones(1,T);
ParamPath.tau_a=Par.tau_a*ones(1,T);
ParamPath.xi_a=Par.xi_a*ones(1,T);
% B is constant, so can skip putting path on it as long as value in Par is correct one.

% Same as before, except that the government budget now has to be explicit that it is last period gov debt.
GeneralEqmEqns_TransPath.capitalmarket=GeneralEqmEqns.capitalmarket;
% GeneralEqmEqns_TransPath.govbudget=@(r,B,B_tminus1,G,TaxRevenue) (1+r)*B_tminus1+G-(B+TaxRevenue);
% Note: As discussed at start of this script BM2022 are presumably solving for a path where B is constant over time, but the commented out line above still writes out B_tminus1 so that it would still work for other setups.
GeneralEqmEqns_TransPath.govbudget=@(r,B,G,TaxRevenue) (1+r)*B+G-(B+TaxRevenue); % Take advantage of constant B

transpathoptions.GEnewprice=3;
% Need to explain to transpathoptions how to use the GeneralEqmEqns to update the general eqm transition prices (in PricePath).
transpathoptions.GEnewprice3.howtoupdate=... % a row is: GEcondn, price, add, factor
    {'capitalmarket','r',0,0.3;...  % captialMarket GE condition will be positive if r is too big, so subtract
    'govbudget','iota',0,0.05;... % govbudget GE condition will be positive if iota is too big, so subtract [iota is subtracted from the tax, so bigger iota means smaller TaxRevenue; note, units of iota are essentially the same as units of TaxRevenue, so don't need to think too hard about the size of the factor]
    };
% Note: the update is essentially new_price=price+factor*add*GEcondn_value-factor*(1-add)*GEcondn_value
% Notice that this adds factor*GEcondn_value when add=1 and subtracts it what add=0
% A small 'factor' will make the convergence to solution take longer, but too large a value will make it 
% unstable (fail to converge). Technically this is the damping factor in a shooting algorithm.


% For the transition path, turn on divide and conquer
vfoptions.divideandconquer=1;
vfoptions.level1n=11; % might be slightly faster or slower with a higher/lower value (do a tic-toc on ValueFnOnTransPath_Case1() to find the fastest if you want to find out what to set this to)

transpathoptions.verbose=1;
tic;
PricePath=TransitionPath_Case1(PricePath0, ParamPath, T, V_final, StationaryDist_init, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Par, DiscountFactorParamNames, transpathoptions, vfoptions, simoptions);
tpathtime=toc

[VPath,PolicyPath]=ValueFnOnTransPath_Case1(PricePath, ParamPath, T, V_final, Policy_final, Par, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);

AgentDistPath=AgentDistOnTransPath_Case1(StationaryDist_init, PolicyPath,n_d,n_a,n_z,pi_z,T,simoptions);

AggVarsPath=EvalFnOnTransPath_AggVars_Case1(FnsToEvaluate2,AgentDistPath,PolicyPath,PricePath,ParamPath, Par, T, n_d, n_a, n_z, pi_z, d_grid, a_grid,z_grid,simoptions);
% AllStatsPath=EvalFnOnTransPath_AllStats_Case1()

wagepath=((1-Par.alpha)*((PricePath.r+Par.delta)/Par.alpha).^(Par.alpha/(Par.alpha-1)))';
outputpath=(AggVarsPath.K.Mean.^(Par.alpha)).*(AggVarsPath.L.Mean.^(1-Par.alpha));

C_t=AggVarsPath.Consumption.Mean;
C_tplus1=[AggVarsPath.Consumption.Mean(2:end),AggVarsPath.Consumption.Mean(end)].^(-Par.theta); % C_t+1 becomes constant in period T, so I just duplicate this
laborwedgepath=wagepath.*(C_t.^(-Par.theta))./(AggVarsPath.L.Mean.^Par.gamma); % Rearrange eqn at bottom of page 80 for varthetabar
savingswedgepath=Par.beta*(1+Par.r)*(C_t.^(-Par.theta))./(C_tplus1.^(-Par.theta)); % Rearrange eqn (4) on pg 81 to get zetabar_t, the aggregate capital wedge 

save BM2022.mat

figure_c=figure_c+1;
figure(figure_c);
subplot(2,4,1); plot(0:1:T, [laborwedge_init,laborwedgepath]) 
title('Labor wedge')
subplot(2,4,2); plot(0:1:T, [savingswedge_init,savingswedgepath]) 
title('Savings wedge')
subplot(2,4,3); plot(0:1:T, [output_init,outputpath]) 
title('Output')
subplot(2,4,4); plot(0:1:T, [AggVars_init.K.Mean,AggVarsPath.K.Mean]) 
title('Capital')
subplot(2,4,5); plot(0:1:T, [AggVars_init.L.Mean,AggVarsPath.L.Mean]) 
title('Labor')
subplot(2,4,6); plot(0:1:T, [p_eqm_init.r,PricePath.r']) 
title('Interest rate')
subplot(2,4,7); plot(0:1:T, [wage_init,wagepath])
title('Wage')
% subplot(2,4,8); plot(1:1:T,[AllStats_init.K.Gini,AllStatsPath.K.Gini])
% title('Gini wealth')


%% Welfare calculations
% Paper doesn't seem to specify which period welfare is defined in, but presumably 
% the period in which the reform is revealed. In which case we already have V, 
% it is the first time period in VPath, and we already have the agent
% distribution, it is StationaryDist_init (which is also the first time
% period in StationaryDistPath)

% Social welfare preferences
Par.Delta=1; 
  % 0: average welfare
  % 1: utilitarian
  % Inf: Rawlsian

% Eqn on page 83 defining omega, with infinite sum algebra (as omega is
% constant and beta is less than 1) gives us
omegamod=VPath(:,:,1)*(1-Par.beta);
if Par.theta==1
    omega=exp(omegamod); % as (omega^(1-theta))/(1-theta) becomes ln(omega)
else
    omega=(omegamod*(1-Par.theta))^(1/(1-Par.theta));
end

% Evaluate the social welfare function (from pg 83)
if ~isfinite(Par.Delta) % Rawlsian
    % Not clear how they actually did Rawlsian, taking the min seems extreme but 
    % they do not specify a value of Delta, only say limit as Delta goes to infinity
    temp=omega(StationaryDist_init>0);
    SocialWelfare=min(temp);
else
    temp=(omega.^(1-Par.Delta)).*StationaryDist_init;
    temp(isnan(temp))=0; % in case there are any -Inf*0, which would give nan
    SocialWelfare=sum(sum(temp)).^(1/(1-Par.Delta));
end

% Repeat the welfare calculation, but now for the initial welfare level in
% the initial stationary general eqm
% Eqn on page 83 defining omega, with infinite sum algebra (as omega is
% constant and beta is less than 1) gives us
omegamod0=V_init*(1-Par.beta);
if Par.theta==1
    omega0=exp(omegamod0); % as (omega^(1-theta))/(1-theta) becomes ln(omega)
else
    omega0=(omegamod0*(1-Par.theta))^(1/(1-Par.theta));
end
% Evaluate the social welfare function (from pg 83)
if ~isfinite(Par.Delta) % Rawlsian
    % Not clear how they actually did Rawlsian, taking the min seems extreme but 
    % they do not specify a value of Delta, only say limit as Delta goes to infinity
    temp=omega0(StationaryDist_init>0);
    SocialWelfare0=min(temp);
else
    temp=(omega0.^(1-Par.Delta)).*StationaryDist_init;
    temp(isnan(temp))=0; % in case there are any -Inf*0, which would give nan
    SocialWelfare0=sum(sum(temp)).^(1/(1-Par.Delta));
end

%% Plot some other things from paper. First, the tax functions like in Figure 4

income_grid=linspace(0,10,1000);
marginalincometaxrate=1-(1-Par.tau)*(income_grid.^(-Par.xi)); % derivative of tax fn on pg 80 w.r.t. income
averageincometaxrate=(income_grid-(1-Par.tau)*((income_grid.^(1-Par.xi))/(1-Par.xi))-Par.iota)./income_grid; % taxes paid divided by income
wealth_grid=linspace(0,10,1000);
marginalwealthtaxrate=1-(1-Par.tau_a)*(wealth_grid.^(-Par.xi_a)); % derivative of tax fn on pg 80 w.r.t. income
averagewealthtaxrate=(wealth_grid-(1-Par.tau_a)*((wealth_grid.^(1-Par.xi_a))/(1-Par.xi_a)))./wealth_grid; % taxes paid divided by income

figure_c=figure_c+1;
figure(figure_c);
yyaxis left
subplot(1,2,1); plot(income_grid,averageincometaxrate,income_grid,marginalincometaxrate)
hold on
yyaxis right
% subplot(1,2,1); plot(income_grid,income_cdf,'.')
ylim([0,1])
hold off
xlabel('income')
ylabel('%')
legend('average income tax', 'marginal income tax','initial cdf of income')
yyaxis left
subplot(1,2,2); plot(wealth_grid,averagewealthtaxrate,wealth_grid,marginalwealthtaxrate)
hold on
yyaxis right
subplot(1,2,2); plot(a_grid,asset_cdf,'.')
ylim([0,1])
hold off
xlabel('wealth')
ylabel('%')
legend('average wealth tax', 'marginal wealth tax','initial cdf of wealth')
% I modified this figure so it also shows the pdf of agent distribution
% Otherwise you cannot tell where anyone actually is

%% Calculate welfare for the terciles of welfare distribution (like what gets plotted in Fig 2)


%% Calculate tax revenue for the terciles of the ??? paper calls them 'richest third', but is it wealth or income, or different for top and bottom panels?? (like what gets plotted in Fig 5)



