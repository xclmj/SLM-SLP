disp('Initilization_Condition')

clear
linsolve = @mldivide;
try
    require deckformat incomp
catch
    mrstModule add deckformat incomp
end

clear

mrstModule add incomp

%% Define geometry and rock parameters
%% Initial pressure and water saturation states
p0 = 300*barsa();
s0=[0.2, 0.8];
gravity off

% weighting coefficient for two facies
gammaF=0.8;  % weighting coefficient for facies

%% Create a fluid model: viscosity, density, and Corey-type relative
% permeability functions
fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                       'rho', [1014, 859]*kilogram/meter^3, ...
                       'n'  , [   2,   2]                 , ...
                       'sr' , [ 0.2, 0.2]                 , ...
                       'kwm', [   1,   1]);
%% Create a fluid model: viscosity, density, and Corey-type relative
% permeability functions
% three uncertain parameters
% fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
%                        'rho', [1014, 859]*kilogram/meter^3, ...
%                        'n'  , [   4,   4]                 , ...
%                        'sr' , [ 0.2, 0.2]                 , ...
%                        'kwm', [   1,   1]);
% n_min=2; n_max=4; n_mean=3; n_std=0.5;
% S_min=0.1; S_max=0.3; S_mean=0.2; S_std=0.1;
% kr_min=0.8; kr_max=1; kr_mean=0.9; kr_std=0.1;

%% define time schedule
linsolve_p = @mldivide;  % Pressure
linsolve_t = @mldivide;  % Transport (implicit)

%history matching 
Tend           = 5*year();          % end time
dT             = 0.1*year();       % simulation timestep size 0.01*year() (30*day())
%dT             = 10*day();         % simulation timestep size 0.01*year() (30*day())
Tsim=15*year();

%predection
nt = ceil(Tend/dT);                 % number of timestep

