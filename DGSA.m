

%% DGSA: Modified from main_DGSA_Reservoir_Sensitivity.m

% Part1: Compute main/conditional effects of Prior models

clear all; close all; fclose('all'); rng('default');
%% 0. Add directories for path
addpath(genpath('/Users/ianpg/Documents/ProjectsLocal/DGSA-master'))

%% 1. Specify inputs for DGSA. 

%load hausdorff matrix

exptdir = '/Users/ianpg/Documents/ProjectsLocal/DelawareSGD/work/homogenous/MC_expt_2018-12-20-20-39';
f = 'hausdorff.mat';
fname = fullfile(exptdir,f);
load(fname)

%This is added for the to fix the messed up hdorf matrix
hausdorff_mat(496:end,:)= [];
hausdorff_mat(:,496:end)= [];
%ParametersValues(:,4:7) = log10(ParametersValues(:,4:7));
%% 2. Load & process input data
D = hausdorff_mat;
N = length(D);
ParametersValues = ParametersValues;
ParametersNames = fieldnames(InputParams)';


DGSA={};
DGSA.D = D;
DGSA.N = N;
DGSA.ParametersValues = ParametersValues;
DGSA.ParametersNames = ParametersNames;

%% 4. Compute & display main effects

% 4.1 Inputs for clustering and display options.
DGSA.Nbcluster=2; % # of clusters
DGSA.MainEffects.Display.ParetoPlotbyCluster=1; % if true, main effects over cluster will be displayed with Pareto plot.
DGSA.MainEffects.Display.StandardizedSensitivity='CI'; 

% if 'CI', the confidence interval will be overlapped on the Pareto plot
% if 'Pareto', Pareto plots will be displayed only
% if 'None' No plot will be generated for standardd main effects
% (default)

% 4.2 Compute main effects from DGSA.
SigLevel_CI=[.95,.9,1]; % This is needed only when you want to display confidence intervals
                        % If you do not need to display confidence intervals,set DGSA.MainEffects.Display.StandardizedSensitivity='Preto' or 'None'     

% 4.3 Perform clustering
DGSA.Clustering=kmedoids(DGSA.D,DGSA.Nbcluster,10); % In this example, K medoid clustering is applied.

% 4.4 Compute Main Effects
DGSA=ComputeMainEffects(DGSA,SigLevel_CI); % If you use Pareto plot or do not want to display main effects, remove SigLevel_CI, otherwise it shows an error.

%% Display cdfs

cdf_MainFactor(DGSA.ParametersValues, DGSA.Clustering, DGSA.ParametersNames); %last arg can be cell of params e.g. {'al','vka'}
%% 5. Compute & Display conditional effects

% 5.1 Specify additional variables to estimate conditional effects.

DGSA.ConditionalEffects.NbBins=2*ones(1,length(DGSA.ParametersNames));

% 5.2 Compute conditional effects

rng('default'); 
DGSA=ComputeConditionalEffects(DGSA);

% 5.3 Display conditional effcts

% Speicify the method to display standardized conditional effects.
DGSA.ConditionalEffects.Display.SensitivityByClusterAndBins=1; % if true, display pareto plots to visualize sensitivities by bins/clusters. 
DGSA.ConditionalEffects.Display.StandardizedSensitivity='Hplot'; % If omitted Pareto plot will be used. However, this is not recommended when there are lots of parameters
%%
% Visualize conditional effects
DisplayConditionalEffects(DGSA,DGSA.ConditionalEffects.Display.StandardizedSensitivity)

% 5.4 Display class condtional CDFs
cdf_ConditionalEffects('hk','vka',DGSA,1)

%% Save all variables for futher application
%save('VariablesSaved/DGSA_Completed.mat');
%save('VariablesSaved/DGSA_Spatial_Completed.mat');