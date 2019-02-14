

%% DGSA: Modified from main_DGSA_Reservoir_Sensitivity.m

% Part1: Compute main/conditional effects of Prior models

clear all; close all; fclose('all'); rng('default');
%% 0. Add directories for path
if ispc()
    addpath(genpath('E:\MATLAB\Packages\DGSA-master_altered'));
else
    addpath(genpath('/Users/ianpg/Documents/ProjectsLocal/DGSA-master'))
end

%% 1. Specify inputs for DGSA. 

%load hausdorff matrix
exptdir = '/Users/ianpg/Dropbox/TempShared/MC_expt_2019-01-08-02-50';
exptdir = 'E:\Projects\DelawareSGD\work\passive_active\MC_expt_2019-02-12-11-55'
finput = 'hausdorff.mat';
printplotsyn = true;
%fhausdorff = 'hausdorff.csv';
load(fullfile(exptdir,finput));
try
    load(fullfile(exptdir,'culled.mat'));
catch
end
%% 
sz = size(hausdorff_mat);
if sz(1)~=sz(2)
    hausdorff_mat = squareform(hausdorff_mat);
end
%% 2. Load & process input data
D = hausdorff_mat;
N = length(D);
ParametersValues = ParametersValues;
ParametersNames = fieldnames(InputParams)';
ParametersNames = names2latex(ParametersNames);
%%
ParametersValues2 = ParametersValues;
logvalues = [4,5,6,7,8,9,10,12];
for i=1:12
    if ismember(i,logvalues)
        ParametersValues2(:,i) = log10(ParametersValues(:,i));
    end
end


DGSA={};
DGSA.D = D;
DGSA.N = N;
DGSA.ParametersValues = ParametersValues;
DGSA.ParametersNames = ParametersNames;

%%

zeroResp = D(:,1)==0;
logvalues = [4,5,6,7,8,9,10,12];
log10_str = @(x) log10(str2num(x));

fig = figure;
for i=1:12
    ax = subplot(4,3,i);
    if ismember(i,logvalues)
        val = log10(ParametersValues(:,i));
    else
        val = ParametersValues(:,i);
    end
    h = histogram(val(zeroResp),'DisplayName','zero-value',...
        'Normalization','Probability');
    binedge = h.BinEdges;
    hold on
    h2 = histogram(val(~zeroResp),'DisplayName','value',...
                'Normalization','Probability');
    %h2.BinEdges = binedge;
    title(ParametersNames{i});
    if i==1
        legend;
    end
end
%% 4. Cluster and Compute Main Effects

% 4.1 Inputs for clustering and display options.
DGSA.Nbcluster=3; % # of clusters
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
medoids = DGSA.Clustering.medoids;
save(fullfile(exptdir,'medoids.mat'),'medoids')
med_mats = culled_conc_mat(medoids,:,:,:);
med_mats(med_mats>=1e30)=NaN;
%%
%Plot MDS results and representative plots
plotreps=true;
rowslice= 2;
printplotsyn=true;

data2pos = @(data,datalims,axlims) (data-datalims(1))/diff(datalims)*diff(axlims) + axlims(1)

figure;
[ClusterColor,pts]=DisplayMDSplot(DGSA.D, DGSA.Clustering);
title('MDS Plot: Clustering by response')
DGSA.Clustering.ClusterColor = ClusterColor;
ax= gca;
h=.8;
ax.OuterPosition = [0,0,1,h];
ax_pos = ax.Position;
pts.MarkerFaceColor = 'flat';
pts.MarkerEdgeColor = 'flat';


if plotreps
    hold on
    w_h = .15;
    pos = [.15,.5,w_h,w_h;
           .25,.7,w_h,w_h;
           .45,.75,w_h,w_h;
           .65,.7,w_h,w_h;
           .8,.6,w_h,w_h];
    pos = [.05,h+.02,w_h,w_h;
           .25,h+.02,w_h,w_h;
           .45,h+.02,w_h,w_h;
           .65,h+.02,w_h,w_h;
           .83,h+.02,w_h,w_h];   
    pos = linspace(0.05,.95,DGSA.Nbcluster+1)';
    w_h = min([(pos(2)-pos(1))*.75,1-(h+.02)]);
    pos = horzcat(pos(1:end-1),repmat(h+.02,DGSA.Nbcluster,1),repmat(w_h,DGSA.Nbcluster,1),repmat(w_h,DGSA.Nbcluster,1));
       
    [~,sortind] = sort(pts.XData(medoids));
    for i=1:length(medoids)
        ax2 = axes('Position',pos(i,:),'Color','none');
        imagesc(squeeze(med_mats(sortind(i),:,rowslice,:)));
        ax2.XTickLabel='';
        ax2.YTickLabel='';
        x=[data2pos(pts.XData(medoids(sortind(i))),ax.XLim,[ax.Position(1),ax.Position(1)+ax.Position(3)]),...
            ax2.Position(1)+.5*ax2.Position(4)];
        y=[data2pos(pts.YData(medoids(sortind(i))),ax.YLim,[ax.Position(2),ax.Position(2)+ax.Position(4)]),...
            ax2.Position(2)];
        annotation('arrow',x,y,'LineStyle','--','HeadStyle','vback3')
        colormap(ax2,jet(50))
    end
end
if printplotsyn 
    if plotreps
        namext='wreps';
    else
        namext='';
    end
    namext = strcat(namext,num2str(DGSA.Nbcluster));
    print(fullfile(exptdir,strcat('clustering',namext)),'-dtiff','-r200')
end

%% Plot 4 corners of MDS plot as well
rowslice= 10;

[~,minx] = min(pts.XData);
[~,maxx] = max(pts.XData);
[~,miny] = min(pts.YData);
[~,maxy] = max(pts.YData);
vals = horzcat(minx,maxx,miny,maxy);
figure;
for i=1:length(vals)
    disp(i)
    mat = squeeze(culled_conc_mat(vals(i),:,rowslice,:));
    mat(mat>=1e30)=NaN;
    subplot(2,2,i);
    imagesc(mat)
    if i==1
        title('Min X Value');
    elseif i==2
        title('Max X Value');
    elseif i==3
        title('Min Y Value');
    elseif i==4
        title('Max Y Value');
    end
        
end


%% Sort and plot figures to see how they match up

numplots=19;
startind=1;
every=10;
rowslice= 10;
[~,sortind] = sort(pts.XData);

figure;
ind=1;
for i=startind:every:startind+numplots*every
    mat = squeeze(culled_conc_mat(sortind(i),:,rowslice,:));
    mat(mat>=1e30)=NaN;
    ax=subplot(ceil(sqrt(numplots)),ceil(numplots/floor(sqrt(numplots))),ind);
    imagesc(mat);
    ax.XTickLabel='';
    ax.YTickLabel='';
    colormap(ax,jet(50))
    ind= ind+1;
    title(sprintf('X-value %d',i))
end


%% 4.4 Compute Main Effects
DGSA=ComputeMainEffects(DGSA,SigLevel_CI); % If you use Pareto plot or do not want to display main effects, remove SigLevel_CI, otherwise it shows an error.
if printplotsyn 
    print(fullfile(exptdir,'maineffects'),'-dtiff','-r150')
end
%% 4.5 Display cdfs

cdf_MainFactor(DGSA.ParametersValues, DGSA.Clustering, DGSA.ParametersNames,{'al','hk1','hk2'}); %last arg can be cell of params e.g. {'al','vka'}
if printplotsyn 
    %print(strcat(exptdir,'/cdfs'),'-dtiff','-r150')
end

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
%%
% 5.4 Display class condtional CDFs
cdf_ConditionalEffects('hk','vka',DGSA,1)

%% Save all variables for futher application
%save('VariablesSaved/DGSA_Completed.mat');
%save('VariablesSaved/DGSA_Spatial_Completed.mat');



%% Custom functions

%Change names with underscores to be read by latex as a subscript
function out_cellarray = names2latex(name_cellarray)
   out_cellarray = name_cellarray;
    for i=1:length(name_cellarray)
        sp = split(name_cellarray{i},'_');
        if length(sp)>1
            out_cellarray{i} = strcat(sp{1},'_{',sp{2:end},'}');
        else
            out_cellarray{i} = sp{1};
        end
        
        %out_cellarray{i} = strcat('$$',out_cellarray{i},'$$');
    end
end 