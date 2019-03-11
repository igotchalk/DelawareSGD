

%% DGSA: Modified from main_DGSA_Reservoir_Sensitivity.m

% Part1: Compute main/conditional effects of Prior models

clear all; close all; fclose('all'); rng('default');
%% 0. Add directories for path
if ispc()
    addpath(genpath('E:\MATLAB\Packages\DGSA-master_altered'));
    addpath(genpath('E:\Projects\DelawareSGD'));
else
    addpath(genpath('/Users/ianpg/Documents/ProjectsLocal/DGSA-master'))
    addpath(genpath('/Users/ianpg/Documents/ProjectsLocal/DelawareSGD'))
end

%% 1. Specify inputs for DGSA. 

%load hausdorff matrix
%exptdir = '/Users/ianpg/Dropbox/TempShared/MC_expt_2019-01-08-02-50';
%exptdir = 'E:\Projects\DelawareSGD\work\passive_active\MC_expt_2019-02-12-11-55';
%exptdir ='E:/Projects/DelawareSGD/work/passive_active/MC_expt_2019-02-14-14-26'; %200 it, heterogenous, changed head_inland_sum
%exptdir ='E:/Projects/DelawareSGD/work/passive_active/MC_expt_2019-02-18-17-52/'; %500it,heterogenous,riv_cond updated (was way too high) 
exptdir = 'E:\Projects\DelawareSGD\work\mps\MC_expt_2019-03-02-18-49\'; %500it, BC prior should be identical to MC_expt_2019-02-18-17-52
addpath(exptdir);
totims = {'2340.0','4860.0','7200.0'};
tim = totims(end);
finput = char(strcat('hausdorff',char(tim),'.mat'));
finput_hk = 'hk_mat.mat';

% exptdir = 'E:\Projects\DelawareSGD\work\mps\MC_expt_2019-03-05-18-08'; %250it, heterogenous, Lt=40yrs
% pers = (180:180:40*360);
% totims = pers(1:10:end);
% tim = totims(end);
% finput = char(strcat('hausdorff',num2str(tim),'.mat'));

% finput = char(strcat('hausdorff',num2str(totims(1)),'.mat'));
printplotsyn = true;
%fhausdorff = 'hausdorff.csv';
load(fullfile(exptdir,finput));
try
    load(fullfile(exptdir,'culled.mat'));
catch
end

sz = size(hausdorff_mat);
if sz(1)~=sz(2)
    hausdorff_mat = squareform(hausdorff_mat);
end

%% Load spatial matrices and calculate rank

load(finput_hk);
hk_flat = reshape(hk_mat,size(hk_mat,1),[])';

%run the KPCASOM_rank script
%we used Gaussian radial basis function kernel
%For Gaussian RBF, kernel_para is for setting up the sigma value in RBF function 
% sigma= kernel_para * mean distance between models 
kernel_para=5;
[Y_eig,ndim,rank]=KPCASOM_rank(hk_flat,'gaussian',kernel_para);

%% 2. Load & process input data
D = hausdorff_mat;
N = length(D);
ParametersNames = fieldnames(InputParams)';
ParametersNames = names2latex(ParametersNames);

ParametersValues = horzcat(ParametersValues,rank');
ParametersNames{end+1} = 'hk_spatial';



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


%% 4. Cluster and Compute Main Effects

% 4.1 Inputs for clustering and display options.
DGSA.Nbcluster=4; % # of clusters
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
rowslice= 15;
printplotsyn=true;

data2pos = @(data,datalims,axlims) (data-datalims(1))/diff(datalims)*diff(axlims) + axlims(1);

figure;
[ClusterColor,pts]=DisplayMDSplot(DGSA.D, DGSA.Clustering);
ttl = sprintf('Clustering by response: \nt= %s days',string(tim));


title(ttl);
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
    namext = strcat(namext,num2str(DGSA.Nbcluster),'_t',string(tim));
    print(fullfile(exptdir,strcat('clustering',namext)),'-dtiff','-r300')
end

%% Compute DGSA sensitivity with time
% totims = {'2340.0','4860.0','7200.0'};
tgrid = totims(2:end);
hausdorff_mat_cat = concat_hausdorff_mats(exptdir,tgrid);

m=zeros(length(tgrid),size(ParametersNames,2));

sens = DGSA_over_time(DGSA,hausdorff_mat_cat,tgrid);
%% Plot sens with time
[senss,inds] = sort(sum(sens,2),'descend');
numplots = 5;
figure;
for k=1:numplots
    plot([0,tgrid],[0,sens(inds(k),:)],'LineWidth',2); hold on
end

legend(ParametersNames(inds(1:numplots)),'Location','EastOutside');
xlim([0 tgrid(end)]);ylim([0 6]);xlabel('Time(Day)','Fontsize',12,'Fontweight','bold'); 
ylabel('Sensitivity','Fontsize',12,'Fontweight','bold');
title(sprintf('Sensitivity with time \ntop %d most influential parameters',numplots))
set(gca,'Fontsize',12);
hold off

print(fullfile(exptdir,strcat('sens_wtime_top',num2str(numplots))),'-dtiff','-r300')

%% Plot 4 corners of MDS plot as well
rowslice= 1;

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

numplots=4;
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

cdf_MainFactor(DGSA.ParametersValues, DGSA.Clustering, DGSA.ParametersNames,{'al','hk1','hk2','hk_spatial'}); %last arg can be cell of params e.g. {'al','vka'}
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
%%
% Speicify the method to display standardized conditional effects.
DGSA.ConditionalEffects.Display.SensitivityByClusterAndBins=1; % if true, display pareto plots to visualize sensitivities by bins/clusters. 
DGSA.ConditionalEffects.Display.StandardizedSensitivity='Hplot'; % If omitted Pareto plot will be used. However, this is not recommended when there are lots of parameters

% Visualize conditional effects
[H,X] = DisplayConditionalEffects(DGSA,DGSA.ConditionalEffects.Display.StandardizedSensitivity)
%%
% 5.4 Display class condtional CDFs
cdf_ConditionalEffects('vka','hk2',DGSA,5)

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