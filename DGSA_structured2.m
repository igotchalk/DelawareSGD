

%% DGSA: Modified from main_DGSA_Reservoir_Sensitivity.m

% Part1: Compute main/conditional effects of Prior models

clear all; close all; fclose('all'); rng('default');
%% 0. Add directories for path
if ispc()
    addpath(genpath('E:\MATLAB\Packages\DGSA-master_altered'));
    addpath(genpath('E:\Projects\DelawareSGD'));
else
    addpath(genpath('/Users/ianpg/Documents/ProjectsLocal/DGSA-master_altered'))
    addpath(genpath('/Users/ianpg/Documents/ProjectsLocal/DelawareSGD'))
end

%% 1. Specify inputs for DGSA. 

%load hausdorff matrix
% exptdir = '/Users/ianpg/Dropbox/TempShared/MC_expt_2019-01-08-02-50';
% exptdir ='E:\Projects\DelawareSGD\work\homogenous\MC_expt_2018-12-21-23-17'; %homogenous
%exptdir = 'E:\Projects\DelawareSGD\work\homogenous\MC_expt_2019-01-08-02-50'; %homogenous
%exptdir = 'E:\Projects\DelawareSGD\work\homogenous_500it\MC_expt_2019-04-04-19-30'; %homogenous500, presented to Peter,Jef
% exptdir = 'E:\Projects\DelawareSGD\work\heterogenous_500it\MC_expt_2019-04-04-21-09'; %heterogenous500 presented to peter, Jef
%exptdir = 'E:\Projects\DelawareSGD\work\passive_active\MC_expt_2019-02-12-11-55';
%exptdir = 'E:\Projects\DelawareSGD\work\passive_active\MC_expt_2019-02-18-17-52';
%exptdir ='E:/Projects/DelawareSGD/work/passive_active/MC_expt_2019-02-14-14-26'; %200 it, heterogenous, changed head_inland_sum
%exptdir ='E:/Projects/DelawareSGD/work/passive_active/MC_expt_2019-02-18-17-52/'; %500it,heterogenous,riv_cond updated (was way too high) 
%exptdir = 'E:\Projects\DelawareSGD\work\mps\MC_expt_2019-03-02-18-49\'; %500it, BC prior should be identical to MC_expt_2019-02-18-17-52


% exptdir = "E:\Projects\DelawareSGD\work\mps\MC_expt_2019-03-02-18-49\conc_mat.mat";
exptdir = '/Users/ianpg/Documents/ProjectsLocal/SWIsmall/work/heterog_1000/MC_expt_2020-02-01-18-15'


addpath(exptdir);

distmat = 'hdorf';
printplotsyn = true;


%Load conc_mat structure
conc_mat = load(fullfile(exptdir,'conc_mat.mat'));
varlist_full  = load(fullfile(exptdir,'varlist_final.mat'));
% CF_mat = load(fullfile(exptdir,'CF_mat_remake.mat'));
CF_mat = load(fullfile(exptdir,'hk_mat_remake.mat'));

% conc_mat_names = sort(fieldnames(S));
% culled_conc_mat =  conc_mat.(string(conc_mat_names(1))); %load an example
% culled_conc_mat = culled_conc_mat(1,:,:,:);
totims= conc_mat.times; 
conc_mat = rmfield(conc_mat,'times');


varlist_full_names = fieldnames(varlist_full);
varlist = varlist_full;

add_vars = {'wel0','wel1','wel2','wel3'};
for i=1:length(add_vars)
    varlist.(string(add_vars(i))) = varlist.wel(i,:);
end

rm_vars = {'it','success','seed','hk_mean','wel'};
for i=1:length(rm_vars)
    varlist =rmfield(varlist,rm_vars(i));
end
varlist_names = fieldnames(varlist);


% %Sort the conc mats
% conc_mat_names = fieldnames(conc_mat);
% its_unsorted = [];
% for i=1:length(conc_mat_names)
%     nam = conc_mat_names{i};
%     its_unsorted = [its_unsorted,str2num(nam(9:end))];
% end
% [~,its_ind] = sort(its_unsorted);


conc_mat_names = fieldnames(conc_mat);
CF_mat_names = fieldnames(CF_mat);
totims = [1800, 3600, 5400, 7200];
[its_unsorted_conc,its_ind_conc] = get_it_from_struct(conc_mat,'conc_mat');
[its_unsorted_CF,its_ind_CF] = get_it_from_struct(CF_mat,'CF');


%Remove any members of CF_mat that aren't in conc_mat
rm_CF = [];
for i=1:length(its_unsorted_CF)
    if ~ismember(its_unsorted_CF(i),its_unsorted_conc)
        rm_CF = [its_unsorted_CF(i),rm_CF];
    end
end

for i=1:length(rm_CF)
    CF_mat = rmfield(CF_mat,strcat('CF',string(rm_CF(i))));
end

%Get the new field names and inds
CF_mat_names = fieldnames(CF_mat);
[its_unsorted_CF,its_ind_CF] = get_it_from_struct(CF_mat,'CF');

comp_inds = its_unsorted_CF(its_ind_CF)==its_unsorted_conc(its_ind_conc);
if ~all(comp_inds)
    fprintf('Indicies dont line up!!! Index %d \n',find(comp_inds==0))
end


%Get the shape of the individual matricies
matshape = size(conc_mat.(string(conc_mat_names(1))));
matshape = matshape(2:end);

%Pull out all conc_mats from a specific time
tim = length(totims);
conc_mat_equitemp = zeros(numel(conc_mat_names),...
    matshape(1),matshape(2),matshape(3));
CF_mat_equitemp = zeros(size(conc_mat_equitemp));

for i=1:length(conc_mat_names)
    conc_mat_it = conc_mat.(string(conc_mat_names(its_ind_conc(i))));
    conc_mat_equitemp(i,:,:,:) = conc_mat_it(tim,:,:,:);
    CF_mat_equitemp(i,:,:,:) = CF_mat.(string(CF_mat_names(its_ind_CF(i))));
end
DGSA={};

%%
% ParametersValues = inputStruct.ParametersValues

hausdorff_mat = load(fullfile(exptdir,'hdorf_conc.mat'));
hausdorff_mat = hausdorff_mat.D;




%% Check that conc_mat and CF_mat are lined up
figure;
LU = [.05,.95];

conc_mat_temp = conc_mat_equitemp;
conc_mat_temp((conc_mat_temp < LU(1)*35) | (conc_mat_temp > LU(2)*35))=nan;

for i=1:8
    subplot(4,4,i)
    imagesc(log10(squeeze(CF_mat_equitemp(i*10,:,10,:))))
    caxis([1 2.5]);
    subplot(4,4,i+8)
    imagesc(squeeze(conc_mat_temp(i*10,:,10,:)))


    
    
end

%%
% inputStruct = load(fullfile(exptdir,'InputParams_filt.mat'));
InputParams = varlist;
ParametersNames = names2latex(fieldnames(InputParams)');
ParametersValues = zeros(length(its_unsorted_conc),length(varlist_names));
for i=1:length(varlist_names)
    tmp = varlist.(string(varlist_names(i)));
    ParametersValues(:,i) = tmp(its_unsorted_conc+1);
end

% %Do KPCA on the conc_mats
% % LU = [.05,.95];
 LU = [.05,.95];
conc_mat_temp = conc_mat_equitemp;
conc_mat_temp((conc_mat_temp < LU(1)*35) | (conc_mat_temp > LU(2)*35))=0;
conc_flat = reshape(conc_mat_temp,size(conc_mat_temp,1),[])';
[Y, eigVector,eigValue]=kpca_process(conc_flat','poly',1);

%plot the percetange of variance
ndim=ploteigvalue(eigValue,0.8);
figure;
scatter3(Y(:,1),Y(:,2),Y(:,3))
hausdorff_mat = squareform(pdist(real(Y)));    

%Assign some variables
D = hausdorff_mat;
% D = hausdorff_mat/max(max(hausdorff_mat));
N = length(D);

DGSA.D = D;
DGSA.N = N;


%% Load spatial matrices and calculate rank
%run the KPCASOM_rank script
%we used Gaussian radial basis function kernel
%For Gaussian RBF, kernel_para is for setting up the sigma value in RBF function 
% sigma= kernel_para * mean distance between models 


CF_flat = reshape(CF_mat_equitemp,size(CF_mat_equitemp,1),[])';
kernel_para=7;
[Y_eig,ndim,rank]=KPCASOM_rank(CF_flat,'gaussian',kernel_para);
%Add spatial heterogeneity to list if relevant 
ParametersNames{end+1} = 'CF_{spatial}';
ParametersValues = horzcat(ParametersValues,rank');


DGSA.ParametersValues = ParametersValues;
DGSA.ParametersNames = ParametersNames;





%% 4. Cluster and Compute Main Effects

% 4.1 Inputs for clustering and display options.
DGSA.Nbcluster=10; % # of clusters
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
med_mats = CF_mat_equitemp(medoids,:,:,:);
% med_mats = conc_mat_equitemp(medoids,:,:,:);
UL = [0.05 .95]; %Pct seawater that classifies the transition zone
% med_mats((med_mats < UL(1)*35) | (med_mats > UL(2)*35))=NaN;

med_mats(med_mats>=1e30)=NaN;
%%
%Plot MDS results and representative plots
plotreps=true;
rowslice= 10;
printplotsyn=false;

data2pos = @(data,datalims,axlims) (data-datalims(1))/diff(datalims)*diff(axlims) + axlims(1);

figure;
[ClusterColor,pts]=DisplayMDSplot(DGSA.D, DGSA.Clustering);
ttl = sprintf('Clustering by response: \nt= %s days',string(totims(tim)));


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
    pos = linspace(0.05,.95,DGSA.Nbcluster+1)';
    w_h = min([(pos(2)-pos(1))*.75,1-(h+.02)]);
    pos = horzcat(pos(1:end-1),repmat(h+.02,DGSA.Nbcluster,1),repmat(w_h,DGSA.Nbcluster,1),repmat(w_h,DGSA.Nbcluster,1));
       
    [~,sortind] = sort(pts.XData(medoids));
    for i=1:length(medoids)
        ax2 = axes('Position',pos(i,:),'Color','none');
        data=  squeeze(med_mats(sortind(i),:,rowslice,:));
        imagesc(data,'AlphaData',~isnan(data));
        ax2.XTickLabel='';
        ax2.YTickLabel='';
        x=[data2pos(pts.XData(medoids(sortind(i))),ax.XLim,[ax.Position(1),ax.Position(1)+ax.Position(3)]),...
            ax2.Position(1)+.5*ax2.Position(4)];
        y=[data2pos(pts.YData(medoids(sortind(i))),ax.YLim,[ax.Position(2),ax.Position(2)+ax.Position(4)]),...
            ax2.Position(2)];
        annotation('arrow',x,y,'LineStyle','--','HeadStyle','vback3')
        colormap(ax2,jet(50))
        caxis([1,2.5]);
    end
end
if printplotsyn 
    if plotreps
        namext='wreps';
    else
        namext='';
    end
    namext = strcat(namext,num2str(DGSA.Nbcluster),'_t',string(tim),'.tif');
    print(fullfile(exptdir,strcat('clusteringConc',namext)),'-dtiff','-r300')
end

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
    mat = squeeze(conc_mat_equitemp(vals(i),:,rowslice,:));
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

numplots=40;
startind=1;
every=15;
rowslice= 10;
[~,sortind] = sort(pts.XData);

figure;
ind=1;
for i=startind:every:startind+numplots*every
    mat = squeeze(conc_mat_equitemp(sortind(i),:,rowslice,:));
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
    print(fullfile(exptdir,'maineffectsConc'),'-dtiff','-r150')
end
%% 4.5 Display cdfs
figure;
cdf_MainFactor(DGSA.ParametersValues, DGSA.Clustering, DGSA.ParametersNames,{'CF_{var}','hk_{var}','corr_{lenyx}'}); %last arg can be cell of params e.g. {'al','vka'}
if printplotsyn 
%     print(fullfile(exptdir,sprintf('CDF_%s%s',),'-dtiff','-r300')
end

%% Compute DGSA sensitivity with time
% totims = {'2340.0','4860.0','7200.0'};
% tgrid = totims;
%hausdorff_mat_cat = concat_hausdorff_mats(exptdir,tgrid);
hausdorff_mat_cat = zeros(length(conc_mat_names),size(hausdorff_mat,1),size(hausdorff_mat,2));
tgrid = [];
for i=1:length(conc_mat_names)
    tgrid(i) = totims(i);
    hausdorff_mat_cat(i,:,:) = kpca_concmat(S.(conc_mat_names{i}));
end
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
DisplayConditionalEffects(DGSA,DGSA.ConditionalEffects.Display.StandardizedSensitivity)
%%
% 5.4 Display class condtional CDFs
fact1='hk_{var}';
fact2='CF_{spatial}';
classno=1;
for classno=1:4
    cdf_ConditionalEffects(fact1,fact2,DGSA,classno)
    print(fullfile(exptdir,sprintf('CDFcond_%s%s%d',fact1,fact2,classno)),'-dtiff','-r150')
end
%% Save all variables for futher application
%save('VariablesSaved/DGSA_Completed.mat');
%save('VariablesSaved/DGSA_Spatial_Completed.mat');



%% Custom functions

function [names_unsorted,sort_ind] = get_it_from_struct(S,prefix)

    %Sort the AEM mats
    S_names = fieldnames(S);
    names_unsorted = zeros(1,length(S_names));
    for i=1:length(S_names)
        nam = S_names{i};
        names_unsorted(i) = str2num(nam(length(prefix)+1:end));
    end
    [~,sort_ind] = sort(names_unsorted);
end


function D = kpca_concmat(conc_mat)
    LU = [.05,.95];
    conc_mat((conc_mat < LU(1)*35) | (conc_mat > LU(2)*35))=0;
    conc_flat = reshape(conc_mat,size(conc_mat,1),[])';
    [Y, eigVector,eigValue]=kpca_process(conc_flat','gaussian',7);
    D = squareform(pdist(Y));    
end

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