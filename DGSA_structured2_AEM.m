

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


% exptdir = '/Users/ianpg/Documents/ProjectsLocal/SWIsmall/work/heterog_1000/MC_expt_2020-02-01-18-15'
% exptdir = '/Users/ianpg/Documents/ProjectsLocal/SWIsmall/work/heterog_1000/MC_expt_2020-02-01-18-15/export_1000'
exptdir = '/Volumes/Samsung_T5/ProjectsLocal/SWIsmall/work/heterog_1000/MC_expt_2020-02-01-18-15/export_1000'


addpath(exptdir);

distmat = 'kpca';
printplotsyn = true;
AEMactive = true;

%Load conc_mat structure
conc_mat = load(fullfile(exptdir,'conc_mat.mat'));
AEM_mat = load(fullfile(exptdir,'AEM_mat.mat'));
varlist_full  = load(fullfile(exptdir,'varlist_final.mat'));
CF_mat = load(fullfile(exptdir,'CF_mat2.mat'));
% conc_mat_names = sort(fieldnames(S));
% culled_conc_mat =  conc_mat.(string(conc_mat_names(1))); %load an example
% culled_conc_mat = culled_conc_mat(1,:,:,:);
totims= conc_mat.times; 
conc_mat = rmfield(conc_mat,'times');


%Remove non-random variables
varlist_full_names = fieldnames(varlist_full);
varlist = varlist_full;

add_vars = {'wel0','wel1','wel2','wel3'};
for i=1:length(add_vars)
    varlist.(string(add_vars(i))) = varlist.wel(i,:);
end

rm_vars = {'it','success','seed','hk_mean','wel','clay_lyryn'};
for i=1:length(rm_vars)
    try
        varlist =rmfield(varlist,rm_vars(i));
    catch
        continue
    end
end
varlist_names = fieldnames(varlist);

%% 
[its_unsorted_conc,its_ind_con] = get_it_from_struct(conc_mat,'conc_mat');
[its_unsorted_AEM,its_ind_AEM] = get_it_from_struct(AEM_mat,'AEM_mat');
[its_unsorted_CF,its_ind_CF] = get_it_from_struct(CF_mat,'CF');
%%

CF_mat_names = fieldnames(CF_mat);
conc_mat_names = fieldnames(conc_mat);
AEM_mat_names = fieldnames(AEM_mat);

matshape = size(conc_mat.(string(conc_mat_names(1))));
matshape = matshape(2:end);
matshapeAEM = size(AEM_mat.(string(AEM_mat_names(1))));
matshapeAEM = matshapeAEM(2);
matshapeCF = size(CF_mat.(string(CF_mat_names(1))));
%%

tim_ind = length(totims);
tim = totims(tim_ind);

% Align indicies of conc_mat, CF_mat, AEM_mat

%Remove any members of CF_mat that aren't in AEM_mat
rm_CF = [];
for i=1:length(its_unsorted_CF)
    if ~ismember(its_unsorted_CF(i),its_unsorted_AEM) || ...
            ~ismember(its_unsorted_CF(i),its_unsorted_conc)
        rm_CF = [its_unsorted_CF(i),rm_CF];
    end
end
%%
for i=1:length(rm_CF)
    try
        CF_mat = rmfield(CF_mat,strcat('CF',string(rm_CF(i))));
    catch
        fprintf('it %d already not in CF_mat\n',rm_CF(i))
    end
    
    try
        conc_mat = rmfield(conc_mat,strcat('conc_mat',string(rm_CF(i))));
    catch
        fprintf('it %d already not in conc_mat\n',rm_CF(i))
    end
    
    try
        AEM_mat = rmfield(AEM_mat,strcat('AEM_mat',string(rm_CF(i))));
    catch
        fprintf('it %d already not in AEM_mat\n',rm_CF(i))
    end
end

%%
%Get the new field names and inds
[its_unsorted_conc,its_ind_conc] = get_it_from_struct(conc_mat,'conc_mat');
[its_unsorted_AEM,its_ind_AEM] = get_it_from_struct(AEM_mat,'AEM_mat');
[its_unsorted_CF,its_ind_CF] = get_it_from_struct(CF_mat,'CF');

%%
CF_mat_names = fieldnames(CF_mat);
conc_mat_names = fieldnames(conc_mat);
AEM_mat_names = fieldnames(AEM_mat);

comp_inds = its_unsorted_CF(its_ind_CF)==its_unsorted_AEM(its_ind_AEM) & ...
    its_unsorted_CF(its_ind_CF)==its_unsorted_conc(its_ind_conc) & ...
    its_unsorted_conc(its_ind_conc)==its_unsorted_AEM(its_ind_AEM) ;
if ~all(comp_inds)
    fprintf('Indicies dont line up!!! Index %d \n',find(comp_inds==0))
else
    fprintf('All indicies line up for conc_mat, AEM_mat, CF_mat!\n')
end

%% Pull out data for either 0-500 or 501-1000
export_csvs=true;
export_DGSA=true;
savehists = false;
load_rank=false;
norm_ext ='_norm';

% iter=2;
% tim_ind=4;
exptdir_DGSA =  fullfile(exptdir,'DGSA');

for tim_ind=1:4
    for iter=1:4
%         if iter~=4
%             continue
%         end
%         if tim_ind~=4
%             continue
%         end
    tim_ind
    iter
    if iter==1
        section_name='0_500';
        test_var = 'conc';
    elseif iter==2
        section_name='0_500';
        test_var = 'AEM';
    elseif iter==3
        section_name='501_1000';
        test_var = 'conc';
    elseif iter==4
        section_name='501_1000';
        test_var = 'AEM';
    end        



    if strcmp(section_name,'0_500')
        section = its_unsorted_AEM(its_unsorted_AEM <=500);
    elseif strcmp(section_name,'501_1000')
        section = its_unsorted_AEM(its_unsorted_AEM > 500);
    end
    [~,section_ind] = sort(section);

    AEM_mat_equitemp = zeros(length(section),matshapeAEM);
    CF_mat_equitemp = zeros(length(section),matshapeCF(1),matshapeCF(2),matshapeCF(3));
    conc_mat_equitemp = zeros(length(section),matshape(1),matshape(2),matshape(3));


    for i=1:length(section)
    %     fprintf('section %d\n',section_AEM(section_ind_AEM(i)))
        AEM_mat_it = AEM_mat.(strcat('AEM_mat',string(section(section_ind(i)))));
        conc_mat_it = conc_mat.(strcat('conc_mat',string(section(section_ind(i)))));

        if any(any(AEM_mat_it>0))
            fprintf('it %d has negative values\n',i)
        end
        AEM_mat_equitemp(i,:) = log10(abs(AEM_mat_it(tim_ind,:)));
        CF_mat_equitemp(i,:,:,:) = CF_mat.(strcat('CF',string(section(section_ind(i)))));
        conc_mat_equitemp(i,:,:,:) = conc_mat_it(tim_ind,:,:,:);
    end
    fprintf('Done!\n')




    if load_rank
        rank = readmatrix(fullfile(exptdir,sprintf('rank_%s.csv',section_name)))';
        norm_ext=''; %set to blank in case it wasnt done earlier
    else
        % Load spatial matrices and calculate rank
        CF_flat = reshape(CF_mat_equitemp,size(CF_mat_equitemp,1),[])';
        
        if strcmp(norm_ext,'')
            kernel_para=3; %25,17--> 15,11 <--
            [Y_eig,ndim,rank]=KPCASOM_rank(CF_flat,'gaussian',kernel_para);
        else
            CF_flat_norm = (CF_flat - repmat(mean(CF_flat,1),size(CF_flat,1),1))./repmat(std(CF_flat,1,1),size(CF_flat,1),1);
%             kernel_para=11; % <-- Use this
%             [Y_eig,ndim,rank]=KPCASOM_rank(CF_flat_norm,'gaussian',kernel_para); %<--Use this
            kernel_para=3; %25,17--> 15,11 <-- 
            [Y_eig,ndim,rank]=KPCASOM_rank(CF_flat_norm,'gaussian',kernel_para);

        end
    end
    % z=kPCA_PreImage(Y(:,end),eigVector,CF_flat_norm',para);
    % z_reconstruct = reshape(z,26,    20,   100);
    % imagesc(squeeze(z_reconstruct(:,10,:)))
    % figure;imagesc(squeeze(CF_mat_equitemp(1,:,10,:)))

    %plot the percetange of variance
%     ndim=ploteigvalue(eigValue,0.8);
%     fprintf('kernel %d ndim %d',kernel_para,ndim)
% 
%     figure;
%     scatter3(Y(:,1),Y(:,2),Y(:,3))

    % 
    % 
    figure;
    hold on
    sorted_rank = sort(rank);
    for i=1:16
    %     if ismember(i,1:4)
    %     ind  = i;
    %     elseif ismember(i,5:8)
    %     ind  = i + 100;
    %     elseif ismember(i,9:12)
    %     ind  = i + 200;
    %     elseif ismember(i,13:16)
    %     ind  = size(CF_mat_equitemp,1) - 5 -13 + i;
    %     end
        ind = i;
        subplot(4,4,i)
        imagesc(1:30:100*30,1:3:26*3, squeeze(CF_mat_equitemp(sorted_rank(ind),:,10,:)))
%         caxis([-.5,1.5])
        title(num2str(ind))
    %     axis equal
        daspect([10 1 1])
    end

    figure;
    hold on
    sorted_rank = sort(rank);
    for i=1:16
    %     if ismember(i,1:4)
    %     ind  = i;
    %     elseif ismember(i,5:8)
    %     ind  = i + 100;
    %     elseif ismember(i,9:12)
    %     ind  = i + 200;
    %     elseif ismember(i,13:16)
    %     ind  = size(CF_mat_equitemp,1) - 5 -13 + i;
    %     end
        ind = 499-17+i;
        subplot(4,4,i)
        imagesc(1:30:100*30,1:3:26*3, squeeze(CF_mat_equitemp(sorted_rank(ind),:,10,:)))
%         caxis([-.5,1.5])
        title(num2str(ind))
    %     axis equal
        daspect([10 1 1])
    end


    if strcmp(test_var,'conc')
        fprintf('Clustering on concentration data!')
        %CONCENTRATION
        LU = [.05,.95];
    %     LU = [.001,.99];
        conc_flat = reshape(conc_mat_equitemp,size(conc_mat_equitemp,1),[])';
        conc_flat((conc_flat < LU(1)*35) | (conc_flat > LU(2)*35))=0;
        [Y, eigVector,eigValue]=kpca_process(conc_flat','poly',3);
    elseif strcmp(test_var,'AEM')
        fprintf('Clustering on AEM data!')

        conc_flat = reshape(AEM_mat_equitemp,size(AEM_mat_equitemp,1),[])';
        [Y, eigVector,eigValue]=kpca_process(conc_flat','gaussian',7);
    end

    %plot the percetange of variance
    ndim=ploteigvalue(eigValue,0.8);
    figure;
    scatter3(Y(:,1),Y(:,2),Y(:,3))
    hausdorff_mat = squareform(pdist(Y));    


    % Create InputParams
    ParametersNames = names2latex(fieldnames(varlist)');
    ParametersValues = zeros(length(section),length(varlist_names));
    for i=1:length(varlist_names)
        tmp = varlist.(string(varlist_names(i)));
        ParametersValues(:,i) = tmp(section(section_ind)+1);
    end

    %Assign some variables
    D = hausdorff_mat;
    N = length(D);
    ParametersNames{end+1} = 'CF_{spatial}';
    ParametersValues = horzcat(ParametersValues,rank');


    DGSA={};
    DGSA.D = D;
    DGSA.N = N;
    DGSA.ParametersValues = ParametersValues;
    DGSA.ParametersNames = ParametersNames;



    % 4. Cluster and Compute Main Effects

    % 4.1 Inputs for clustering and display options.
    DGSA.Nbcluster=5; % # of clusters
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


    if strcmp(test_var,'conc')
        med_mats = conc_mat_equitemp(medoids,:,:,:);
        med_mats(med_mats>=1e30)=NaN;
    elseif strcmp(test_var,'AEM')
        med_mats = AEM_mat_equitemp(medoids,:);
        % med_mats(med_mats>=1e30)=NaN;
    end

    % Sort by x-location
    % Make the MDS Plot
    [Xd_, e_] = cmdscale(DGSA.D);
    % Reduction of the dimension of the MDS space, here 2D
    dims = 2;
    Xd = Xd_(:,1:dims);
    [~,sortind] = sort(Xd(DGSA.Clustering.medoids,1));


    newT = zeros(size(DGSA.Clustering.T));

    for i=1:length(DGSA.Clustering.medoids)
        newT(DGSA.Clustering.T==sortind(i))=i;
    end
    newMedoids=DGSA.Clustering.medoids(sortind);
    newWeights=DGSA.Clustering.weights(sortind);
    % newColor= DGSA.Clustering.ClusterColor(sortind,:);

    DGSA.Clustering.T=newT;
    DGSA.Clustering.labels=newT;
    DGSA.Clustering.medoids=newMedoids;
    DGSA.Clustering.weights=newWeights;
    % DGSA.Clustering.ClusterColor=newColor;


    times_hm_312 = [6.646000e-05, 7.347000e-05, 8.246000e-05, 9.396000e-05,...
           1.084600e-04, 1.264600e-04, 1.489500e-04, 1.769500e-04, ...
           2.124500e-04, 2.579500e-04, 3.149500e-04, 3.869500e-04, ...
           4.779500e-04, 5.919500e-04, 7.359500e-04, 9.174500e-04, ...
           1.146250e-03, 1.434250e-03, 1.798250e-03, 2.257250e-03, ...
           2.835250e-03, 3.554250e-03, 4.434250e-03, 5.511250e-03, ...
           6.829250e-03, 8.443250e-03, 1.041825e-02];


    times_lm_312 = [1.09700e-05, 1.79800e-05, 2.69700e-05, 3.84700e-05, 5.29700e-05, ...
           7.09700e-05, 9.34600e-05, 1.21460e-04, 1.56960e-04, 2.02460e-04, ...
           2.59460e-04, 3.31460e-04, 4.22460e-04, 5.36460e-04, 6.80460e-04, ...
           8.61960e-04, 1.09076e-03, 1.37876e-03];

    % Write matricies for R
    tim_ext = sprintf('_%d',tim_ind);

    if export_csvs
        writematrix(DGSA.Clustering.T',fullfile(exptdir_DGSA,sprintf('Cluster_%s_%s%s.csv',test_var,section_name,tim_ext)));
        writecell([DGSA.ParametersNames;num2cell(DGSA.ParametersValues)],fullfile(exptdir_DGSA,sprintf('ParametersValues_%s_%s%s%s.csv',test_var,section_name,tim_ext,norm_ext)));
        writematrix(DGSA.D,fullfile(exptdir_DGSA,sprintf('D_%s_%s%s.csv',test_var,section_name,tim_ext)));
    end

    if export_DGSA
        save(fullfile(exptdir_DGSA,sprintf('DGSA_%s_%s%s%s.mat',test_var,section_name,tim_ext,norm_ext)),'DGSA');
    end



    C= [31, 119, 180; 255, 127, 14; 44, 160, 44; 214, 39, 40; 148, 103, 189]/255;
    %Plot hists
    Tconc = readmatrix(fullfile(exptdir,sprintf('Cluster_conc_%s.csv',section_name)));
    DGSA.Clustering.Tconc = Tconc';
    set(0, 'DefaultLineLineWidth', 1);
    set(0,'defaultAxesFontSize',11)
    for i=1:DGSA.Nbcluster
        fig = figure('units','inch','position',[0,0,1.4,1.25],'color','w');

        %background data
        cluster_tot= zeros(5,1);
        cross_vals= zeros(5,1);
        hist_data = [];

        for j=1:5
            cluster_tot(j) = sum(DGSA.Clustering.T==j);
            cross_vals(j) = sum(DGSA.Clustering.T(DGSA.Clustering.Tconc==i)==j);
            hist_data = horzcat(hist_data, j*ones(1,floor(cross_vals(j)/cluster_tot(j)*1000)));
        end

        histogram(hist_data,'FaceColor',C(i,:),'LineWidth',1,'Normalization','Probability')
    %         histogram(DGSA.Clustering.T(DGSA.Clustering.Tconc==i),'FaceColor',C(i,:),'LineWidth',1)
        xticks(1:DGSA.Nbcluster)
        title(sprintf('%d',i),'FontSize',12)
        xlabel('AEM data cluster','FontSize',11)
        ylabel('Prob.','FontSize',11);
        ylim([0,.6]);
        set(gca,'fontname','Corbel')  % shows you what you are using.
        set(gca,'LineWidth',1)
        if savehists
            print(fullfile(exptdir,sprintf('clusteringConcHist_%s_%d',section_name,i)),'-dpng','-r500')
        end
    end

    close all
    end
end


%%
%Plot MDS results and representative plots
plotreps=true;
plotconc=false;
rowslice= 20;
set(0, 'DefaultLineLineWidth', 1);
UL = [0.05 .95]; %Pct seawater that classifies the transition zone
printplotsyn=false;

data2pos = @(data,datalims,axlims) (data-datalims(1))/diff(datalims)*diff(axlims) + axlims(1);
if plotreps
    fig = figure('units','inch','position',[0,0,7.48,4.5],'color','w');
else
    fig = figure('units','inch','position',[0,0,7.48,3],'color','w');
end


% Make the MDS Plot
% [Xd_, e_] = cmdscale(DGSA.D);
% % Reduction of the dimension of the MDS space, here 2D
% dims = 2;
% Xd = Xd_(:,1:dims);
% DGSA.Clustering.label=DGSA.Clustering.T;
%plotcmdmap
% C  = fake_parula(100);
% C = C(floor(linspace(1,size(C,1),max(DGSA.Clustering.T))),:);
colors = C(DGSA.Clustering.T,:);
colors_conc = C(DGSA.Clustering.Tconc,:);
% colors=DGSA.Clustering.T;
%figure
if plotconc
    pts = scatter(Xd(:,1), Xd(:,2),20,colors_conc,'filled','LineWidth', 1.3);
else
    pts = scatter(Xd(:,1), Xd(:,2),20,colors,'filled','LineWidth', 1.3);
end
hold on
if plotreps
    for i_med=1:DGSA.Nbcluster
        x=Xd(DGSA.Clustering.medoids(i_med),1);
        y=Xd(DGSA.Clustering.medoids(i_med),2);
        scatter(x,y,70,C(i_med,:),...
            'LineWidth', 1,'MarkerFaceColor','flat', ...
            'MarkerEdgeColor','k','Marker','diamond');
    %     text(x+.2,y+.2,sprintf('%d',i_med),'FontWeight','bold')
    end
end
ax = gca;
ax.XTickLabel = [];
ax.YTickLabel = [];
ax.LineWidth=1;
% 
% set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
% set(gca, 'XTickLabelMode', 'manual', 'YTickLabel', []);

set(gca, 'color',[.9,.9,.9,1]);
box on
xlabTxt=['Variance explained:',num2str(round(e_(1)/sum(e_)*100)),'%'];
ylabTxt=['Variance explained:',num2str(round(sum(e_(2))/sum(e_)*100)),'%'];    
xlabel(xlabTxt,'Fontsize',13); ylabel(ylabTxt,'Fontsize',13);
%



DGSA.Clustering.ClusterColor = C;
set(gca,'fontname','Corbel')  % shows you what you are using.

if plotreps
    h=.735;
else
    h=.95;
end
ax.OuterPosition = [0,0,1,h];
ax_pos = ax.Position;
pts.MarkerFaceColor = 'flat';
pts.MarkerEdgeColor = 'flat';
itime = 30;

if plotreps
    hold on
    pos = linspace(0.085,.97,DGSA.Nbcluster+1)';
    w_h = min([(pos(2)-pos(1))*.75,1-(h+.02)]);
    if strcmp(test_var,'conc')
        pos = horzcat(pos(1:end-1),repmat(h+.02,DGSA.Nbcluster,1),repmat(w_h,DGSA.Nbcluster,1),repmat(w_h,DGSA.Nbcluster,1));
    elseif strcmp(test_var,'AEM')
        pos = horzcat(pos(1:end-1),repmat(h+.02,DGSA.Nbcluster,1),repmat(w_h,DGSA.Nbcluster,1),repmat(w_h,DGSA.Nbcluster,1));
    end   
    [~,sortind] = sort(pts.XData(medoids));
    for i=1:length(medoids)
        ax2 = axes('Position',pos(i,:));
        if strcmp(test_var,'conc')
            data=  squeeze(med_mats(sortind(i),:,rowslice,:));
            data((data < UL(1)*35) | (data > UL(2)*35))=NaN;
            imagesc(data,'AlphaData',~isnan(data));
%             if i==1
                ylabel('Elev.','FontSize',12)
%             end
            xlabel('Dist.                 ')

        elseif strcmp(test_var,'AEM')
            data = med_mats(sortind(i),:);
            data = reshape(data,numel(times_hm_312)+numel(times_lm_312), []);
            DATA_HM = data(1:numel(times_hm_312),:);
            DATA_LM = data(numel(times_hm_312)+1:end,:);
            for j=itime-20:itime+20
                hold on;
                plot(times_hm_312,DATA_HM(:,j),'Color',[0,0,0,0.2],'LineWidth',.5);
                plot(times_lm_312,DATA_LM(:,j),'Color',[0,0,0,0.2],'LineWidth',.5);

            end
            set(ax2,'Xscale','log')
            box on
            if i==1
    %             xlabel('time')
                ylabel('db/dt','FontSize',12)
            end
        end

        set(gca,'fontname','Corbel');
        ax2.LineWidth=1;
        ax2.XTickLabel='';
        ax2.YTickLabel='';
        x=[data2pos(pts.XData(medoids(sortind(i))),ax.XLim,[ax.Position(1),ax.Position(1)+ax.Position(3)]),...
            ax2.Position(1)+.5*ax2.Position(4)];
        y=[data2pos(pts.YData(medoids(sortind(i))),ax.YLim,[ax.Position(2),ax.Position(2)+ax.Position(4)]),...
            ax2.Position(2)];
        annotation('arrow',x,y,'LineStyle','--','Linewidth',1,'HeadStyle','vback2')
        colormap(ax2,jet(50))
        set(ax2, 'color',[.9,.9,.9,1]);
        title(sprintf('%d',i),'FontSize',12)
    end
end

if strcmp(test_var,'conc')
    varname_disp = 'Concentration';
elseif strcmp(test_var,'AEM')
    varname_disp = 'AEM data';
end
if strcmp(section_name,'0_500')
    section_disp = 'unconfined';
elseif strcmp(section_name,'501_1000')
    section_disp = 'confined';
end
ttl = sprintf('%s, %s',varname_disp,section_disp);

axes( 'Position', [0, 0.93, 1, 0.05] ) ;
axis off
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
text( 0.5, 0, ttl, 'FontSize', 14', 'FontWeight', 'Bold', ...
  'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom','FontName','Corbel' ) ;
% suptitle(ttl);

if printplotsyn 
    if plotreps
        namext='wreps';
    else
        namext='';
    end
    namext = strcat(namext,num2str(DGSA.Nbcluster),'_t',string(tim));
    print(fullfile(exptdir,sprintf('clustering_%s_%s_%s',test_var,section_name,namext)),'-dpng','-r500')
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

numplots=40;
startind=1;
every=15;
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
    print(fullfile(exptdir,'maineffectsAEM'),'-dtiff','-r300')
end
%% 4.5 Display cdfs

cdf_MainFactor(DGSA.ParametersValues, DGSA.Clustering, DGSA.ParametersNames,{'corr_{lenyx}','CF_{spatial}'}); %last arg can be cell of params e.g. {'al','vka'}
if printplotsyn 
    %print(strcat(exptdir,'/cdfs'),'-dtiff','-r150')
end

%% Compute DGSA sensitivity with time
% totims = {'2340.0','4860.0','7200.0'};
% tgrid = totims;
%hausdorff_mat_cat = concat_hausdorff_mats(exptdir,tgrid);
hausdorff_mat_cat = zeros(length(totims),size(hausdorff_mat,1),size(hausdorff_mat,2));

for j=1:length(totims)
    for i=1:length(its_unsorted_AEM)
        AEM_mat_it = AEM_mat.(strcat('AEM_mat',string(its_unsorted_AEM(its_ind_AEM(i)))));
        AEM_mat_equitemp(i,:) = log10(abs(AEM_mat_it(j,:)));
    end
    AEM_flat = reshape(AEM_mat_equitemp,size(AEM_mat_equitemp,1),[])';
    [Y, eigVector,eigValue]=kpca_process(AEM_flat','gaussian',7);
    D = squareform(pdist(Y));    
    hausdorff_mat_cat(j,:,:) = D;
end

% 
% for j=1:length(totims)
%     for i=1:length(its_unsorted_AEM)
%         AEM_mat_it = AEM_mat.(string(its_unsorted_AEM(its_ind_AEM(i))));
%         AEM_mat_equitemp(i,:) = AEM_mat_it(tim_ind,:);
%     end
%     AEM_flat = reshape(AEM_mat_equitemp,size(AEM_mat_equitemp,1),[])';
%     [Y, eigVector,eigValue]=kpca_process(AEM_flat','gaussian',7);
%     D = squareform(pdist(Y));    
%     hausdorff_mat_cat(j,:,:) = D;
% end
tgrid = [0,totims];

m=zeros(length(tgrid),size(ParametersNames,2));

sens = DGSA_over_time(DGSA,[],tgrid,hausdorff_mat_cat);
%% Plot sens with time
[senss,inds] = sort(sum(sens,2),'descend');
numplots = 5;
figure;
for k=1:numplots
    plot(tgrid,[0,sens(inds(k),:)],'Marker','x','LineWidth',1.5); hold on
end

legend(DGSA.ParametersNames(inds(1:numplots)),'Location','EastOutside');
xlim([mean(0,tgrid(2)), tgrid(end)]);xlabel('Time(Day)','Fontsize',12,'Fontweight','bold'); 
ylim([.4,.8])
% ylim([5.5,6.6])

ylabel('Sensitivity','Fontsize',12,'Fontweight','bold');
title(sprintf('Sensitivity with time \ntop %d most influential parameters',numplots))
set(gca,'Fontsize',12);
hold off

print(fullfile(exptdir,strcat('sens_wtimeAEM_top',num2str(numplots))),'-dtiff','-r300')


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
cdf_ConditionalEffects('wel0','CF_{spatial}',DGSA,4)

%% Save all variables for futher application
save('VariablesSaved/DGSA_Completed_AEM.mat','DGSA');
% save('VariablesSaved/DGSA_Spatial_Completed.mat');



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