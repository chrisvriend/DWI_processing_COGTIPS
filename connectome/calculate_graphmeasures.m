%% Extracting graph measures from structual connectivity matrices created using MRtrix3
% relies on the Brain Connectivity Toolbox (BCT)
%
% Assumptions: Input matrices are correctly pre processed. Input matrices
% are expected to be parcelated according to the BNA atlas and reshufled to
% have all left brain areas first and Right areas after.
%
% Author:  Chris Vriend
% Date: Jul 2020
% _____________________________________________________

clear

%% initialization
% Add BCT to path
addpath('/path/to/BCT')

%% input variables
datadir='/path/to/data';
outputdir=strcat(datadir,filesep,'graphmeasures');
conv_list_YEO=dlmread(strcat(datadir,filesep,'YEO_BNA_subnetwork_labels_mod2sample.txt'));
samplename='cogtips';
% BNA Atlas #roi = 224 - excluded nodes
nroi=222;

%% begin processing

if exist(outputdir,'dir')~=7
    mkdir(outputdir)
end


cd(datadir)


% select all connectivity matrices
a=dir('*_connectome.BNA_reshuffled.mat');
disp(['there are ' num2str(length(a)) ' files to process'])

% preallocation
nsubjects=length(a);
data_raw = nan(nroi, nroi, nsubjects);
data_length = nan(nroi, nroi, nsubjects);
subjects=cell(nsubjects,1);

for i = 1: length(a)

    subjects{i,1}=extractBefore(a(i).name,'_connectome');
    file=matfile(a(i).name);
    smatrix=file.connmatrix;
    if size(smatrix,1)~=nroi
        error('matrix dimension not consistent with number of nodes')

    end

    smatrix_temp=weight_conversion(smatrix,'autofix'); % autofix to ensure that matrix is truly symmetrical,
    nmatrix=weight_conversion(smatrix_temp, 'normalize');
    %  lmatrix = weight_conversion(nmatrix,'lengths');

    % Save matrices
    data_raw(:,:, i) = nmatrix;
    %    data_length(:,:,i)=lmatrix;
    clear smatrix* nmatrix
    %  clear lmatrix

end


%% Definition of subnetworks
% subnetwork delineated according to Yeo et al 2011 Journal of
% Neurophysiology 7Network solution

% list of original BNA node labels to Yeo subnetwork labels
% column 1 = BNA nodes (orig labels), column 2 = Yeo network (1-7),
% subcortical (8) or cerebellum (9)
% row index = row/col index of connectivity matrix

if size(conv_list_YEO,1)~=nroi
    error('array dimension of YEO subnetwork list not consistent with number of rois')
end

if size(conv_list_YEO,2)~=3
    conv_list_YEO=[transpose(1:length(conv_list_YEO)), conv_list_YEO];
end

dmn  = conv_list_YEO(conv_list_YEO(:,3)==7,1);
fpn  = conv_list_YEO(conv_list_YEO(:,3)==6,1);
%limb = conv_list_YEO(conv_list_YEO(:,3)==5,1);
van  = conv_list_YEO(conv_list_YEO(:,3)==4,1);
dan  = conv_list_YEO(conv_list_YEO(:,3)==3,1);
%smn  = conv_list_YEO(conv_list_YEO(:,3)==2,1);
%vis  = conv_list_YEO(conv_list_YEO(:,3)==1,1);

modular_community_structure = conv_list_YEO(:, 3);

subnetworks={dmn fpn dan van };
subnetworks_name={'DMN' 'FPN' 'DAN' 'VAN'};
% only the DMN, FPN, DAN and VAN are used below


%% Graph measures calculation
% Pre-allocation
% One per subject
efficiency_global = nan(nsubjects, 1);
q_optimal = nan(nsubjects, 1);
average_clustering = nan(nsubjects, 1);

% One per node
clustering = nan(nsubjects, nroi);
modularity = nan(nsubjects, nroi);


% Whole subnetwork measures
efficiency_global_subnetwork = nan(nsubjects, length(subnetworks));
clustering_subnetwork = nan(nsubjects, length(subnetworks));
between_subnetwork_DMN = nan(nsubjects, length(subnetworks));
between_subnetwork_FPN = nan(nsubjects, length(subnetworks));
between_subnetwork_DAN = nan(nsubjects, length(subnetworks));
between_subnetwork_VAN = nan(nsubjects, length(subnetworks));

%% CALCULATE NETWORK MEASURES / SUBJECT

% code is optimized for the parallel toolbox. use the for instead of parfor loop if not available
%for nn = 1:nsubjects
parfor nn = 1:nsubjects

    disp(nn);
    % Get subject nn data
    subjmatrix=data_raw(:,:,nn);

    %length_matrix = data_length(:,:,nn);

    %% Global Measures

    % global efficiency
    efficiency_global(nn) = efficiency_wei(subjmatrix);

    % Newman Spectral Community detection
    [modularity(nn, :), q_optimal(nn, 1)] = modularity_und(subjmatrix);

    %% Nodal measures

    % Strenght
    % strenght(nn, :) = strengths_und(subjmatrix)./(nroi-1);

    % Weighted clustering coefficient (nodal) & Average CC (global)
    clustering(nn, :) = clustering_coef_wu(subjmatrix)';
    average_clustering(nn) = mean(clustering(nn, :)); % = GLOBAL MEASURE

    clust = clustering(nn, :); % Auxiliar variable to subnetwork calculations


    % Subnetwork Measures - local variables to parfor
    eff_temp = nan(length(subnetworks), 1);
    acc_temp = nan(length(subnetworks), 1);
    bc1_temp = nan(length(subnetworks), 1);
    bc2_temp = nan(length(subnetworks), 1);
    bc3_temp = nan(length(subnetworks), 1);
    bc4_temp = nan(length(subnetworks), 1);
    %   bc5_temp = nan(length(subnetworks), 1);
    %   bc6_temp = nan(length(subnetworks), 1);
    %   bc7_temp = nan(length(subnetworks), 1);

    % For each subnetwork - whole-network measures
    for sub_nn = 1:length(subnetworks)
        subnetwork = subjmatrix(subnetworks{sub_nn}, subnetworks{sub_nn});
        %  subnetwork_lenght = length_matrix(subnetworks{sub_nn}, subnetworks{sub_nn});

        % calculate efficiency of subnetwork
        eff_temp(sub_nn, 1) = efficiency_wei(subnetwork);


        % calculate clustering coefficient
        % by averaging the nodal measures of all nodes that belong to a
        % particular subnetwork

        acc_temp(sub_nn) = mean(clust(subnetworks{sub_nn}));

        % Between Connectivity measures
        % order = 1= DMN , 2= FPN, 3= DAN, 4= VAN
        bc1_temp(sub_nn) = mean(mean(subjmatrix(subnetworks{sub_nn}, subnetworks{1})));
        bc2_temp(sub_nn) = mean(mean(subjmatrix(subnetworks{sub_nn}, subnetworks{2})));
        bc3_temp(sub_nn) = mean(mean(subjmatrix(subnetworks{sub_nn}, subnetworks{3})));
        bc4_temp(sub_nn) = mean(mean(subjmatrix(subnetworks{sub_nn}, subnetworks{4})));
        %  bc5_temp(sub_nn) = mean(mean(subjmatrix(subnetworks{sub_nn}, subnetworks{5})));
        %  bc6_temp(sub_nn) = mean(mean(subjmatrix(subnetworks{sub_nn}, subnetworks{6})));
        %  bc7_temp(sub_nn) = mean(mean(subjmatrix(subnetworks{sub_nn}, subnetworks{7})));
    end


    efficiency_global_subnetwork(nn, :) = eff_temp;
    clustering_subnetwork(nn, :) = acc_temp;

    between_subnetwork_DMN(nn, :) = bc1_temp;
    between_subnetwork_FPN(nn, :) = bc2_temp;
    between_subnetwork_DAN(nn, :) = bc3_temp;
    between_subnetwork_VAN(nn, :) = bc4_temp;
    %    between_subnetwork_SMN(nn, :) = bc3_temp;
    %    between_subnetwork_LIMB(nn, :) = bc6_temp;
    %    between_subnetwork_VIS(nn, :) = bc7_temp;
end

%% Save global whole-brain measures

% global measures
tab = table(subjects, efficiency_global,...
    q_optimal, average_clustering);

% Save subnetwork measures

tab1 = table(subjects,efficiency_global_subnetwork ,...
    clustering_subnetwork, between_subnetwork_DMN,...
    between_subnetwork_FPN, between_subnetwork_DAN,...
    between_subnetwork_VAN );

writetable(tab,strcat(outputdir, filesep,'graph_global_BNA_', samplename, '.csv'),'Delimiter',',','WriteVariableNames', 1);
writetable(tab1,strcat(outputdir,filesep, 'graph_subnetwork_BNA_', samplename, '.csv'),'Delimiter',',','WriteVariableNames', 1);
