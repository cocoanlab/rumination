%% Overview
%
% This is for 
% "A Dorsomedial Prefrontal Cortex-based Dynamic Functional Connectivity Model of Rumination"
% on Nat. Comm. 2023
%
% Dependencies
%   https://github.com/cocoanlab/CocoanCore
%   https://github.com/canlab/CanlabCore
%   https://github.com/spm/spm12

% TESTED ON
%   MATLAB version 2017b
%   Mac OS Mojave 
% 
% Details in README.md
%% Prerequisites.
cd('~/rumination'); % direction to github repo.
% make sure every path to dependencies are added.
% addpath(genpath('~/PATH/TO/COCOANLAB/'))
% addpath(genpath('~/PATH/TO/CANLAB/'))
% addpath(genpath('~/PATH/TO/SPM12/'))
addpath(genpath('functions'));

%% Fig.1
clear;clc
load('Datasets.mat');
load('dMPFCmodel.mat') % load model
mdl_weight = dMPFC_model.other_output{1};
mdl_intcpt = dMPFC_model.other_output{2};
 
X2 = Study2.X.dMPFC;    X3 = Study3.X.dMPFC;
y2 = Study2.Y.D;        y3 = Study3.Y.D;

y2fit = (X2 * mdl_weight) + mdl_intcpt;
y3fit = (X3 * mdl_weight) + mdl_intcpt;

% Fig.1b-2
figure(1);
scatter(y2, y2fit, 250, 'fill', 'MarkerFaceAlpha', 1,...
        'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', [0 0 0] ,'LineWidth', 2);
xlim([min(y2)-2, max(y2)+2]);
ylim([min(y2fit)-2, max(y2fit)+2]);
l = lsline;set(l, 'LineWidth', 3, 'Color', [0.1 0.1 0.1]);
set(gcf, 'Color', 'white')
set(gca, 'TickDir', 'out', 'LineWidth', 4, 'fontsize', 40, 'fontname', 'Helvetica')

% Fig.1b-3
figure(2);
scatter(y3, y3fit, 250, 'fill', 'MarkerFaceAlpha', 1,...
        'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', [0 0 0] ,'LineWidth', 2);
xlim([min(y3)-2, max(y3)+2]);
ylim([min(y3fit)-2, max(y3fit)+2]);
l = lsline;set(l, 'LineWidth', 3, 'Color', [0.1 0.1 0.1]);
set(gcf, 'Color', 'white')
set(gca, 'TickDir', 'out', 'LineWidth', 4, 'fontsize', 40, 'fontname', 'Helvetica')

%% Fig.2
clear;clc
load('Datasets.mat');
load('dMPFCmodel.mat') % load model
wobj = dMPFC_model.weight_obj;
temp = fmri_data(which('Fan_et_al_atlas_r280.nii'));

for i = 1:numel(wobj.dat)
    temp.dat(temp.dat == i) = wobj.dat(i);
end

pos_wobj = temp;
neg_wobj = temp;

pos_wobj.dat = pos_wobj.dat .* (pos_wobj.dat > 0);
neg_wobj.dat = neg_wobj.dat .* (neg_wobj.dat < 0);   

% Fig.2a
orthviews(pos_wobj) % same with "positive_weights.nii "
orthviews(neg_wobj) % same with "negative_weights.nii "

%% Fig.3
% Virtual lesion analysis.
clear;clc;
load('Datasets.mat')           % load model
load('dMPFCmodel.mat')           % load model
load('cluster_Fan_Net_r280.mat') % load parcel info.

mdl_weight  = dMPFC_model.other_output{1};
nonzeroidx  = find(mdl_weight);
weights84   = mdl_weight(nonzeroidx);
Fan_names   = cluster_Fan_Net.names;
Fan_names84 = Fan_names(nonzeroidx, :);

X2 = Study2.X.dMPFC;    X3 = Study3.X.dMPFC;
y2 = Study2.Y.D;        y3 = Study3.Y.D;

orig_y2r = corr(X2 * mdl_weight, y2); 
orig_y3r = corr(X3 * mdl_weight, y3);

for i = 1:numel(nonzeroidx)
    lesion_weights = mdl_weight;
    lesion_weights(nonzeroidx(i), :) = 0;
    
    lesion_y2fit = Study2.X.dMPFC * lesion_weights;
    study2_delta_corr(i, :) = (orig_y2r - corr(lesion_y2fit, y2));
    
    lesion_y3fit = Study3.X.dMPFC * lesion_weights;
    study3_delta_corr(i, :) = (orig_y3r - corr(lesion_y3fit, y3));
end

% Results for Fig. 3a 
corr(study2_delta_corr, study3_delta_corr);
% 0.6283

both_numidx = (study2_delta_corr > 0) & (study3_delta_corr > 0);
Fan_names21 = Fan_names84(both_numidx, :);
weights21   = round(weights84(both_numidx, :), 4);

study2_delta_corr21 = round(study2_delta_corr(both_numidx, :), 4);
study3_delta_corr21 = round(study3_delta_corr(both_numidx, :), 4);
[~, sortidx] = sort(round(weights84(both_numidx, :), 4));

% Results for Fig. 3b
T= table;
T.('regions') = Fan_names21(sortidx, :);
T.('orig_weights') = weights21(sortidx, :);
T.('Study2_delta_r') = study2_delta_corr21(sortidx, :);
T.('Study3_delta_r') = study3_delta_corr21(sortidx, :);
disp(T)


%% Fig.4
clear;clc
load('dMPFCmodel.mat') % load model
load('summaries.mat') % load model

idx84       = summaries.beta ~= 0;
weights84   = summaries.beta(idx84);
LSN84       = summaries.large_scale_network_brainnetomeprovide(idx84);
names84     = summaries.region_names(idx84);
meanDCC84   = summaries.mean_DCC(idx84);
idx21       = (summaries.delta_r_study2 > 0) & (summaries.delta_r_study3 > 0);
idx21_84    = idx21(idx84);

[LSN84, sortidx] = sort(LSN84);
[~, ~, gidx] = unique(LSN84);
gidx = [1;gidx+1];
weights84 = weights84(sortidx);
names84 = names84(sortidx);
meanDCC84 = meanDCC84(sortidx);
idx21_84 = idx21_84(sortidx);

names85 = [{'dMPFC'};names84];
names21 = names84(idx21_84);
posidx = find(meanDCC84 > 0);
posDCC = meanDCC84(posidx);
negidx = find(meanDCC84 < 0);
negDCC = meanDCC84(negidx);

poscmap = cbrewer('seq', 'Reds', 256);
negcmap = cbrewer('seq', 'Blues', 256);
poscmap = poscmap(round(rescale(posDCC, 1, 256)), :);
negcmap = negcmap(round(rescale(negDCC, 1, 256)), :);

ctemplate = zeros(84, 3);
ctemplate(posidx, :) = poscmap;
ctemplate(negidx, :) = negcmap;
cmap_patch = [0 0 0;ctemplate];

Amtrx = [0 weights84';zeros(84, 85)];

texcols = zeros(84, 3);
texcols(~idx21_84, :) = repmat([0.8 0.8 0.8], sum(~idx21_84), 1);
texcols = [0 0 0;texcols];

secondpatchcol = [166,206,227
                31,120,180
                178,223,138
                51,160,44
                251,154,153
                227,26,28
                253,191,111
                255,127,0
                202,178,214
                106,61,154
                191,91,23]./255;

cidx = num2logidx_ycgosu(gidx);
cols = zeros(85, 3);
for i = 1:size(cidx, 2)
    cols(find(cidx(:, i)), :) = repmat(secondpatchcol(i, :), numel(find(cidx(:, i))), 1);
end
cols(1, :) = [0 0 0];
idx21_85 = [true;idx21_84];
gidx2 = gidx(idx21_85);
emptyg = find(~histcounts(gidx2));
gidx2(emptyg(1) < gidx2 & emptyg(2) > gidx2) = gidx2(emptyg(1) < gidx2 & emptyg(2) > gidx2) - 1;
gidx2(emptyg(2) < gidx2) = gidx2(emptyg(2) < gidx2)-2;
cols2 = cols(idx21_85, :);
cmap_patch2 = cmap_patch(idx21_85, :);
layerinfo2 = {'layer', gidx2 ./ (max(gidx2) + 1), 'color', cmap_patch2};
Amtrx2 = Amtrx(idx21_85, idx21_85);

circos_multilayer_yc(Amtrx2, 'group', gidx2, 'each_patch_color', cols2, 'region_names', names85(idx21_85), ...
    'rotate', 185, 'patch_edge_color', [1 1 1],...
    'conn_width', rescale(abs(weights84(idx21_84)), 1, 3), 'conn_alpha', rescale(abs(weights84(idx21_84)), 0.25, 1), ...
    'text_color', texcols(idx21_85, :), 'add_layer', layerinfo2)
set(gca, 'xlim', [-2 2], 'ylim', [-2 2])
set(gcf, 'Position', [470   301   539   496])



%% Additional notes.
%
% DCC variance is calculated by "DCC" function which results in
% temporal connectivity information across parcels.
% "https://github.com/canlab/Lindquist_Dynamic_Correlation"
%
% Model is from "predict" function in 
% "https://github.com/canlab/CanlabCore"












