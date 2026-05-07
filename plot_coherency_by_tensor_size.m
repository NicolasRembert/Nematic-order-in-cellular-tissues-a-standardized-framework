clear all; close all; clc

% ============================================================
% GLOBAL STYLE
% ============================================================
set(groot,'defaultAxesFontWeight','bold')
set(groot,'defaultTextFontWeight','bold')

% ============================================================
% PATH
% ============================================================
mainFolder = 'Path to \Coherency data';
searchSubfolders = false;

% ============================================================
% FAMILIES + COLORS
% ============================================================
families = struct( ...
    'Myoblasts',   {{'BAOSMC','C2C12','H9C2','L6','hBEC'}}, ...
    'Fibroblasts', {{'CAF','NAF','NIH3T3'}}, ...
    'Epithelial',  {{'MDCK','RPE1','U2OS'}}, ...
    'Other',       {{'SKMEL'}} ...
);

baseColors.Myoblasts   = [0 0.4470 0.7410];
baseColors.Epithelial  = [0.8500 0.3250 0.0980];
baseColors.Fibroblasts = [0.4660 0.6740 0.1880];
baseColors.Other       = [0.4940 0.1840 0.5560];

colorMap = createColorMap(families, baseColors);

% ============================================================
% COHERENCY VS FEATURE SIZE
% ============================================================
cellTypeFolders = dir(mainFolder);
cellTypeFolders = cellTypeFolders([cellTypeFolders.isdir] & ~ismember({cellTypeFolders.name},{'.','..'}));

figure('Color','w'); hold on;

legendHandles = [];
legendNames = {};

for i = 1:length(cellTypeFolders)

    cellTypeName = cellTypeFolders(i).name;
    cellTypePath = fullfile(mainFolder, cellTypeName);

    if searchSubfolders
        csvFiles = dir(fullfile(cellTypePath, '**', '*.csv'));
    else
        csvFiles = dir(fullfile(cellTypePath, '*.csv'));
    end

    % clé = feature size, valeur = vecteur de moyennes de replicates
    featureData = containers.Map('KeyType','double','ValueType','any');

    for j = 1:length(csvFiles)

        fname = csvFiles(j).name;
        fpath = fullfile(csvFiles(j).folder, fname);

        % ===== EXTRAIRE FEATURE SIZE =====
        tokens = regexp(fname, 'FeatureSize_(\d+)', 'tokens','once');
        if isempty(tokens), continue; end
        featureSize = str2double(tokens{1});

        data = readmatrix(fpath);
        if size(data,2) < 7, continue; end

        coherency = data(:,7);

        % ===== MOYENNE PAR REPLICATE =====
        repMean = mean(coherency,'omitnan');

        if isKey(featureData,featureSize)
            featureData(featureSize) = [featureData(featureSize); repMean];
        else
            featureData(featureSize) = repMean;
        end
    end

    if isempty(featureData), continue; end

    featureSizes = sort(cell2mat(keys(featureData)));
    meanC = zeros(size(featureSizes));
    semC  = zeros(size(featureSizes));

    for k = 1:length(featureSizes)

        vals = featureData(featureSizes(k));

        meanC(k) = mean(vals,'omitnan');
        semC(k)  = std(vals,'omitnan')/sqrt(numel(vals));

        % DEBUG
        fprintf('%s | FeatureSize %d → nRep = %d\n', ...
            cellTypeName, featureSizes(k), numel(vals));
    end

    key = lower(cellTypeName);
    if isKey(colorMap,key)
        col = colorMap(key);
    else
        col = [0 0 0];
    end

    % ===== SHADOW SEM =====
    x = featureSizes;
    y = meanC;
    e = semC;

    fill([x fliplr(x)], ...
         [y+e fliplr(y-e)], ...
         col, ...
         'FaceAlpha',0.2, ...
         'EdgeColor','none');

    % ===== LIGNE MOYENNE =====
    h = plot(x,y,'-','Color',col,'LineWidth',2);

    legendHandles(end+1) = h;
    legendNames{end+1} = cellTypeName;
end

uistack(findobj(gca,'Type','line'),'top')

% ============================================================
% AXES
% ============================================================
xlabel('Feature size (px)','FontSize',18)
ylabel('Average coherency','FontSize',18)

xlim([0 100])
xticks(0:25:100)

ylim([0 0.55])
yticks(linspace(0,0.55,5))

set(gca,'FontSize',18,'LineWidth',2,'Box','off')

pbaspect([1 1 1])
grid off

% ============================================================
% LEGEND
% ============================================================
legend(legendHandles, legendNames, ...
    'Location','eastoutside', ...
    'FontSize',12, ...
    'Box','off');

% ============================================================
% FUNCTION
% ============================================================
function cmap = createColorMap(families, baseColors)

    famNames = fieldnames(families);
    cmap = containers.Map();

    for f = 1:numel(famNames)

        members = families.(famNames{f});
        base = baseColors.(famNames{f});
        N = numel(members);

        tVals = linspace(0.4,0.7,N);

        for k = 1:N
            shade = base*(1-tVals(k)) + tVals(k)*[1 1 1];
            cmap(lower(members{k})) = shade;
        end
    end
end
