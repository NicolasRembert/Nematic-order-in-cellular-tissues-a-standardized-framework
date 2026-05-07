clear all; close all; clc

% ============================================================
% GLOBAL STYLE
% ============================================================
set(groot,'defaultAxesFontWeight','bold')
set(groot,'defaultTextFontWeight','bold')

% ============================================================
% PATH
% ============================================================
mainFolder = 'path to \Coherency data';
searchSubfolders = false; % Change to True if your csv are regrouped by subfolders example A02 -B02-C02 , to false if they are all in the cell type folder 

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
% COHERENCY VS TENSOR SIZE
% ============================================================
cellTypeFolders = dir(mainFolder);
cellTypeFolders = cellTypeFolders([cellTypeFolders.isdir] & ~ismember({cellTypeFolders.name},{'.','..'}));

figure('Color','w'); hold on;

% stockage pour légende
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

    % clé = tensor, valeur = vecteur de moyennes de replicates
    tensorData = containers.Map('KeyType','double','ValueType','any');

    for j = 1:length(csvFiles)

        fname = csvFiles(j).name;
        fpath = fullfile(csvFiles(j).folder, fname);

        tokens = regexp(fname, 'Tensor_(\d+)', 'tokens','once');
        if isempty(tokens), continue; end
        tensorSize = str2double(tokens{1});

        data = readmatrix(fpath);
        if size(data,2) < 7, continue; end

        coherency = data(:,7);

        % ===== MOYENNE PAR REPLICATE =====
        repMean = mean(coherency,'omitnan');

        if isKey(tensorData,tensorSize)
            tensorData(tensorSize) = [tensorData(tensorSize); repMean];
        else
            tensorData(tensorSize) = repMean;
        end
    end

    if isempty(tensorData), continue; end

    tensorSizes = sort(cell2mat(keys(tensorData)));
    meanC = zeros(size(tensorSizes));
    semC  = zeros(size(tensorSizes));

    for k = 1:length(tensorSizes)

        vals = tensorData(tensorSizes(k)); % = replicates

        meanC(k) = mean(vals,'omitnan');
        semC(k)  = std(vals,'omitnan')/sqrt(numel(vals));

        % DEBUG
        fprintf('%s | Feature size %d → nRep = %d\n', ...
            cellTypeName, tensorSizes(k), numel(vals));
    end

    key = lower(cellTypeName);
    if isKey(colorMap,key)
        col = colorMap(key);
    else
        col = [0 0 0];
    end

    % ===== SHADOW SEM =====
    x = tensorSizes;
    y = meanC;
    e = semC;

    fill([x fliplr(x)], ...
         [y+e fliplr(y-e)], ...
         col, ...
         'FaceAlpha',0.2, ...
         'EdgeColor','none');

    % ligne moyenne uniquement
    h = plot(x,y,'-','Color',col,'LineWidth',2);

    % stockage pour la légende
    legendHandles(end+1) = h;
    legendNames{end+1} = cellTypeName;

end

% remettre lignes au-dessus
uistack(findobj(gca,'Type','line'),'top')

% ============================================================
% AXES FIXES
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
legend(legendHandles, legendNames,'Location','eastoutside')

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
