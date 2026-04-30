% ==============================
%  Master Script: Defect Analysis + Plotting
%  ==============================

clear; clc; close all;


addpath("Path to \Functions")
% -------------------------------
% Main directories, subfolder naming must be the same in case of different
% maindirs 
% -------------------------------
MainDir = "Path to \Cell data";          % no_corr_length files
MainDir_PIV = 'Path to \PIV Data';  % <<<< PIV directory


output_dir = "Path to \Figures velocity";

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fig_output=output_dir;

% -------------------------------
% Ignore list
% -------------------------------
ignoreCells = {'MDCK','U2OS','hBEC','SKMEL'}; % Add cells you don't want to analyse

% -------------------------------
% Families & Colors
% -------------------------------
families = struct( ...
    'Epithelial',  {{'MDCK','RPE1','U2OS'}}, ...
    'Myoblasts',   {{'BAOSMC','C2C12','H9C2','L6'}}, ...
    'Fibroblasts', {{'CAF','NAF','NIH3T3'}}, ...
    'Other',       {{'hBEC','SKMEL'}} ...
);

baseColors.Myoblasts   = [0 0.4470 0.7410];     
baseColors.Epithelial  = [0.8500 0.3250 0.0980];
baseColors.Fibroblasts = [0.4660 0.6740 0.1880];
baseColors.Other       = [0.4940 0.1840 0.5560];

colorMap = createColorMap(families, baseColors);


% ==============================
%  Phase 1: Analysis
% ==============================

Sub_maindir = dir(MainDir);
Sub_maindir = Sub_maindir([Sub_maindir.isdir]);
Sub_maindir = Sub_maindir(~ismember({Sub_maindir.name},{'.','..','outputdir'}));

S = struct();

 typeanal_list={'pos','neg'};


 for jk=1:length(typeanal_list)

     typeanal= typeanal_list{jk};


for sub = 1:length(Sub_maindir)

    cellType = Sub_maindir(sub).name;

    if ismember(cellType,ignoreCells)
        fprintf('Skipping ignored cell type: %s\n', cellType);
        continue;
    end

    

   
    % ======================== INITIALIZE CELLS ========================
    All_Magn_all = {};
    All_Div_all = {};
    All_u_ori_all = {};
    All_v_ori_all = {};
    All_u_vel_all = {};
    All_v_vel_all = {};
    All_xBox_vel_all = {};
    All_yBox_vel_all = {};
    All_xBox_ori_all = {};
    All_yBox_ori_all = {};

    All_RadialDist_Pos = {};
    All_RadialMag_Pos  = {};
    All_RadialDist_Neg = {};
    All_RadialMag_Neg  = {};

    All_NemVelAngle_all = {};

    fprintf('\n=== Processing Cell type: %s ===\n', cellType);

    subPath = fullfile(MainDir, cellType);        % no_corr_length
    subPath_PIV = fullfile(MainDir_PIV, cellType); % PIV directory

    % -------------------------------
    % Find subfolders
    % -------------------------------
    subfolders = dir(subPath);
    subfolders = subfolders([subfolders.isdir]);
    subfolders = subfolders(~ismember({subfolders.name},{'.','..'}));

    % Only keep folders like A1 B2 etc
    pattern = '^[A-Za-z]\d+';
    isValid = ~cellfun(@isempty, regexp({subfolders.name}, pattern));
    subfolders = subfolders(isValid);

    matFiles = [];
    pivFiles = [];

    % -------------------------------
    % Search files
    % -------------------------------
    for k = 1:length(subfolders)

        f = subfolders(k);

        % ---- no_corr_length ----
        files_noCorr = dir(fullfile(subPath, f.name, '*_defects.mat'));

        % ---- PIV files in second directory ----
        files_piv = dir(fullfile(subPath_PIV, f.name, '*PIVlab*batch*.mat'));

        matFiles = [matFiles; files_noCorr];
        pivFiles = [pivFiles; files_piv];

    end

    % ==============================
    % LOAD LOOP
    % ==============================

    for m = 1:length(matFiles)

        movieIdx = m;

        % ---- load no_corr_length ----
        matPath = fullfile(matFiles(m).folder, matFiles(m).name);
        fprintf('   Found corr file: %s\n', matPath);

        data = load(matPath);

        
        


        if ~isfield(data,'Qmat_list')
            fprintf('   ⚠️ Missing Qmat_list in %s → skipping.\n', matPath);
            continue;
        end

        % ---- find corresponding PIV file ----
        currentSubfolder = subfolders(m).name;

        pivInSameFolder = pivFiles(contains({pivFiles.folder}, currentSubfolder));

        if isempty(pivInSameFolder)

            fprintf('   ⚠️ No PIVlab batch file for %s → skipping.\n', currentSubfolder);
            continue;

        end

        pivPath = fullfile(pivInSameFolder(1).folder, pivInSameFolder(1).name);

        fprintf('   Found PIV file: %s\n', pivPath);

        pivData = load(pivPath);

        close all

        % -------------------------------------------------
        % Your analysis continues here
        % -------------------------------------------------

    



 

    fields = fieldnames(pivData);
        for k = 1:numel(fields)
            data.(fields{k}) = pivData.(fields{k});
        end

        % Extract variables
        u_component = data.u_component;
        v_component = data.v_component;
        divergence = data.divergence;
        velocity_magnitude = data.velocity_magnitude;
        OrientationMat_list = data.OrientationMat_list;
        OrientationMat = data.OrientationMat;
        Qparam_list = data.Qparam_list;
        Qmat_list = data.Qmat_list;
        Phi_list = data.Phi_list;

        Grid_size=data.Grid_size;
        pixelSize=data.Pixel_size;

        SizeImagePixels=length(OrientationMat)* Grid_size;
        ImageSize= SizeImagePixels^2;   
        timePerPoint=data.timePerPoint;


        % --- Prepare full cell arrays (as fournis par les .mat) ---
U_comp_list = u_component;           % cell array
V_comp_list = v_component;           % cell array
Div_list    = divergence;            % cell array
Mag_list    = velocity_magnitude;    % cell array (will multiply by factor below)
Ori_list    = OrientationMat_list;   % cell array
Qparam_list = Qparam_list;           % already present
Qmat_list   = Qmat_list;             % cell array
Phi_list_all= Phi_list;              % cell array

numFrames = length(Ori_list);
if numFrames == 0
    warning('No orientation frames found; skipping pair.');
    continue;
end

% --- Compute x_0,y_0 and scaleFactor once from first frame (assumed constant) ---
Uvalues0 = U_comp_list{1};
Orientation0 = Ori_list{1};

SizeNematic0  = length(Orientation0);
SizeVelocity0 = length(Uvalues0);

stepVel = SizeImagePixels / SizeVelocity0;
y_0 = repmat(stepVel:stepVel:SizeImagePixels,1,SizeVelocity0)'; % column vector
x_0 = repelem(stepVel:stepVel:SizeImagePixels,SizeVelocity0)';

stepNem = SizeImagePixels / SizeNematic0;
y = repmat(stepNem:stepNem:SizeImagePixels,1,SizeNematic0)';
x = repelem(stepNem:stepNem:SizeImagePixels,SizeNematic0)';

% ========================
% Nematic–Velocity difference (ALL time steps)
% ========================

NemVelAngle_list = cell(1, numFrames);


for t = 1:numFrames

    % --- Nematic orientation (on nematic grid) ---
    Ori_t = Ori_list{t};   % angle (rad)
    
    % interpolate nematic angle on velocity grid
    Ori_interp = interpolateValues(x_0(:), y_0(:), x(:), y(:), Ori_t);
    
    % reshape to velocity grid
    nxy_vel = sqrt(numel(U_comp_list{t}));
    Ori_interp = reshape(Ori_interp, nxy_vel, nxy_vel);

    % --- Nematic unit vector ---
    nx = cos(Ori_interp);
    ny = sin(Ori_interp);

    % --- Velocity field ---
    U = U_comp_list{t};
    V = V_comp_list{t};
    Vnorm = sqrt(U.^2 + V.^2);

    % avoid division by zero
    Vnorm(Vnorm == 0) = NaN;

    % --- Angle between velocity and nematic ---
    cosTheta = (U .* nx + V .* ny) ./ Vnorm;
    cosTheta = max(min(cosTheta,1),-1);  % numerical safety

    DeltaTheta = acos(cosTheta);  % [0, pi]

    % optional: nematic symmetry
    DeltaTheta = min(DeltaTheta, pi - DeltaTheta);

    % --- Projection (optional but useful) ---
    ProjVN = (U .* nx + V .* ny);

    % --- Store ---
    NemVelAngle_list{t} = DeltaTheta;
    

end



% (optional) compute interpolations for visualization / debugging on first frame
testinterp = interpolateValues(x_0,y_0,x,y,Ori_list{1});
InterpOrder = interpolateValues(x_0,y_0,x,y,Qparam_list{1});
n = ceil(sqrt(numel(testinterp)));
testinterp = reshape(testinterp,n,n);
InterpOrder = reshape(InterpOrder,n,n);

% reshape x_0,y_0 to matrices (as in original)
nxy = ceil(sqrt(numel(x_0)));
x_0 = reshape(x_0,nxy,nxy);
y_0 = reshape(y_0,nxy,nxy);

% compute scale factor (assume constant across frames)
oldSize = length(Orientation0);
newSize = length(Uvalues0);
scaleFactor = newSize / oldSize;

% ======================== IDENTIFY DEFECTS FOR ALL FRAMES ========================
% We'll store positions per frame in cell arrays so analyze_defects2 can process them
posPOSx = cell(1, numFrames);
posPOSy = cell(1, numFrames);
posNEGx = cell(1, numFrames);
posNEGy = cell(1, numFrames);

% Also compute radial magnitude per frame and accumulate in *_All lists
All_RadialDist_Pos = {};
All_RadialMag_Pos  = {};
All_RadialDist_Neg = {};
All_RadialMag_Neg  = {};

edgeThreshold = 10; % keep same threshold per-frame

% Here we will filter out the defects that are less than 10 nodes away to the edges of the image
margin = 10; % number of nodes to exclude from each border

for t = 1:length(Qmat_list)
    Q = Qmat_list{t};
    [rows, cols] = size(Q);

    % zero out defects too close to the border
    Q(1:margin,:)    = 0; % top edge
    Q(end-margin+1:end,:) = 0; % bottom edge
    Q(:,1:margin)    = 0; % left edge
    Q(:,end-margin+1:end) = 0; % right edge

    Qmat_list{t} = Q;
end

      

            for t = 1:numFrames
    Qmat = Qmat_list{t};
    % ensure Qmat is square
    matrixSize = size(Qmat,1);

    tmp_pos_x = [];
    tmp_pos_y = [];
    tmp_neg_x = [];
    tmp_neg_y = [];

    for yy = 1:matrixSize
        for xx = 1:matrixSize
            charge = Qmat(xx,yy); % follow your indexing convention
            if xx <= edgeThreshold || xx > matrixSize-edgeThreshold || ...
               yy <= edgeThreshold || yy > matrixSize-edgeThreshold
                continue;
            end
            if charge == 0.5
                tmp_pos_x(end+1) = xx; %
                tmp_pos_y(end+1) = yy; %
            elseif charge == -0.5
                tmp_neg_x(end+1) = xx; %
                tmp_neg_y(end+1) = yy; %
            end
        end
    end

    % store raw positions (in Q-matrix coordinates) for this frame
    posPOSx{t} = tmp_pos_x;
    posPOSy{t} = tmp_pos_y;
    posNEGx{t} = tmp_neg_x;
    posNEGy{t} = tmp_neg_y;

    % --- compute scaled positions for radial magnitude function (uses velocity/magnitude resolution) ---
    % Use U_comp_list{t} to get newSize for this frame (should match newSize computed above)
    Uvalues_t = U_comp_list{t};
    oldSize_t = length(Ori_list{t});
    newSize_t = length(Uvalues_t);
    scaleFactor_t = newSize_t / oldSize_t; % for safety if it differs frame-to-frame

    posPOSx_new_t = round(tmp_pos_x * scaleFactor_t);
    posPOSy_new_t = round(tmp_pos_y * scaleFactor_t);
    posNEGx_new_t = round(tmp_neg_x * scaleFactor_t);
    posNEGy_new_t = round(tmp_neg_y * scaleFactor_t);

    % --- Magnitude & Div for this frame (with unit conversion as in original) ---
    MagnitudeData_t = Mag_list{t} .* 3.6e9;
    Divdata_t      = Div_list{t};

    % compute radial magnitude for this frame and append to global lists
    if ~isempty(posPOSx_new_t)
        [distMatPos, magMatPos, ~, ~] = compute_defect_radial_magnitude( ...
            MagnitudeData_t, posPOSx_new_t, posPOSy_new_t, [], [], SizeImagePixels, pixelSize);
        All_RadialDist_Pos{end+1} = distMatPos; %
        All_RadialMag_Pos{end+1}  = magMatPos;  %
    end
    if ~isempty(posNEGx_new_t)
        [~,~, distMatNeg, magMatNeg] = compute_defect_radial_magnitude( ...
            MagnitudeData_t, [], [], posNEGx_new_t, posNEGy_new_t, SizeImagePixels, pixelSize);
        All_RadialDist_Neg{end+1} = distMatNeg; %
        All_RadialMag_Neg{end+1}  = magMatNeg;  %
    end

            end


% ======================== CALL analyze_defects2 (it loops internally) ========================
% Pass cell arrays of positions and fields; analyze_defects2 will loop over frames,
% concatenate results into single lists as you specified earlier.
[All_Magn, All_Div, All_u_ori, All_v_ori, All_u_vel, All_v_vel, ...
    All_xBox_vel, All_yBox_vel, All_xBox_ori, All_yBox_ori] = ...
    analyze_defects2(x_0, y_0, posPOSx, posPOSy, posNEGx, posNEGy, ...
    Mag_list, Div_list, scaleFactor, Ori_list, U_comp_list, V_comp_list, Phi_list_all, typeanal,fig_output,cellType);

% ======================== APPEND to per-main-directory accumulators ========================
% (these 'All_*_all' were initialized at the top of the mainDir loop)
All_Magn_all  = [All_Magn_all,  All_Magn];
All_Div_all   = [All_Div_all,   All_Div];
All_u_ori_all = [All_u_ori_all, All_u_ori];
All_v_ori_all = [All_v_ori_all, All_v_ori];
All_u_vel_all = [All_u_vel_all, All_u_vel];
All_v_vel_all = [All_v_vel_all, All_v_vel];
All_xBox_vel_all = [All_xBox_vel_all, All_xBox_vel];
All_yBox_vel_all = [All_yBox_vel_all, All_yBox_vel];
All_xBox_ori_all = [All_xBox_ori_all, All_xBox_ori];
All_yBox_ori_all = [All_yBox_ori_all, All_yBox_ori];
All_NemVelAngle_all = [All_NemVelAngle_all, NemVelAngle_list];




All_RadialDist_Pos_all = All_RadialDist_Pos; %# this is per-pair; adjust if needed
All_RadialMag_Pos_all  = All_RadialMag_Pos;
All_RadialDist_Neg_all = All_RadialDist_Neg;
All_RadialMag_Neg_all  = All_RadialMag_Neg;




if isempty(All_Magn), continue; end


        mainFolderName = Sub_maindir(sub).name;
       

        
end


if ~isempty(All_RadialDist_Pos)
        allDistPos = vertcat(All_RadialDist_Pos{:});
        allMagPos  = vertcat(All_RadialMag_Pos{:});
        allDistNeg = vertcat(All_RadialDist_Neg{:});
        allMagNeg  = vertcat(All_RadialMag_Neg{:});

       [avgDistPos, avgMagPos, stdMagPos, avgDistNeg, avgMagNeg, stdMagNeg] = average_radial_magnitudes(allDistPos, allMagPos, allDistNeg, allMagNeg);
    
    %     [Xq,Yq,nematicMat, velMat, magnMat, divMat] = average_fields(...
    % All_u_ori, All_v_ori, All_u_vel, All_v_vel, All_Magn, All_Div, ...
    % All_xBox_ori, All_yBox_ori, All_xBox_vel, All_yBox_vel);

    % [Xq,Yq,nematicMat, velMat, magnMat, divMat,stepoo] =average_fields_regular(...
    % All_u_ori, All_v_ori, All_u_vel, All_v_vel, All_Magn, All_Div, ...
    % All_xBox_ori, All_yBox_ori, All_xBox_vel, All_yBox_vel,1);

      [Xq,Yq,nematicMat, velMat, magnMat, divMat,stepoo,numdefanal] =average_fields_regular(...
    All_u_ori_all, All_v_ori_all, All_u_vel_all, All_v_vel_all, All_Magn_all, All_Div_all, ...
    All_xBox_ori_all, All_yBox_ori_all, All_xBox_vel_all, All_yBox_vel_all,1);

if typeanal=='pos'
% === Save averaged radial magnitudes ===

[ThetaValues, PhiValues, ThetaNeg, PhiNeg] = AngleDistributionv2_withCoords(nematicMat,Xq,Yq,0,0);

saveName = fullfile(output_dir, sprintf('%s_radialMagnitude_pos_avg.mat', mainFolderName));
save(saveName,'Xq','Yq','nematicMat','velMat','magnMat','divMat' ,'avgDistPos', 'avgMagPos', 'avgDistNeg', 'avgMagNeg','stdMagPos','stdMagNeg','ThetaValues','PhiValues','numdefanal');
fprintf('✅ Averaged radial magnitudes saved: %s\n', saveName);
end

if typeanal=='neg'
% === Save averaged radial magnitudes ===
saveName = fullfile(output_dir, sprintf('%s_radialMagnitude_neg_avg.mat', mainFolderName));
save(saveName,'Xq','Yq','nematicMat','velMat','magnMat','divMat' ,'avgDistPos', 'avgMagPos', 'avgDistNeg', 'avgMagNeg','stdMagPos','stdMagNeg','numdefanal');
fprintf('✅ Averaged radial magnitudes saved: %s\n', saveName);
end



darkGray = [0.3 0.3 0.3];
end

%Plotting the magnitude around defects, velocity vectors and divergence
    %around 

    if typeanal=='pos'
    
    % === Plot velocity field colored by magnitude ===
output_name_vel = sprintf('%s_pos_nematicandvelocity.svg', mainFolderName);

U=cos(velMat);
V=sin(velMat);

nColors = 256;
cmap = parula(nColors);

cMin = min(magnMat(:)); % overall min
cMax = max(magnMat(:)); % overall max

colorIdx = round(1 + (nColors - 1) * (magnMat - cMin) / (cMax - cMin));
colorIdx = min(max(colorIdx, 1), nColors); % clamp element-wise




% --- Scaling factor setup ---
userScale = 0.8;   % adjust this manually (like AutoScaleFactor)
scaleFactor = userScale / (cMax - cMin);  % inversely proportional to range

fig1=figure; hold on; axis equal;
% for k = 1:size(Xq,1)
%     for j = 1:size(Xq,1)
%         rgb = cmap(colorIdx(k,j), :);
% 
%         % scale arrow by normalized magnitude
%         magNorm = (magnMat(k,j) - cMin);
%         quiver(Xq(k,j), Yq(k,j), magNorm * U(k,j) * scaleFactor, ...
%                magNorm * V(k,j) * scaleFactor, ...
%                "ShowArrowHead", "on", "Color", rgb,'LineWidth',1.5);
%     end
% end
% --- Dimensions ---
[nRows, nCols] = size(nematicMat);

% Centre de la matrice
midR = floor(nRows/2)+1;
midC = floor(nCols/2)+1;

% Indices du bloc 4x4 centré
rows = midR-4 : midR+4;
cols = midC-4 : midC+4;

% Extraction des blocs
Xq4 = Xq(rows, cols);
Yq4 = Yq(rows, cols);

magn4 = magnMat(rows, cols);       % magnitude vitesse
vel4  = velMat(rows, cols);        % orientation vitesse (angles)
nem4  = nematicMat(rows, cols);    % angle nématique (directeur)

% --- QUIVER ---
hold on;

% VECTEURS ROUGES : champ vitesse pondéré
quiver(Xq4, Yq4, magn4.*cos(vel4), magn4.*sin(vel4), 0.8, ...
       'Color','red','LineWidth',1.2);
lightGray = [0.7 0.7 0.7];
% VECTEURS NOIRS : champ nématique
quiver( Xq4, Yq4, cos(nem4), sin(nem4), ...
        1, 'ShowArrowHead','off', 'Color',lightGray, 'LineWidth',1.5 );

set(gca,'XTickLabelRotation',0)
set(gca,'linew',2,'fontsize',20,'box','off','InnerPosition',[0.3 0.3 0.5 0.5],'Units','centimeters');
set(gcf,'OuterPosition',[300 100 650 550])
pbaspect([1 1 1]);
box off
grid off
axis off




saveas(fig1, fullfile(output_dir,output_name_vel),'svg'); 


% === Plot divergence map with nematic overlay ===
output_name_nem = sprintf('%s_%s_defect_pos_divcolormap.svg', mainFolderName);

fig2=figure; hold on;

% --- Background divergence map
pcolor(Xq, Yq, divMat);
shading interp; 
set(gca, 'YDir', 'normal');
colormap(parula);
cb = colorbar;
cb.Label.String = 'Divergence';
axis equal tight;

quiver(Xq,Yq,cos(nematicMat),sin(nematicMat),1,"ShowArrowHead","off",'color',darkGray,'LineWidth',1.5);


saveas(fig2, fullfile(output_dir,output_name_nem),'svg'); 


    end

if typeanal=='neg'

    % === Plot velocity field colored by magnitude ===
output_name_vel = sprintf('%s_neg_nematicandvelocity.svg', mainFolderName);

U=cos(velMat);
V=sin(velMat);

nColors = 256;
cmap = parula(nColors);

cMin = min(magnMat(:)); % overall min
cMax = max(magnMat(:)); % overall max

colorIdx = round(1 + (nColors - 1) * (magnMat - cMin) / (cMax - cMin));
colorIdx = min(max(colorIdx, 1), nColors); % clamp element-wise




% --- Scaling factor setup ---
userScale = 0.8;   % adjust this manually (like AutoScaleFactor)
scaleFactor = userScale / (cMax - cMin);  % inversely proportional to range

fig1=figure; hold on; axis equal;
% for k = 1:size(Xq,1)
%     for j = 1:size(Xq,1)
%         rgb = cmap(colorIdx(k,j), :);
% 
%         % scale arrow by normalized magnitude
%         magNorm = (magnMat(k,j) - cMin);
%         quiver(Xq(k,j), Yq(k,j), magNorm * U(k,j) * scaleFactor, ...
%                magNorm * V(k,j) * scaleFactor, ...
%                "ShowArrowHead", "on", "Color", rgb,'LineWidth',1.5);
%     end
% end
[nRows, nCols] = size(nematicMat);

% Centre de la matrice
midR = floor(nRows/2)+1;
midC = floor(nCols/2)+1;

% Indices du bloc 4x4 centré
rows = midR-4 : midR+4;
cols = midC-4 : midC+4;

% Extraction des blocs
Xq4 = Xq(rows, cols);
Yq4 = Yq(rows, cols);

magn4 = magnMat(rows, cols);       % magnitude vitesse
vel4  = velMat(rows, cols);        % orientation vitesse (angles)
nem4  = nematicMat(rows, cols);    % angle nématique (directeur)

% --- QUIVER ---
hold on;

% VECTEURS BLEUS : champ vitesse pondéré
quiver(Xq4, Yq4, magn4.*cos(vel4), magn4.*sin(vel4), 0.8, ...
       'Color','blue','LineWidth',1.2);
lightGray = [0.7 0.7 0.7];

% VECTEURS NOIRS : champ nématique
quiver( Xq4, Yq4, cos(nem4), sin(nem4), ...
        1, 'ShowArrowHead','off', 'Color',lightGray, 'LineWidth',1.5 );

set(gca,'XTickLabelRotation',0)
set(gca,'linew',2,'fontsize',20,'box','off','InnerPosition',[0.3 0.3 0.5 0.5],'Units','centimeters');
set(gcf,'OuterPosition',[300 100 650 550])
pbaspect([1 1 1]);
box off
grid off
axis off

saveas(fig1, fullfile(output_dir,output_name_vel),'svg'); 


% === Plot divergence map with nematic overlay ===
output_name_nem = sprintf('%s_defect_neg_divcolormap.svg', mainFolderName);

fig2=figure; hold on;

% --- Background divergence map
pcolor(Xq, Yq, divMat);
shading interp; 
set(gca, 'YDir', 'normal');
colormap(parula);
cb = colorbar;
cb.Label.String = 'Divergence';
axis equal tight;

quiver(Xq,Yq,cos(nematicMat),sin(nematicMat),1,"ShowArrowHead","off",'color',darkGray,'LineWidth',1.5);


saveas(fig2, fullfile(output_dir,output_name_nem),'svg'); 

end
end
end
% Functions 

function [All_Magn,All_Div, All_u_ori, All_v_ori,All_u_vel,All_v_vel,All_xBox_vel,All_yBox_vel,All_xBox_ori,All_yBox_ori] = ...
    analyze_defects2(xMat,yMat,posPOSx, posPOSy, posNEGx, posNEGy, MagnitudeData,Divdata, Scale_factor,OrientationMat,UMat,VMat,Philist,type_analysis,fig_output,cellType)

    % ----------------------------------------------
    %  ❗ Détection du multi-frame 
    % ----------------------------------------------
    if iscell(OrientationMat)
        numFrames = length(OrientationMat);
    else
        numFrames = 1;
    end

    size_patch=6;

    % ----------------------------------------------
    %  ❗ Initialisation des outputs (une seule liste)
    % ----------------------------------------------
    All_Magn  = {};
    All_Div   = {};
    All_u_ori = {};
    All_v_ori = {};
    All_u_vel = {};
    All_v_vel = {};
    All_xBox_vel = {};
    All_yBox_vel = {};
    All_xBox_ori = {};
    All_yBox_ori = {};

    % ----------------------------------------------
    %  ❗ Boucle temporelle interne
    % ----------------------------------------------
    for t = 1:numFrames

        % ------- Extraire la frame t si cell array -------
        if iscell(MagnitudeData), Magn_t = MagnitudeData{t}; else, Magn_t = MagnitudeData; end
        if iscell(Divdata),       Div_t = Divdata{t};       else, Div_t = Divdata;       end
        if iscell(OrientationMat),Ori_t = OrientationMat{t};else, Ori_t = OrientationMat;end
        if iscell(UMat),          U_t = UMat{t};            else, U_t = UMat;            end
        if iscell(VMat),          V_t = VMat{t};            else, V_t = VMat;            end
        if iscell(Philist),       Phi_t = Philist{t};       else, Phi_t = Philist;       end

        % ------- Positions des défauts -------
        if iscell(posPOSx), posPOSx_t = posPOSx{t}; else, posPOSx_t = posPOSx; end
        if iscell(posPOSy), posPOSy_t = posPOSy{t}; else, posPOSy_t = posPOSy; end
        if iscell(posNEGx), posNEGx_t = posNEGx{t}; else, posNEGx_t = posNEGx; end
        if iscell(posNEGy), posNEGy_t = posNEGy{t}; else, posNEGy_t = posNEGy; end

        % ------- Défaillances possibles -------
        posPOSx_t = posPOSx_t(:)';
        posPOSy_t = posPOSy_t(:)';
        posNEGx_t = posNEGx_t(:)';
        posNEGy_t = posNEGy_t(:)';

        % Process POS or NEG depending on type
        if strcmp(type_analysis,"pos")
            posX_list = posPOSx_t;
            posY_list = posPOSy_t;
        elseif strcmp(type_analysis,"neg")
            posX_list = posNEGx_t;
            posY_list = posNEGy_t;
        else
            error('type_analysis must be "pos" or "neg"');
        end

        % ----------------------------------------------
        %  ❗ Boucle défauts dans cette frame
        % ----------------------------------------------
        num_defects = length(posX_list);
        pos_new_x = round(posX_list * Scale_factor);
        pos_new_y = round(posY_list * Scale_factor);

        for i = 1:num_defects

            x_ori = posX_list(i);
            y_ori = posY_list(i);
            x = pos_new_x(i);
            y = pos_new_y(i);

            % Angle PHI
            Angle_def = Phi_t(x_ori, y_ori);

            % Boîtes (orientation : coordonnées originales)
            OrientBox = extractSquareBox(Ori_t', x_ori, y_ori, size_patch);

           

            % Boîtes magnitude/div/velocity (coords rescalées)
            MagnBox = extractSquareBox(Magn_t, x, y, size_patch);
            DivBox  = extractSquareBox(Div_t,  x, y, size_patch);
            Ubox    = extractSquareBox(U_t,    x, y, size_patch);
            Vbox    = extractSquareBox(V_t,    x, y, size_patch);

            % Rotation orientation et vitesse
            [Urot, Vrot, xrot, yrot] = rotateMatrix(cos(OrientBox), sin(OrientBox), pi - Angle_def);
            [U_rot_interp, V_rot_interp, xu, yu] = rotateMatrix(Ubox, Vbox, pi - Angle_def);

            
            [ThetaValues, PhiValues, ThetaNeg, PhiNeg] = AngleDistributionv2_withCoords(atan2(V_rot_interp,U_rot_interp),xrot,yrot,0,0);


           

            % ----------------------------------------------
            %  ❗ Append aux outputs globaux
            % ----------------------------------------------
            All_Magn{end+1}  = MagnBox;
            All_Div{end+1}   = DivBox;

            All_u_ori{end+1} = Urot;
            All_v_ori{end+1} = Vrot;

            All_u_vel{end+1} = U_rot_interp;
            All_v_vel{end+1} = V_rot_interp;

            All_xBox_vel{end+1} = xu;
            All_yBox_vel{end+1} = yu;

            All_xBox_ori{end+1} = xrot;
            All_yBox_ori{end+1} = yrot;

           
        end % boucle défauts

    end % boucle frames

end


function [Xq, Yq, nematicMat, velMat, magnMat, divMat, gridStep,count_ori] = average_fields_regular( ...
    All_u_ori, All_v_ori, All_u_vel, All_v_vel, ...
    All_Magn, All_Div, ...
    All_xBox_ori, All_yBox_ori, All_xBox_vel, All_yBox_vel, ...
    gridStep)
% AVERAGE_FIELDS_REGULAR Average nematic orientation, velocity, magnitude, and divergence onto a regular grid.
%
% Note: nematic orientation is handled in double-angle space (cos(2θ), sin(2θ)).
%       Per your request, the double-angle components are SUMMED across datasets
%       (no division by count for the nematic components). Cells with zero
%       contributions are set to NaN.
%
% Inputs:
%   All_u_ori, All_v_ori   - cell arrays of orientation vector components (u,v)
%   All_u_vel, All_v_vel   - cell arrays of velocity vector components (u,v)
%   All_Magn, All_Div      - cell arrays of scalar fields (magnitude, divergence)
%   All_xBox_ori, All_yBox_ori - cell arrays of coordinates for orientation vectors
%   All_xBox_vel, All_yBox_vel - cell arrays of coordinates for velocity/magn/div
%   gridStep (optional)    - scalar grid step; if empty, provide a value externally
%
% Outputs:
%   Xq, Yq     - meshgrid for the regular grid
%   nematicMat - averaged nematic orientation (radians) computed from summed double-angle components
%   velMat     - averaged velocity orientation (radians) (usual mean)
%   magnMat    - averaged magnitude
%   divMat     - averaged divergence
%   gridStep   - step used

    % --- 1. Collect all coordinates to determine bounds
    allX = []; allY = [];
    for k = 1:numel(All_xBox_ori)
        allX = [allX; All_xBox_ori{k}(:); All_xBox_vel{k}(:)];
        allY = [allY; All_yBox_ori{k}(:); All_yBox_vel{k}(:)];
    end
    xMin = min(allX(:)); xMax = max(allX(:));
    yMin = min(allY(:)); yMax = max(allY(:));

    % If gridStep is empty or <= 0, pick a default (e.g., 30)
    if nargin < 11 || isempty(gridStep) || gridStep <= 0
        gridStep = 30;
    end

    % --- 3. Build regular grid
    xq = xMin:gridStep:xMax;
    yq = yMin:gridStep:yMax;
    [Xq, Yq] = meshgrid(xq, yq);
    sz = size(Xq);

    % --- 4. Initialize accumulators
    % For nematic orientation (we will SUM double-angle components across datasets)
    sumU_ori = zeros(sz); 
    sumV_ori = zeros(sz); 
    count_ori = zeros(sz);

    % For velocity (average normally)
    sumU_vel = zeros(sz); 
    sumV_vel = zeros(sz); 
    count_vel = zeros(sz);

    % For magnitude and divergence (average normally)
    sumMagn = zeros(sz); countMagn = zeros(sz);
    sumDiv  = zeros(sz); countDiv  = zeros(sz);

    % --- 5. Loop over all datasets
    nDatasets = numel(All_u_ori);
    for k = 1:nDatasets
        % fetch fields for dataset k
        uOri = All_u_ori{k};
        vOri = All_v_ori{k};
        uVel = All_u_vel{k};
        vVel = All_v_vel{k};
        magn = All_Magn{k};
        divv = All_Div{k};

        Xo = All_xBox_ori{k}; Yo = All_yBox_ori{k};
        Xv = All_xBox_vel{k}; Yv = All_yBox_vel{k};

        % ---- NEMATIC: compute orientation angle then double-angle components
        Ori = atan2(vOri, uOri);          % director (radians), head-tail ambiguous
        U2 = cos(2* Ori);              % double-angle X component
        V2 = sin(2 * Ori);              % double-angle Y component

        % Interpolate double-angle components onto regular grid
        U2_i = griddata(Xo, Yo, U2, Xq, Yq, 'natural');
        V2_i = griddata(Xo, Yo, V2, Xq, Yq, 'natural');

        % Accumulate SUM of double-angle components (no normalization here)
        validOri = ~isnan(U2_i) & ~isnan(V2_i);
        sumU_ori(validOri) = sumU_ori(validOri) + U2_i(validOri);
        sumV_ori(validOri) = sumV_ori(validOri) + V2_i(validOri);
        count_ori(validOri) = count_ori(validOri) + 1;  % keep count for masking only

        % ---- Velocity: interpolate and accumulate (normal averaging)
        uVel_i = griddata(Xv, Yv, uVel, Xq, Yq, 'natural');
        vVel_i = griddata(Xv, Yv, vVel, Xq, Yq, 'natural');
        validVel = ~isnan(uVel_i) & ~isnan(vVel_i);
        sumU_vel(validVel) = sumU_vel(validVel) + uVel_i(validVel);
        sumV_vel(validVel) = sumV_vel(validVel) + vVel_i(validVel);
        count_vel(validVel) = count_vel(validVel) + 1;

        % ---- Magnitude: interpolate and accumulate
        magn_i = griddata(Xv, Yv, magn, Xq, Yq, 'natural');
        validMagn = ~isnan(magn_i);
        sumMagn(validMagn) = sumMagn(validMagn) + magn_i(validMagn);
        countMagn(validMagn) = countMagn(validMagn) + 1;

        % ---- Divergence: interpolate and accumulate
        div_i = griddata(Xv, Yv, divv, Xq, Yq, 'natural');
        validDiv = ~isnan(div_i);
        sumDiv(validDiv) = sumDiv(validDiv) + div_i(validDiv);
        countDiv(validDiv) = countDiv(validDiv) + 1;
    end

    % --- 6. Compute averaged matrices (outputs)
    nematicMat = NaN(sz);
    velMat     = NaN(sz);
    magnMat    = NaN(sz);
    divMat     = NaN(sz);

    % NEMATIC: reconstruct angle from SUMMED double-angle components (no division)
    % Set summed components to NaN where there was no contribution (count_ori==0)
    U_mean = sumU_ori;
    V_mean = sumV_ori;
    U_mean(count_ori == 0) = NaN;
    V_mean(count_ori == 0) = NaN;

    % director angle from summed double-angle components:
    % theta = 0.5 * atan2(sum_sin2, sum_cos2)
    validOriMask = (count_ori > 0);
    nematicMat(validOriMask) = 0.5 * atan2(V_mean(validOriMask), U_mean(validOriMask));

    % VELOCITY: normal averaging (sum / count)
    validVelMask = (count_vel > 0);
    Uv_mean = NaN(sz); Vv_mean = NaN(sz);
    Uv_mean(validVelMask) = sumU_vel(validVelMask);
    Vv_mean(validVelMask) = sumV_vel(validVelMask) ;
    velMat(validVelMask) = atan2(Vv_mean(validVelMask), Uv_mean(validVelMask));

    % MAGNITUDE: average
    validMagnMask = (countMagn > 0);
    magnMat(validMagnMask) = sumMagn(validMagnMask) ./ countMagn(validMagnMask);

    % DIVERGENCE: average
    validDivMask = (countDiv > 0);
    divMat(validDivMask) = sumDiv(validDivMask) ./ countDiv(validDivMask);

    fprintf('Averaged %d datasets onto regular grid (%dx%d) with gridStep=%g\n', ...
        nDatasets, size(Xq,1), size(Xq,2), gridStep);
end

% helper to build colormap
function cmap = createColorMap(families, baseColors)
    famNames = fieldnames(families);
    cmap = containers.Map();
    for f = 1:numel(famNames)
        fam = famNames{f};
        members = families.(fam);
        base = baseColors.(fam);
        N = numel(members);
        tVals = linspace(0.3,0.7,N);
        for k = 1:N
            shade = base*(1-tVals(k)) + tVals(k)*[1 1 1];
            cmap(lower(members{k})) = max(0,min(1,shade));
        end
    end
end