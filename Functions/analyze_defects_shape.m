function [All_u_ori, All_v_ori, All_xBox_ori, All_yBox_ori] = ...
analyze_defects_shape(posPOSx,posPOSy,posNEGx,posNEGy,OrientationMat,Philist,type_analysis)
%ANALYZE_DEFECTS_SHAPE Extract and align orientation patches around defects
%
% [All_u_ori, All_v_ori, All_xBox_ori, All_yBox_ori] = analyze_defects_shape( ...
%     posPOSx, posPOSy, posNEGx, posNEGy, OrientationMat, Philist, type_analysis)
%
% Inputs
%   posPOSx, posPOSy   - x and y coordinates of positive defects (+1/2)
%   posNEGx, posNEGy   - x and y coordinates of negative defects (-1/2)
%   OrientationMat     - orientation field (matrix or cell array for multi-frame data)
%   Philist            - local defect orientation field (same format as OrientationMat)
%   type_analysis      - defect type to analyze: "pos" or "neg"
%
% Outputs
%   All_u_ori          - cell array containing rotated u components of the orientation field
%   All_v_ori          - cell array containing rotated v components of the orientation field
%   All_xBox_ori       - cell array of rotated x coordinates of the extracted patches
%   All_yBox_ori       - cell array of rotated y coordinates of the extracted patches
%
% The function extracts orientation patches around defects and rotates them
% into the local defect reference frame.

% Detection multi-frame

if iscell(OrientationMat)
    numFrames = length(OrientationMat);
else
    numFrames = 1;
end

size_patch = 6;

All_u_ori = {};
All_v_ori = {};
All_xBox_ori = {};
All_yBox_ori = {};

for t = 1:numFrames

    % Extraction frame
    if iscell(OrientationMat), Ori_t = OrientationMat{t}; else, Ori_t = OrientationMat; end
    if iscell(Philist), Phi_t = Philist{t}; else, Phi_t = Philist; end

    % Positions défauts
    if iscell(posPOSx), posPOSx_t = posPOSx{t}; else, posPOSx_t = posPOSx; end
    if iscell(posPOSy), posPOSy_t = posPOSy{t}; else, posPOSy_t = posPOSy; end
    if iscell(posNEGx), posNEGx_t = posNEGx{t}; else, posNEGx_t = posNEGx; end
    if iscell(posNEGy), posNEGy_t = posNEGy{t}; else, posNEGy_t = posNEGy; end

    posPOSx_t = posPOSx_t(:)';
    posPOSy_t = posPOSy_t(:)';
    posNEGx_t = posNEGx_t(:)';
    posNEGy_t = posNEGy_t(:)';

    % Choix type défaut
    if strcmp(type_analysis,"pos")
        posX_list = posPOSx_t;
        posY_list = posPOSy_t;
    elseif strcmp(type_analysis,"neg")
        posX_list = posNEGx_t;
        posY_list = posNEGy_t;
    else
        error('type_analysis must be pos or neg')
    end

    num_defects = length(posX_list);

    for i = 1:num_defects

        x_ori = posX_list(i);
        y_ori = posY_list(i);

        Angle_def = Phi_t(x_ori,y_ori);

        % Extraction patch orientation
        OrientBox = extractSquareBox(Ori_t',x_ori,y_ori,size_patch);

        % Rotation dans le repère du défaut
        [Urot,Vrot,xrot,yrot] = rotateMatrix(...
            cos(OrientBox),...
            sin(OrientBox),...
            pi - Angle_def);

        % Sauvegarde
        All_u_ori{end+1} = Urot;
        All_v_ori{end+1} = Vrot;

        All_xBox_ori{end+1} = xrot;
        All_yBox_ori{end+1} = yrot;

    end

end

end