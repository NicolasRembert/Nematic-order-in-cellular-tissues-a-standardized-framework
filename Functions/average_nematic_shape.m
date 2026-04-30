function [Xq,Yq,nematicMat,gridStep,count_ori] = ...
average_nematic_shape(All_u_ori,All_v_ori,All_xBox_ori,All_yBox_ori,gridStep)
%AVERAGE_NEMATIC_SHAPE Average nematic orientation fields on a common grid
%
% [Xq,Yq,nematicMat,gridStep,count_ori] = average_nematic_shape( ...
%     All_u_ori, All_v_ori, All_xBox_ori, All_yBox_ori, gridStep)
%
% Inputs
%   All_u_ori      - cell array containing u components of orientation fields
%   All_v_ori      - cell array containing v components of orientation fields
%   All_xBox_ori   - cell array of x coordinates for each orientation patch
%   All_yBox_ori   - cell array of y coordinates for each orientation patch
%   gridStep       - scalar defining the spacing of the interpolation grid
%
% Outputs
%   Xq, Yq         - meshgrid coordinates of the common interpolation grid
%   nematicMat     - averaged nematic orientation field (angle in radians)
%   gridStep       - grid spacing used for interpolation
%   count_ori      - number of datasets contributing to each grid point
%
% The function interpolates multiple orientation patches onto a common grid
% and computes the average nematic orientation field using the 2θ
% representation.


allX = [];
allY = [];

for k = 1:numel(All_xBox_ori)
    allX = [allX; All_xBox_ori{k}(:)];
    allY = [allY; All_yBox_ori{k}(:)];
end

xMin = min(allX);
xMax = max(allX);
yMin = min(allY);
yMax = max(allY);


xq = xMin:gridStep:xMax;
yq = yMin:gridStep:yMax;

[Xq,Yq] = meshgrid(xq,yq);

sz = size(Xq);

sumU = zeros(sz);
sumV = zeros(sz);
count_ori = zeros(sz);

nDatasets = numel(All_u_ori);

for k = 1:nDatasets

    uOri = All_u_ori{k};
    vOri = All_v_ori{k};

    Xo = All_xBox_ori{k};
    Yo = All_yBox_ori{k};

    Ori = atan2(vOri,uOri);

    U2 = cos(2*Ori);
    V2 = sin(2*Ori);

    U2_i = griddata(Xo,Yo,U2,Xq,Yq,'natural');
    V2_i = griddata(Xo,Yo,V2,Xq,Yq,'natural');

    valid = ~isnan(U2_i);

    sumU(valid) = sumU(valid) + U2_i(valid);
    sumV(valid) = sumV(valid) + V2_i(valid);

    count_ori(valid) = count_ori(valid) + 1;

end

nematicMat = NaN(sz);

valid = count_ori>0;

nematicMat(valid) = 0.5*atan2(sumV(valid),sumU(valid));

fprintf('Averaged %d nematic datasets onto grid (%dx%d)\n',...
nDatasets,size(Xq,1),size(Xq,2))

end