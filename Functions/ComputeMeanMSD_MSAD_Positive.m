function [ ...
    meanMSD_par, ...
    meanMSD_perp, ...
    meanMSAD, ...
    coeffs] = ComputeMeanMSD_MSAD_Positive(defects, dt)
% Moyennes MSD parallèle / perpendiculaire et MSAD
% pour les défauts positifs uniquement (+1/2)
%
% INPUT
%   defects : structure array
%   dt      : pas de temps (ex : 15 min)
%
% OUTPUT
%   meanMSD_par
%   meanMSD_perp
%   meanMSAD
%   coeffs.MSD_par
%   coeffs.MSD_perp
%   coeffs.MSAD

if nargin < 2
    dt = 1;
end

%% ---------- défauts positifs ----------
idxPos = find([defects.charge] > 0);
defects = defects(idxPos);
nDef = numel(defects);

%% ---------- longueurs ----------
len_par  = arrayfun(@(d) numel(d.MSD_par),  defects);
len_perp = arrayfun(@(d) numel(d.MSD_perp), defects);
len_ang  = arrayfun(@(d) numel(d.AngularMSD),     defects);

maxLen = max([len_par, len_perp, len_ang]);


MSDpar_all  = nan(nDef, maxLen);
MSDperp_all = nan(nDef, maxLen);
MSAD_all    = nan(nDef, maxLen);

for i = 1:nDef
    MSDpar_all(i, 1:len_par(i))   = defects(i).MSD_par(:)';
    MSDperp_all(i,1:len_perp(i))  = defects(i).MSD_perp(:)';
    MSAD_all(i,  1:len_ang(i))    = defects(i).AngularMSD(:)';
end

%% ---------- moyennes ----------
meanMSD_par  = nanmean(MSDpar_all,  1);
meanMSD_perp = nanmean(MSDperp_all, 1);
meanMSAD     = nanmean(MSAD_all,    1);

%% ---------- coefficients diffusifs ----------
t = (1:maxLen) * dt;

valid = ~isnan(meanMSD_par);
coeffs.MSD_par = polyfit(t(valid), meanMSD_par(valid), 1);

valid = ~isnan(meanMSD_perp);
coeffs.MSD_perp = polyfit(t(valid), meanMSD_perp(valid), 1);

valid = ~isnan(meanMSAD);
coeffs.MSAD = polyfit(t(valid), meanMSAD(valid), 1);

end
