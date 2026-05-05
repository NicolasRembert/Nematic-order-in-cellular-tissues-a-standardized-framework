function DefectStructure = CalculateMSD_ParPerp(DefectStructure, pixelSize, nodeSpacing)
% CalculateMSD_ParPerp
% Calcule la MSD parallèle et perpendiculaire à l’axe du défaut
% uniquement pour les défauts positifs (+1/2)

for def = 1:numel(DefectStructure)

    % ---------- uniquement défauts positifs ----------
    if DefectStructure(def).charge <= 0
        continue
    end

    positions = DefectStructure(def).positions;
    phi       = DefectStructure(def).phi(:);

    % ---------- se caler sur la longueur minimale ----------
    nPos = size(positions,1);
    nPhi = numel(phi);
    n    = min(nPos, nPhi);

    positions = positions(1:n,:);
    phi       = phi(1:n);

    % ---------- filtrage ----------
    valid = positions(:,1) ~= 0 & positions(:,2) ~= 0 & ~isnan(phi);

    positions = positions(valid,:);
    phi       = phi(valid);

    if size(positions,1) < 2
        DefectStructure(def).MSD_par  = [];
        DefectStructure(def).MSD_perp = [];
        continue
    end

    % ---------- paramètres ----------
    lifetime = size(positions,1);
    numberOfDeltaT = floor(lifetime/2);

    MSD_par  = zeros(numberOfDeltaT,1);
    MSD_perp = zeros(numberOfDeltaT,1);

    scale = pixelSize * nodeSpacing;

    % ---------- calcul MSD ----------
    for dt = 1:numberOfDeltaT

        dx = positions(1+dt:end,1) - positions(1:end-dt,1);
        dy = positions(1+dt:end,2) - positions(1:end-dt,2);

        phi0 = phi(1:end-dt);

        % vecteurs unitaires
        n_par  = [cos(phi0), sin(phi0)];
        n_perp = [-sin(phi0), cos(phi0)];

        % projections
        dr_par  = dx .* n_par(:,1)  + dy .* n_par(:,2);
        dr_perp = dx .* n_perp(:,1) + dy .* n_perp(:,2);

        MSD_par(dt)  = mean((dr_par  * scale).^2);
        MSD_perp(dt) = mean((dr_perp * scale).^2);
    end

    % ---------- stockage ----------
    DefectStructure(def).MSD_par  = MSD_par;
    DefectStructure(def).MSD_perp = MSD_perp;

end
end
