function pairs = pairOppositeDefects_phiproduct(defects, tol)

if nargin < 2
    tol = 1e-6;
end

posIdxAll = find([defects.charge] > 0);
negIdxAll = find([defects.charge] < 0);

% nombre max de pas de temps
nT = 0;
for d = 1:numel(defects)
    nT = max(nT, size(defects(d).positions,1));
end

pairs = struct('time', cell(nT,1), 'list', cell(nT,1));

for t = 1:nT

    list = struct('idxPos', {}, 'idxNeg', {}, ...
                  'distance', {}, 'rpn', {});

    for iPos = 1:numel(posIdxAll)

        idxPos = posIdxAll(iPos);

        % ---- vérifications robustes ----
        if size(defects(idxPos).positions,1) < t
            continue
        end

        pos = defects(idxPos).positions(t,:);
        if all(pos == 0)
            continue
        end

        phi = defects(idxPos).phi;
        if isempty(phi) || t > numel(phi) || isnan(phi(t))
            continue
        end

        % ---- recherche du - le plus proche ----
        rmin = Inf;
        idxNegBest = NaN;

        for iNeg = 1:numel(negIdxAll)

            idxNeg = negIdxAll(iNeg);

            if size(defects(idxNeg).positions,1) < t
                continue
            end

            neg = defects(idxNeg).positions(t,:);
            if all(neg == 0)
                continue
            end

            d = norm(neg - pos);

            if d < rmin
                rmin = d;
                idxNegBest = idxNeg;
            end
        end

        if isnan(idxNegBest) || rmin < tol
            continue
        end

        % ---- calcul r_hat ----
        neg = defects(idxNegBest).positions(t,:);
        r = neg - pos;
        r_hat = r / norm(r);

        % ---- polarité du + ----
        phi_t = phi(t);
        p_vec = [cos(phi_t), sin(phi_t)];

        rpn = dot(r_hat, p_vec);

        % ---- stockage ----
        list(end+1).idxPos   = idxPos; %#ok<AGROW>
        list(end).idxNeg     = idxNegBest;
        list(end).distance   = rmin * 30 * 0.65;
        list(end).rpn        = rpn;

    end

    pairs(t).time = t;
    pairs(t).list = list;

end
end


% function pairs = pairOppositeDefects_phiproduct(defects, tol)
% % Pair chaque défaut +1/2 avec le -1/2 le plus proche
% % Calcul du produit scalaire entre l’axe phi du +1/2
% % et le vecteur normalisé + -> -
% % Le calcul est fait automatiquement à tous les pas de temps
% 
% if nargin < 2
%     tol = 1e-6;
% end
% 
% %% Indices globaux
% posIdxAll = find([defects.charge] > 0);
% negIdxAll = find([defects.charge] < 0);
% 
% %% Nombre max de pas de temps
% nT = 0;
% for d = 1:numel(defects)
%     nT = max(nT, size(defects(d).positions,1));
% end
% 
% pairs = struct([]);
% 
% %% ===================== boucle temps =====================
% for t = 1:nT
% 
%     list = struct( ...
%         'idxPos', {}, ...
%         'idxNeg', {}, ...
%         'distance', {}, ...
%         'dotPhi', {}, ...
%         'r_hat', {} );
% 
%     %% ---------- défauts positifs valides ----------
%     posIdx = [];
%     for k = 1:numel(posIdxAll)
%         idx = posIdxAll(k);
% 
%         if size(defects(idx).positions,1) < t
%             continue
%         end
%         p = defects(idx).positions(t,:);
%         if p(1)==0 && p(2)==0
%             continue
%         end
% 
%         if numel(defects(idx).phi) < t
%             continue
%         end
%         if isnan(defects(idx).phi(t))
%             continue
%         end
% 
%         posIdx(end+1) = idx; %#ok<AGROW>
%     end
% 
%     %% ---------- défauts négatifs valides ----------
%     negIdx = [];
%     for k = 1:numel(negIdxAll)
%         idx = negIdxAll(k);
% 
%         if size(defects(idx).positions,1) < t
%             continue
%         end
%         p = defects(idx).positions(t,:);
%         if p(1)==0 && p(2)==0
%             continue
%         end
% 
%         negIdx(end+1) = idx; %#ok<AGROW>
%     end
% 
%     if isempty(posIdx) || isempty(negIdx)
%         pairs(t).time = t;
%         pairs(t).list = [];
%         continue
%     end
% 
%     %% Positions
%     Ppos = zeros(numel(posIdx),2);
%     for i = 1:numel(posIdx)
%         Ppos(i,:) = defects(posIdx(i)).positions(t,:);
%     end
% 
%     Pneg = zeros(numel(negIdx),2);
%     for j = 1:numel(negIdx)
%         Pneg(j,:) = defects(negIdx(j)).positions(t,:);
%     end
% 
%     %% ---------- association + -> - ----------
%     for i = 1:numel(posIdx)
% 
%         rmin = Inf;
%         jmin = NaN;
% 
%         for j = 1:numel(negIdx)
%             r = Pneg(j,:) - Ppos(i,:);
%             d = norm(r);
%             if d < rmin
%                 rmin = d;
%                 jmin = j;
%             end
%         end
% 
%         if isnan(jmin) || rmin < tol
%             continue
%         end
% 
%         % vecteur normalisé
%         r = Pneg(jmin,:) - Ppos(i,:);
%         r_hat = r / norm(r);
% 
%         % axe phi du défaut positif
%         phi = defects(posIdx(i)).phi(t);
%         n   = [cos(phi), sin(phi)];
% 
%         dotPhi = dot(n, r_hat);
% 
%         % stockage
%         list(end+1).idxPos   = posIdx(i); %#ok<AGROW>
%         list(end).idxNeg     = negIdx(jmin);
%         list(end).distance   = rmin * 30 * 0.65;
%         list(end).dotPhi     = dotPhi;
%         list(end).r_hat      = r_hat;
%     end
% 
%     %% ---------- sauvegarde par temps ----------
%     pairs(t).time = t;
%     pairs(t).list = list;
% 
% end
% end
