function average_radial_profiles(dist_pos_all, mag_pos_all, dist_neg_all, mag_neg_all, output_dir)
% Averages magnitude vs distance across all defects

binSize = 2; % µm bin size
maxDist = max([dist_pos_all, dist_neg_all]);
bins = 0:binSize:maxDist;

mean_pos = zeros(size(bins));
std_pos  = zeros(size(bins));
mean_neg = zeros(size(bins));
std_neg  = zeros(size(bins));

for b = 1:length(bins)-1
    idx_p = dist_pos_all >= bins(b) & dist_pos_all < bins(b+1);
    idx_n = dist_neg_all >= bins(b) & dist_neg_all < bins(b+1);

    mean_pos(b) = mean(mag_pos_all(idx_p), 'omitnan');
    std_pos(b)  = std(mag_pos_all(idx_p),  'omitnan');
    mean_neg(b) = mean(mag_neg_all(idx_n), 'omitnan');
    std_neg(b)  = std(mag_neg_all(idx_n),  'omitnan');
end

save(fullfile(output_dir, 'Defect_Magnitude_Radial_Distribution.mat'), ...
    'bins', 'mean_pos', 'std_pos', 'mean_neg', 'std_neg', ...
    'dist_pos_all', 'mag_pos_all', 'dist_neg_all', 'mag_neg_all');

fprintf('\n✅ Saved averaged radial magnitude profile in: %s\n', output_dir);
end
