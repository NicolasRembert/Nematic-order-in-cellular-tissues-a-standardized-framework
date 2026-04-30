function DefectStructure = ComputeDefectLifetime(DefectStructure, numTimePoints, Showplot)
%COMPUTEDEFECTLIFETIME Compute lifetime of tracked defects
%
% DefectStructure = ComputeDefectLifetime(DefectStructure, numTimePoints, Showplot)
%
% Inputs
%   DefectStructure - structure array containing tracked defects
%                     (must include fields .positions and .charge)
%   numTimePoints   - total number of frames in the dataset
%   Showplot        - string ("True"/"False") indicating whether to display
%                     the lifetime histogram
%
% Outputs
%   DefectStructure - updated structure array with an additional field
%                     .lifetime for each defect
%
% The function computes the lifetime of each defect from its trajectory
% and optionally plots the lifetime distributions for positive and
% negative defects.

    % Initialize empty lists to store lifetime values
    positiveLifetimeList = [];
    negativeLifetimeList = [];

    for j = 1:numel(DefectStructure)
        defect = DefectStructure(j);

        % Compute lifetime from the first to the last detected position
        positions = defect.positions;
        indices = find(positions(:,1) ~= 0);

        lifetime = max(indices) - min(indices) + 1;

        DefectStructure(j).lifetime = lifetime;

        % Store the lifetime based on the defect charge
        if defect.charge == 0.5
            positiveLifetimeList = [positiveLifetimeList, lifetime];
        elseif defect.charge == -0.5
            negativeLifetimeList = [negativeLifetimeList, lifetime];
        end
    end

    if Showplot == "True"
        % Plot lifetime histograms for positive and negative defects
        figure()

        histogram(positiveLifetimeList, 'BinWidth', 1, 'FaceColor', 'r');
        xlabel('Lifetime of Defects');
        ylabel('Frequency');
        hold on
        histogram(negativeLifetimeList, 'BinWidth', 1, 'FaceColor', 'b');
        legend('Positive Defect Lifetime Distribution','Negative Defect Lifetime Distribution');
        box off
        set(gca,'fontsize',10)
        pbaspect([1 1 1])
    end
end