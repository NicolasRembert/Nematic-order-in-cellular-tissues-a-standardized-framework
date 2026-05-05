function Qmat_list_updated = Remove_Positions(Qmat_list, DefectStructure)
    % Inputs:
    % Qmat_list: A cell array of matrices
    % DefectStructure: A structure containing defect properties, including:
    %   - lifetime: A list of defect lifetimes
    %   - positions: A 2D matrix (Nx2) where each row is [index, y] coordinates
    
    % Initialize the updated Qmat_list as a copy of the original
    Qmat_list_updated = Qmat_list;

    % Loop through each defect in the structure
    for idx = 1:length(DefectStructure)
        % Check if the defect's lifetime is less than 4
        if DefectStructure(idx).lifetime < 4
            % Extract defect positions
            defect_positions = DefectStructure(idx).positions;

            % Loop through all positions
            for pos = 1:size(defect_positions, 1)
                x = defect_positions(pos, 1); % Row index
                y = defect_positions(pos, 2); % Column index

                % Check if x and y are non-zero
                if x ~= 0 && y ~= 0
                    % Validate that the index corresponds to a valid Qmat
                    if pos > 0 && pos <= length(Qmat_list_updated)
                        % Extract the specific Qmat
                        Qmat = Qmat_list_updated{pos};

                        % Validate that x and y are within bounds of Qmat
                        if x > 0 && x <= size(Qmat, 1) && y > 0 && y <= size(Qmat, 2)
                            % Check and update the value if non-zero
                            if Qmat(x, y) ~= 0
                                Qmat(x, y) = 0;
                            end
                        end

                        % Update the Qmat in the cell array
                        Qmat_list_updated{pos} = Qmat;
                    end
                end
            end
        end
    end
end
