function DefectStructure= CalculateMSDnotTA(DefectStructure, pixelSize, nodeSpacing)
    % Calculate Mean Squared Displacement (MSD) for each defect

    for def = 1:numel(DefectStructure)
        positions = DefectStructure(def).positions;

        lifetime=DefectStructure(def).lifetime;
        numberOfDeltaT = floor(lifetime/2); %# for MSD, dt should be up to 1/4 of number of data points
        msd = zeros(numberOfDeltaT,3); %# We'll store [mean, std, n]

        indices=find(positions(:,1)~=0);

        data=[positions(indices,1),positions(indices,2)];
        
        MSD=zeros(numberOfDeltaT,1);
        Distances=[];
        %# calculate msd for all deltaT's
        for dt = 1:numberOfDeltaT

            deltaCoordsx = data(1+dt,1) - data(1,1);
            deltaCoordsy = data(1+dt,2) - data(1,2);

            squaredDisplacement = (deltaCoordsx*pixelSize*nodeSpacing)^2+(deltaCoordsy*nodeSpacing*pixelSize)^2 ;%# dx^2+dy^2+dz^2;

            MSD(dt)= squaredDisplacement;
            
  
        end
    
   
    

        % Store MSD in the DefectStructure
        DefectStructure(def).MSDnotTA = MSD;
    end
end
