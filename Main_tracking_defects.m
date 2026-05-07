%% DEFECT-FINDER
%% 
% Code by Nicolas Rembert (A . Roux Lab, University of Geneva) adapted from "Nematic length, Defect positions, Defect orientations, inter-defect distances & AVERAGE DEFECTS (from OrientationJ matrices) " by A. Codina and P. Guillamat (X. Trepat Lab, Institute for Bioengineering of Catalonia, IBEC)
% in collaboration with C. Blanch-Mercader (CNRS, Institut Curie)


  %Define directories, in this code the main directory contains
  %subdirectories, these subdirectories contain the csv files from
  %Orientation J result tables. This code will batch analyze all the csv
  %files contained in each subdirectory and save all the variables in .mat
  %file of the subdirectory 

  %Example Main_dir/C2C12/A-B-C

  %Beware the name of the csv files should contain W_" Window size" _
  %G_"Grid size"_ as the first elements of their name 

  %If Dynamic analysis is activated, will track defect and stored all variables in a single mat file , defect density, correlation and order will be saved in a different Mat file 



  %%%% ! ! ! ! ! !

  % Run the coherency to set the feature size prior to running this code

  %%%% ! ! ! ! ! ! 

  % -------------------------------------------------------------------------
% OUTPUT STRUCTURE FOR 1ST PART: Results_nematics
%
% This structure contains the results of the nematic analysis pipeline.
% It groups the analysis parameters together with the measured quantities,
% detected defects, and spatial nematic fields.
%
% PARAMETERS (scalars)
%   FeatureSize       - Scalar. Feature size used in OrientationJ to compute
%                       the local orientation field.
%   Downsampling      - Scalar. Spatial downsampling factor applied to the
%                       orientation field.
%   OrderThreshold    - Scalar. Threshold applied to the local nematic order
%                       parameter to filter low-order regions.
%   MinOrderDistance  - Scalar. Minimum distance used in the computation of
%                       the nematic order parameter.
%
% NUMERICAL RESULTS (vectors, one value per analyzed frame)
%   CorrelationLength - Vector. Nematic correlation length extracted from
%                       the spatial decay of orientation correlations.
%   DefectDensity     - Vector. Density of topological defects in the field.
%
% DEFECT COUNTS (vectors, one value per frame)
%   NumPositive       - Vector. Number of +1/2 defects.
%   NumNegative       - Vector. Number of -1/2 defects.
%   NumTotal          - Vector. Total number of detected defects.
%
% FIELDS (cell arrays, one element per analyzed frame)
%   OrderParameter    - Cell array of matrices. Spatial map of the nematic
%                       order parameter.
%   DefectPositions   - Cell array of NxN matrices. Each row contains the
%                       (x,y) coordinates of detected defects.
%   DefectChargeField - Cell array of matrices. Matrix identifying defect
%                       locations. Values are 0 everywhere except at defect
%                       cores where +0.5 and -0.5 correspond to +1/2 and
%                       -1/2 topological defects.
%   OrientationField  - Cell array of matrices. Local nematic orientation
%                       field (angle values).
%   DefectOrientation - Cell array of vectors. Orientation of +1/2 defects.
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% OUTPUT FILE: Data_dynamic.mat
%
% This file contains the results of the defect dynamics analysis.
%
% Variables are saved using the naming convention:
%
%   CellType_variable_replicate
%
% Example:
%   RPE1_AllDist_1
%
% where:
%   CellType   : name of the analyzed cell type
%   variable   : measured quantity
%   replicate  : index of the movie / experiment
%
% -------------------------------------------------------------------------
% DEFECT DISTANCES
%
%   Same        - Distances between defects with the same topological charge
%                 (+1/2 with +1/2, or -1/2 with -1/2).
%
%   Opposite    - Distances between defects with opposite charges
%                 (+1/2 with -1/2).
%
%   AllDist     - All pairwise defect distances regardless of charge.
%
% -------------------------------------------------------------------------
% DEFECT DENSITY
%
%   DefDensity  - Time series of defect density in the field of view.
%
% -------------------------------------------------------------------------
% DEFECT TRAJECTORY STATISTICS
%
%   meanpos     - Mean time-averaged mean squared displacement (TAMSD)
%                 computed over all +1/2 defects of a replicate.
%
%   meanneg     - Mean time-averaged mean squared displacement (TAMSD)
%                 computed over all -1/2 defects of a replicate.
%
% -------------------------------------------------------------------------
% DEFECT MOTION STATISTICS
%
%   meanparMSD  - Mean squared displacement parallel to the orientation
%                 of +1/2 defects.
%
%   meanperpMSD - Mean squared displacement perpendicular to the
%                 orientation of +1/2 defects.
%
%   meanMSAD    - Mean squared angular displacement of defect orientation.
%
% -------------------------------------------------------------------------
% DEFECT VELOCITIES
%
%   posVel      - Instantaneous velocity magnitudes of +1/2 defects.
%
%   negVel      - Instantaneous velocity magnitudes of -1/2 defects.
%
% -------------------------------------------------------------------------
% DEFECT DISPLACEMENTS
%
%   posDx,posDy - Frame-to-frame displacements of +1/2 defects along x and y.
%
%   negDx,negDy - Frame-to-frame displacements of -1/2 defects along x and y.
%
% -------------------------------------------------------------------------
% DEFECT ORIENTATIONS
%
%   posPhi      - Orientation of +1/2 defects.
%
%   negPhi      - Orientation associated with -1/2 defects.
%
% -------------------------------------------------------------------------
% DEFECT STRUCTURE
%
%   DefStruct   - Structure containing the full defect tracking data,
%                 including defect positions, lifetimes, orientations,
%                 velocities, and trajectory statistics.
%
% -------------------------------------------------------------------------
% DEFECT PAIR ANALYSIS
%
%   posPairs            - Pair statistics for +1/2 defects. Computed as the
%                         dot product between defect orientation vectors
%                         aligned with the vector connecting the two defects.
%
%   pos_opposite_Pairs  - Same quantity computed for +1/2 and -1/2 defect pairs.
%
% -------------------------------------------------------------------------
% DEFECT SHAPE
%
%   Defect_shape - Angular structure of the nematic field around defects.
%                  Stored as {ThetaValues, PhiValues}.
%
% -------------------------------------------------------------------------



addpath("Path to your \Functions")

 
Mainsubdir="Path to \ Confluence 48h cell data";



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_threshold=0.5;

Neighbourhood_is_Downsampling =1; %Keep 1 to have the minimal distance to compute order. Change to a desired value to have a different value than downsampling size. Beware the distance to neighbour for order computer must be larger than Downsampling size. 

Charge_box_size=1; 

Merge_factor=4; %Number of nods to connect defects for COM Aggregation

Pixel_size=0.65;

Dynamic_analysis=1; %Change to 1 to analyse several ime points 

outputFile = fullfile(Mainsubdir,'Data_dynamic .mat'); % <<<< CHANGE path and file name for desired outputfile of the dynamic analysis 

%%% Dynamic analysis Parameters %%% 

timePerPoint = 15;      % minutes, to be adjusted 

maxDistance_tracking=10; %Number of nods to look for defects 1 time point after the other
MaxDefect_disappearance=5; %Number of time point a defect is allowed to disappear and still be tracked 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





fold=dir(Mainsubdir);
fold=fold([fold.isdir]);
fold = fold(~ismember({fold.name}, {'.', '..'}));

for foldin=1:length(fold)

    foldname=fold(foldin).name;

MainDir=fullfile(Mainsubdir, foldname);
  



subfolders = dir(MainDir);

subfolders = subfolders([subfolders.isdir]);

% Remove '.' and '..' entries
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));
%%
for k = 1:length(subfolders)
    % Get the subfolder name
    subfolder_name = subfolders(k).name;
    % Create the full path to the subfolder
    subfolder_path = fullfile(MainDir, subfolder_name);

    
    
    
    % get a list of all files in the folder
    files = dir(fullfile(subfolder_path, '*.csv*'));
    
    % extract the file names and paths
    file_names = {files.name}';
    file_paths = strcat(subfolder_path, filesep, file_names);

    % get a list of all TIF files in the folder
    tif_files = dir(fullfile(subfolder_path, '*.tif*'));
    
    % extract the tif names and paths
    tif_names = {tif_files.name}';
    tif_paths = strcat(subfolder_path, filesep, tif_names);
    
    
    tensorlist= zeros(1,length(file_paths));
    gridlist= zeros(1,length(file_paths));
    qtreshlist= zeros(1,length(file_paths));
    npointlist= zeros(1,length(file_paths));
    minorddistlist= zeros(1,length(file_paths));
    linxcutlist= zeros(1,length(file_paths));
    expxcutlist= zeros(1,length(file_paths));
    defectdenslist= zeros(1,length(file_paths));
    numposlist= zeros(1,length(file_paths));
    numneglist= zeros(1,length(file_paths));
    numalllist= zeros(1,length(file_paths));
    meanQlist= zeros(1,length(file_paths));

    meanCohelist= zeros(1,length(file_paths));

    Exp_fit_list=zeros(150,length(file_paths));
    Exp_x_list=zeros(150,length(file_paths));
    Autocorr_list=zeros(150,length(file_paths));
    Qparam_list={};
    Qmat_list={};

    DummyMat_list={};
    DummyMat2_list={};
    Qmat2_list={};


    OrientationMat_list={};
    Phi_list={};
    Results_nematics={};

for numb = 1:length(file_paths)

    close all
    %Do you want to save the data ? If yes Save_data=1
    SAVE=0;
    Save_data=0;

    %Here we remove the temporary variables relative to each field but we
    %keep the list of variables that we need to compare the different
    %csv files analyzed one after the others
    clearvars -except  Results_nematics MaxDefect_disappearance maxDistance_tracking outputFile timePerPoint Merge_factor Neighbourhood_is_Downsampling S_threshold Pixel_size tif_names tif_paths Dynamic_analysis foldin Qmat2_list DummyMat_list DummyMat2_list fold foldname Mainsubdir k Phi_list OrientationMat_list numb file_paths file_names  MainDir subfolders subfolder_name subfolder_path Qparam_list Qmat_list linxcutlist expxcutlist Autocorr_list numalllist meanQlist tensorlist gridlist qtreshlist npointlist numneglist numposlist numalllist defectdenslist minorddistlist   Dir Exp_fit_list Exp_x_list meanCohelist;
    
   
    
    % display the file names and paths
    disp(file_names);
    disp(file_paths);
    %Dir=strcat(num2str(MainDir),'Outputestcut.csv');



    
    Dir=file_paths{numb};
    disp(Dir)
    File_name_1="Figures_"+file_names{numb};
    [~, name, ~] = fileparts(File_name_1);
    File_name= name;
    %Open and read .csv file
    Array_Stitched=csvread(num2str(Dir),1);
    
    
    Image_path = tif_paths{numb};
    Image = imread(Image_path);
    
  
       
  
    
    %%um/px ratio
    ratio=Pixel_size;
    qThreshold=S_threshold; %Defines the order threshold for topological defect detection, basal value is 0.5 
    %% 
    % Defining variables    
    
    x=Array_Stitched(:,1);
    y=Array_Stitched(:,2);
    nx=Array_Stitched(:,4);
    ny=Array_Stitched(:,5);
    Orientation=atan2(ny,nx);
    Coherence=Array_Stitched(:,7);
    x_um=x*ratio;
    y_um=y*ratio;
    nodeDistance=sqrt((x_um(1)-x_um(2))^2+(y_um(1)-y_um(2))^2);
    
    matches = regexp(File_name, '\d+', 'match');


    Tensor_size= str2double(matches{1});
    Grid_size= str2double(matches{2});

    sizeMATS=sqrt(size(x_um,1));
    
    %node interdistance
    distanceThreshold=10*nodeDistance; %This value will determine the radius of filtration of topological defects 
    
    %vectors to matrices:
    xMat=reshape(x_um,[sqrt(size(x,1)) sqrt(size(x,1))]);
    yMat=reshape(y_um,[sqrt(size(y,1)) sqrt(size(y,1))]);
    CoherenceMat=reshape(Coherence,[sqrt(size(x,1)) sqrt(size(x,1))]);
    OrientationMat_RAW=reshape(Orientation,[sqrt(size(Orientation,1)) sqrt(size(Orientation,1))]);
        
           for an=1:size(Orientation,1)
               if Orientation(an,:)<0
                   Orientation(an,:)=Orientation(an,:)+pi(); %!!!!!!!!!!
               end
           end
           
           OrientationMat=reshape(Orientation,[sizeMATS,sizeMATS]);
    
  

    %% =======================
%% Nematic autocorrelation (FFT version)
%% =======================

tic

theta = OrientationMat;
N = size(theta,1);

% Nematic fields
C = cos(2*theta);
S = sin(2*theta);

% FFT-based autocorrelation
FC = fft2(C);
FS = fft2(S);

ACF = real(ifft2(abs(FC).^2 + abs(FS).^2));
ACF = fftshift(ACF);
ACF = ACF / numel(theta);   % normalization

ACFnorm = ACF;  % already nematic-normalized

toc



%% =======================
%% Distance map
%% =======================

[x,y] = meshgrid(1:N,1:N);
cx = (N+1)/2;
cy = (N+1)/2;

ACFdistance = sqrt((x-cx).^2 + (y-cy).^2) * nodeDistance;

%% =======================
%% Radial averaging
%% =======================

maximACF = max(ACFdistance(:));
distanciesPromig = 0:nodeDistance:maximACF;

ACFpromig = zeros(size(distanciesPromig));
ACFdesv   = zeros(size(distanciesPromig));

for j = 1:length(distanciesPromig)-1
    mask = ACFdistance >= distanciesPromig(j) & ...
           ACFdistance <  distanciesPromig(j+1);
    values = ACFnorm(mask);

    ACFpromig(j) = mean(values);
    ACFdesv(j)   = std(values);
end

%% =======================
%% Characteristic length (threshold method)
%% =======================

distVEC = distanciesPromig(:);
ACFvec  = ACFpromig(:);

xdists = linspace(0, max(distVEC), numel(distVEC)*10);
Interp_Data = interp1(distVEC, ACFvec, xdists, 'linear', 'extrap');

idx = find(Interp_Data < 0.37, 1, 'first');
if isempty(idx)
    expxCut = NaN;
else
    expxCut = xdists(idx);
end

%% =======================
%% Nematic length (linear regression on first points, obsolete)
%% =======================

npoints = 5;

xdata = distanciesPromig(1:npoints);
ydata = ACFpromig(1:npoints);

p = polyfit(xdata, ydata, 1);
mReg = p(1);
nReg = p(2);

x = linspace(0, max(distVEC), 1000);
y = mReg*x + nReg;


% Linear cut
xCut = -nReg / mReg;



    %% 
    % % Order parameter

        if Neighbourhood_is_Downsampling == 1

    MinOrderDistance = Grid_size;

elseif Neighbourhood_is_Downsampling < Grid_size

    disp('Distance can not be inferior to Downsampling size, This distance was set to Downsampling size instead')
    MinOrderDistance = Grid_size;

else

    MinOrderDistance = Neighbourhood_is_Downsampling;

end
   
    % --- Paramètres ---
R = 2; %

% --- Champs nématiques ---
C = cos(2*OrientationMat);
S = sin(2*OrientationMat);

% --- Noyau disque ---
[xg,yg] = meshgrid(-R:R,-R:R);
kernel = (xg.^2 + yg.^2) <= R^2;
kernel = kernel / sum(kernel(:));   % normalisation moyenne

% --- Moyennes locales (convolution rapide) ---
Cmean = conv2(C, kernel, 'same');
Smean = conv2(S, kernel, 'same');

% --- Order parameter ---
Qparameter = sqrt(Cmean.^2 + Smean.^2);
    
    
    figure()
    pcolor(Qparameter)
    shading interp
    colormap gray
    c= colorbar;
    c.Label.String=('Q');
    caxis([0 1])
     space = 1; % space far from refPoint
DummyMat  = zeros(size(OrientationMat));
dummy3    = zeros(size(OrientationMat));

numlowspot = 0;

sizeMATS = size(OrientationMat,1); % assuming square
plotRadius = 5; % half of 10x10 neighborhood

for ii = 1+space : sizeMATS-space
    for jj = 1+space : sizeMATS-space

        if Qparameter(ii,jj) < qThreshold

            dAlpha = [];
            numlowspot = numlowspot + 1;

            %% I. Identify the 8 angles surrounding the reference point
            Alpha1 = OrientationMat(ii+space, jj);
            Alpha2 = OrientationMat(ii+space, jj+space);
            Alpha3 = OrientationMat(ii,      jj+space);
            Alpha4 = OrientationMat(ii-space, jj+space);
            Alpha5 = OrientationMat(ii-space, jj);
            Alpha6 = OrientationMat(ii-space, jj-space);
            Alpha7 = OrientationMat(ii,      jj-space);
            Alpha8 = OrientationMat(ii+space, jj-space);

            %% II. Calculate the differences between the angles
            dAlpha = [Alpha2-Alpha1, Alpha3-Alpha2, Alpha4-Alpha3, Alpha5-Alpha4, ...
                      Alpha6-Alpha5, Alpha7-Alpha6, Alpha8-Alpha7, Alpha1-Alpha8];

            %% III. Shortest path adjustment
            for k = 1:length(dAlpha)
                if abs(dAlpha(k)) <= (pi/2)
                    % no change
                elseif dAlpha(k) < -(pi/2)
                    dAlpha(k) = dAlpha(k) + pi;
                elseif dAlpha(k) > (pi/2)
                    dAlpha(k) = dAlpha(k) - pi;
                end
            end

           
           DummyMat(ii,jj)=sum(dAlpha)/(2*pi);

        else
            DummyMat(ii,jj) = 0;
        end
    end
end


    
    for i=1:length(DummyMat)^2
        if DummyMat(i)==0 
            continue
        elseif DummyMat(i)<0
            DummyMat(i)=-0.5;
        elseif DummyMat(i)>0
            DummyMat(i)=0.5;
        end
    end
    
    matriuALLdefTotals=DummyMat;
    posPOSdef=find(matriuALLdefTotals==+0.5);
    posNEGdef=find(matriuALLdefTotals==-0.5);  %%%%%%%%%%%%%%%%%%%%%%%%!!!
    posPOSdef2=[];
    posNEGdef2=[];
    


%% COM Aggregation (fusion of close defects of same charge)

matriuALLdefTotals = DummyMat;

posPOSdef = find(matriuALLdefTotals == +0.5);
posNEGdef = find(matriuALLdefTotals == -0.5);

% fusion des défauts
allDef = [posPOSdef(:); posNEGdef(:)];

% ==========================================================
% Construction positions réelles et charges des défauts
% ==========================================================

if isempty(allDef)

    posPOSdef2 = [];
    posNEGdef2 = [];
    posPOSx = [];
    posPOSy = [];
    posNEGx = [];
    posNEGy = [];

else

    % positions réelles
    xDef = xMat(allDef);
    yDef = yMat(allDef);

    signDef = matriuALLdefTotals(allDef);

    % ==========================================================
    % taille domaine (PBC)
    % ==========================================================

    xmin = min(xMat(:));
    ymin = min(yMat(:));

    Lx = max(xMat(:)) - xmin;
    Ly = max(yMat(:)) - ymin;

    r_merge = Merge_factor * Grid_size;
    r_merge2 = r_merge^2;

    N = length(allDef);
    used = zeros(N,1);

    newPos = [];
    newSign = [];

    for i = 1:N

        if used(i)
            continue
        end

        xi = xDef(i);
        yi = yDef(i);
        si = signDef(i);

        xs = xi;
        ys = yi;
        ss = si;

        used(i) = 1;

        for j = i+1:N

            if used(j)
                continue
            end

            xj = xDef(j);
            yj = yDef(j);
            sj = signDef(j);

            % fusion seulement si même signe
            if sign(si) ~= sign(sj)
                continue
            end

            % distance PBC
            dxi = xj - xi;
            dyi = yj - yi;

            dxi = dxi - Lx * round(dxi/Lx);
            dyi = dyi - Ly * round(dyi/Ly);

            dist2 = dxi^2 + dyi^2;

            if dist2 < r_merge2

                xs(end+1) = xi + dxi;
                ys(end+1) = yi + dyi;
                ss(end+1) = sj;

                used(j) = 1;

            end
        end

        % centre de masse
        x_c = sum(xs .* abs(ss)) / sum(abs(ss));
        y_c = sum(ys .* abs(ss)) / sum(abs(ss));

        % charge moyenne
        s_mean = sum(ss) / length(ss);

        % wrap PBC avec origine correcte
        x_c = mod(x_c - xmin, Lx) + xmin;
        y_c = mod(y_c - ymin, Ly) + ymin;

        % charge finale
        s = 0.5 * sign(s_mean);

        % projection sur la grille
        dist2 = (xMat(:) - x_c).^2 + (yMat(:) - y_c).^2;
        [~, idxCM] = min(dist2);

        newPos(end+1) = idxCM;
        newSign(end+1) = s;

    end

    posPOSdef2 = unique(newPos(newSign > 0));
    posNEGdef2 = unique(newPos(newSign < 0));

    % indices ligne/colonne
    [posPOSx,posPOSy] = ind2sub([sizeMATS sizeMATS],posPOSdef2);
    [posNEGx,posNEGy] = ind2sub([sizeMATS sizeMATS],posNEGdef2);

end

   
    
    %Make graph
    field_defects=figure('visible','off');
    figure(field_defects)
    quiver(x_um,y_um,nx.*Coherence,ny.*Coherence,"k","ShowArrowHead","off")
    hold on
    scatter(x_um(posPOSdef2),y_um(posPOSdef2),"Filled","r","LineWidth",1)
    scatter(x_um(posNEGdef2),y_um(posNEGdef2),"Filled","b","LineWidth",1)

    posmax= x(end)

    xlim([0 posmax])
    ylim([0 posmax])
    axis square
    hold off

%     

    
   
    % Defect densities

      numPOS=length(posPOSdef2);
      numNEG=length(posNEGdef2);
       numALL=numPOS+numNEG;
         Area=max(max(xMat))*max(max(yMat));
         
          DefectDensity=numALL/Area;

    
    
    %% 
    % Defect orientation
    
    %Matrix with values +/-0.5 for positive/negative defects
    
    Qmat=zeros([1,numel(xMat)]);

    %Make sure there are defects here 
    if isempty(find(isnan(posPOSdef2)))

    for i=1:length(posPOSdef2)
        Qmat(posPOSdef2(i))=0.5;
    end

    end
    if isempty(find(isnan(posNEGdef2)))
    for i=1:length(posNEGdef2)
        Qmat(posNEGdef2(i))=-0.5;
    end
    end
    Qmat=reshape(Qmat,sqrt(numel(Qmat)),sqrt(numel(Qmat)));

    for i=1:length(OrientationMat)
        for j=1:length(OrientationMat)
            AlphaREF=OrientationMat(i,j);
            COHREF=CoherenceMat(i,j);
            QxxMat(i,j)=2*COHREF*cos(2*AlphaREF);
            QxyMat(i,j)=2*COHREF*sin(2*AlphaREF);
            QyyMat(i,j)=-2*COHREF*cos(2*AlphaREF);
        end
    end
    
    delta=1;
    
    [dyQxx,dxQxx] = gradient(QxxMat,delta);
    [dyQxy,dxQxy] = gradient(QxyMat,delta);
    [dyQyy,dxQyy] = gradient(QyyMat,delta);
    
    PHI=[];
    
    for i=1:length(OrientationMat)
        for j=1:length(OrientationMat)
            if Qmat(i,j)<0 %% for negative defects
                PHI(i,j)=-atan2((-dxQxy(i,j)-dyQxx(i,j)),(dxQxx(i,j)-dyQxy(i,j)))/3;
            elseif Qmat(i,j)>0 %% for positive defects
                PHI(i,j)=atan2((dxQxy(i,j)-dyQxx(i,j)),(dxQxx(i,j)+dyQxy(i,j)));
            else
                PHI(i,j)=NaN; 
            end
        end
    end
    
    Orientation=figure('visible','on');
    figure(Orientation);
    quiver(x_um,y_um,nx.*Coherence,ny.*Coherence,"k","ShowArrowHead","off")
    hold on
    
    if isempty(find(isnan(posPOSdef2)))
    scatter(x_um(posPOSdef2),y_um(posPOSdef2),"Filled","r","LineWidth",1)
    numPOS=numel(posPOSx); 
    end
      if isempty(find(isnan(posNEGdef2)))
    scatter(x_um(posNEGdef2),y_um(posNEGdef2),"Filled","b","LineWidth",1)
    numNEG=numel(posNEGx);
      end
    
    

    if isempty(find(isnan(posPOSdef2)))

    for p=1:numPOS
            pxi=posPOSx(:,p);
            pyi=posPOSy(:,p);
    
    quiver(xMat(pxi,pyi),yMat(pxi,pyi),cos(PHI(pxi,pyi)),sin(PHI(pxi,pyi)),"r","ShowArrowHead","off","LineWidth",2,"AutoScale","on","AutoScaleFactor",50)
    end
    
    for n=1:numNEG
            nxi=posNEGx(:,n);
            nyi=posNEGy(:,n);
    
    quiver(xMat(nxi,nyi),yMat(nxi,nyi),cos(PHI(nxi,nyi)),sin(PHI(nxi,nyi)),"b","ShowArrowHead","off","LineWidth",2,"AutoScale","on","AutoScaleFactor",50)
    quiver(xMat(nxi,nyi),yMat(nxi,nyi),cos(PHI(nxi,nyi)+2*pi/3),sin(PHI(nxi,nyi)+2*pi/3),"b","ShowArrowHead","off","LineWidth",2,"AutoScale","on","AutoScaleFactor",50)
    quiver(xMat(nxi,nyi),yMat(nxi,nyi),cos(PHI(nxi,nyi)+4*pi/3),sin(PHI(nxi,nyi)+4*pi/3),"b","ShowArrowHead","off","LineWidth",2,"AutoScale","on","AutoScaleFactor",50)
    end
     
    axis square


    %Make graph for PTV ( positive defects position without the quiver)
    positive_defects_PTV=figure('visible','on');
    figure(positive_defects_PTV)
    hold on
    scatter(x_um(posPOSdef2),y_um(posPOSdef2),"Filled","r","LineWidth",1)
    
    numPOS=numel(posPOSx); numNEG=numel(posNEGx);
    hold on
    for p=1:numPOS
            pxi=posPOSx(:,p);
            pyi=posPOSy(:,p);
    
    quiver(xMat(pxi,pyi),yMat(pxi,pyi),cos(PHI(pxi,pyi)),sin(PHI(pxi,pyi)),"r","ShowArrowHead","off","LineWidth",2,"AutoScale","on","AutoScaleFactor",50)
    hold on
    end
    
    

    %posmax= x(end)

    %xlim([0 posmax])
    %ylim([0 posmax])
    set(gca,'box','off','ycolor','w')
    ax1 = gca;                   % gca = get current axis
    ax1.YAxis.Visible = 'off';   % remove y-axis
    ax1.XAxis.Visible = 'off';   % remove x-axis
    
    axis square
    hold off

    %Make graph for PTV ( negative defects position without the quiver)
    negative_defects_PTV=figure('visible','on');
    figure(negative_defects_PTV)
    hold on
    
    scatter(x_um(posNEGdef2),y_um(posNEGdef2),"Filled","b","LineWidth",1)
    numPOS=numel(posPOSx); numNEG=numel(posNEGx);
    hold on
    
    for n=1:numNEG
            nxi=posNEGx(:,n);
            nyi=posNEGy(:,n);
    hold on
    quiver(xMat(nxi,nyi),yMat(nxi,nyi),cos(PHI(nxi,nyi)),sin(PHI(nxi,nyi)),"b","ShowArrowHead","off","LineWidth",2,"AutoScale","on","AutoScaleFactor",50)
    quiver(xMat(nxi,nyi),yMat(nxi,nyi),cos(PHI(nxi,nyi)+2*pi/3),sin(PHI(nxi,nyi)+2*pi/3),"b","ShowArrowHead","off","LineWidth",2,"AutoScale","on","AutoScaleFactor",50)
    quiver(xMat(nxi,nyi),yMat(nxi,nyi),cos(PHI(nxi,nyi)+4*pi/3),sin(PHI(nxi,nyi)+4*pi/3),"b","ShowArrowHead","off","LineWidth",2,"AutoScale","on","AutoScaleFactor",50)
    hold on
    end

    %posmax= x(end)

    %xlim([0 posmax])
    %ylim([0 posmax])
    set(gca,'box','off','ycolor','w')
    ax1 = gca;                   % gca = get current axis
    ax1.YAxis.Visible = 'off';   % remove y-axis
    ax1.XAxis.Visible = 'off';   % remove x-axis
    
    axis square
    hold off
    
    
end
    
    

  
    %% 
    
    %This part finds the order parameter value at the position of the
    %defect


    % Saving figures and variables
    
    meanQ=mean(Qparameter(:));

    tensorlist(numb)= Tensor_size;
    gridlist(numb)= Grid_size;
    qtreshlist(numb)= qThreshold;
   
    minorddistlist(numb)=MinOrderDistance;
    
    defectdenslist(numb)=DefectDensity;
    numposlist(numb)=numPOS;
    numneglist(numb)=numNEG;
    numalllist(numb)=numALL;
    meanQlist(numb)=meanQ;
    meanCohelist(numb)=mean(CoherenceMat(:));

    Qparam_list{numb}=Qparameter;
    Qmat_list{numb}=Qmat;


    DummyMat_list{numb}=DummyMat;

    OrientationMat_list{numb}=OrientationMat;
    Phi_list{numb}=PHI;


    %Structure that saves all the variables 

Results_nematics.Parameters.FeatureSize = tensorlist;
Results_nematics.Parameters.Downsampling = gridlist;
Results_nematics.Parameters.OrderThreshold = qtreshlist;
Results_nematics.Parameters.MinOrderDistance = minorddistlist;

Results_nematics.Statistics.CorrelationLength = expxcut_list;
Results_nematics.Statistics.DefectDensity = defectdenslist;

Results_nematics.Defects.NumPositive = numposlist;
Results_nematics.Defects.NumNegative = numneglist;
Results_nematics.Defects.NumTotal = numalllist;

Results_nematics.Fields.OrderParameter = Qparam_list;
Results_nematics.Fields.DefectChargeField = Qmat_list;
Results_nematics.Fields.OrientationField = OrientationMat_list;
Results_nematics.Fields.DefectOrientation = Phi_list;


    
 
   %To decomment 

    npointlist(numb)= npoints;

    linxcutlist(numb)=xCut;
    expxcutlist(numb)=expxCut;

    Autocorr_list(1:length(distVEC),numb)=ACFvec;

    Exp_fit_list(1:length(Interp_Data),numb)=Interp_Data;

    Exp_x_list(1:length(Interp_Data),numb)=xdists(1:length(Interp_Data));

    % =========================
% LOAD IMAGE
% =========================

img = imread(Image_path);
[img_h,img_w,~] = size(img);

[Q_rows,Q_cols] = size(Qmat);

% =========================
% DEFECT DETECTION
% =========================

[pos_y_def,pos_x_def] = find(Qmat == 0.5 | Qmat == -0.5);

border_margin = 2;

valid = pos_x_def > border_margin & ...
        pos_x_def < Q_cols-border_margin & ...
        pos_y_def > border_margin & ...
        pos_y_def < Q_rows-border_margin;

pos_x_def = pos_x_def(valid);
pos_y_def = pos_y_def(valid);

% =========================
% GRID → IMAGE COORDINATES
% =========================

x_img = linspace(1,img_w,Q_cols);
y_img = linspace(1,img_h,Q_rows);

% =========================
% PLOT IMAGE
% =========================

figure('visible','on','Color','w')
imshow(img)
hold on

arrowScale = 30;


 % --- Champ d'orientation sous-échantillonné ---
            step = 10;

            [X,Y] = meshgrid((1:Q_cols)*Grid_size,(1:Q_rows)*Grid_size);
            U = cos(OrientationMat)';
            V = sin(OrientationMat)';

            Xq = X(1:step:end,1:step:end);
            Yq = Y(1:step:end,1:step:end);
            Uq = U(1:step:end,1:step:end);
            Vq = V(1:step:end,1:step:end);

            quiver(Xq, Yq, Uq, Vq, 0.8, ...
                'Color',[0.8 0.8 0], ...
                'LineWidth',1, ...
                'ShowArrowHead','off');

% =========================
% DEFECTS
% =========================

for i = 1:length(pos_x_def)

    x = x_img(pos_x_def(i));
    y = y_img(pos_y_def(i));

    charge = Qmat(pos_y_def(i),pos_x_def(i));
    phi_def = PHI(pos_y_def(i),pos_x_def(i));

    if charge > 0

        scatter(y,x,220,'r','filled')

        quiver(y,x,...
            cos(phi_def),sin(phi_def),...
            arrowScale,'r','LineWidth',7,...
            'ShowArrowHead','off')

    else

        scatter(y,x,220,'b','filled')

        for m = 0:2

            quiver(y,x,...
                cos(phi_def+m*2*pi/3),...
                sin(phi_def+m*2*pi/3),...
                arrowScale,'b','LineWidth',7,...
                'ShowArrowHead','off')

        end
    end
end

axis off

% =========================
% SAVE IMAGE
% =========================

fig_folder = fullfile(subfolder_path, 'Figures_Defects');

% vérifier si le dossier existe, sinon le créer
if ~exist(fig_folder, 'dir')
    mkdir(fig_folder);
end

save_name = fullfile(fig_folder, strcat('Defects_', file_names{numb}, '.png'));
print(gcf, save_name, '-dpng', '-r300')

close

    
    
end

 % Define the name of the .mat file
    subfolder_name = [subfolder_name '_defects'];
    %subfolder_name = [subfolder_name];
    matfile_name = fullfile(subfolder_path, [subfolder_name, '.mat']);
    
    % Save all variables from the current workspace to the .mat file
    save(matfile_name);
end 
end

%% 

if Dynamic_analysis== 1

% =========================================
%  Master Script: Defect Analysis (tracking) 
%  =========================================

clearvars -except Mainsubdir outputFile Dynamic_analysis
 clc; close all;

% -------------------------------
% Main directory containing cell type subfolders equal to prebious maindir
% -------------------------------
MainDir=Mainsubdir;


% -------------------------------
% Ignore list if needed 
% -------------------------------
 ignoreCells = {'MDCK','U2OS','hBEC','SKMEL','Low density'};
 %ignoreCells = {'MDCK','U2OS','hBEC','SKMEL','NIH3T3','L6','CAF','NAF','H9C2','RPE1','test','C2C12'};

% -------------------------------
% Families & Colors, add any new cell type to its family, name should
% correspond to folder name 
% -------------------------------
families = struct( ...
    'Epithelial',  {{'MDCK','RPE1','U2OS','hBEC'}}, ...
    'Myoblasts',   {{'BAOSMC','C2C12','H9C2','L6'}}, ...
    'Fibroblasts', {{'CAF','NAF','NIH3T3'}}, ...
    'Other',       {{'SKMEL'}} ...
);

baseColors.Myoblasts   = [0 0.4470 0.7410];     
baseColors.Epithelial  = [0.8500 0.3250 0.0980];
baseColors.Fibroblasts = [0.4660 0.6740 0.1880];
baseColors.Other       = [0.4940 0.1840 0.5560];


cmap = createColorMap(families, baseColors);


% -------------------------------
% Parameters to be adjusted 
% -------------------------------
     


NewMat = 0; % For first pass detection, please leave this value to 0 if you want to remove the defects which lifetime is < 3 , otherwise set to 1 

% ==============================
%  Phase 1: Analysis
% ==============================

Sub_maindir = dir(MainDir);
Sub_maindir = Sub_maindir([Sub_maindir.isdir]);
Sub_maindir = Sub_maindir(~ismember({Sub_maindir.name},{'.','..','outputdir'}));

S = struct();

for sub = 1:length(Sub_maindir)
    cellType = Sub_maindir(sub).name;
    if ismember(cellType,ignoreCells)
        fprintf('Skipping ignored cell type: %s\n', cellType);
        continue;
    end

fprintf('\n=== Processing Cell type: %s ===\n', cellType);
subPath = fullfile(MainDir, cellType);

% List all subfolders (first-level subfolders inside the cell type)
subfolders = dir(subPath);
subfolders = subfolders([subfolders.isdir]);
subfolders = subfolders(~ismember({subfolders.name},{'.','..'}));

% Only keep subfolders with a letter followed by a number
pattern = '^[A-Za-z]\d+';
isValid = ~cellfun(@isempty, regexp({subfolders.name}, pattern));
subfolders = subfolders(isValid);

% Now look for .mat files inside only these filtered subfolders
matFiles = [];
for k = 1:length(subfolders)
    f = subfolders(k);
    files = dir(fullfile(f.folder, f.name, '*defects.mat'));
    matFiles = [matFiles; files];
end

    for m = 1:length(matFiles)
        matPath = fullfile(matFiles(m).folder, matFiles(m).name);
        fprintf('   Found file: %s\n', matPath);
        movieIdx = m;
        NewMat = 0;

        data = load(matPath);
        close all
        if ~isfield(data,'Qmat_list')
            fprintf('   ⚠️ Missing Qmat_list in %s → skipping.\n', matPath);
            continue;
        end
        Qmat_list = data.Qmat_list;
        Phi_list= data.Phi_list;
        Grid_size=data.Grid_size;
        Pixel_size=data.Pixel_size;
        timePerPoint=data.timePerPoint;
        OrientationMat_list=data.OrientationMat_list;
        
        
        OrientationMat=OrientationMat_list{1};

        maxDistance_tracking=data.maxDistance_tracking; 
MaxDefect_disappearance=data.MaxDefect_disappearance;


% Here we will filter out the defects that are less than 3 nodes away to the edges of the image and track them using a  Hungarian algorithm 
margin = 3; % number of nodes to exclude from each border, to avoid noise from edges 


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

%In this first phase we track defects and remove from the position matrix
%all defects which lifetime is inferior to 3 
        if NewMat==0
            CleanQmat_list = Qmat_list;
            Phimatrices=Phi_list;
            Cell_of_Matrices = CleanQmat_list;
            numTimePoints = length(Cell_of_Matrices);
         
            DefectData = BuildDefectData(Cell_of_Matrices, Phi_list, maxDistance_tracking, MaxDefect_disappearance); % Here 10 is the number of nod distance to which the cost is computed, 5 is the number of time point the defect can disappeared and is kept alive, hence can still be tracked afterwards 
            defects=DefectData.defects;

            defects = ComputeDefectLifetime(defects,numTimePoints,"False");
            defects = CalculateMSD(defects,Pixel_size,Grid_size);
            defects = ComputeAvgOppositeDistance(defects,Cell_of_Matrices);
            defects = ComputeAvgSameChargeDistance(defects,Cell_of_Matrices);

            CopyQmat_list = Remove_Positions(Qmat_list,defects);
            NewMat = 1; % switch for second pass
        end

        % ---- Phase B: NewMat=1 (final analysis) ----
        if NewMat==1
            CleanQmat_list = CopyQmat_list;
            Cell_of_Matrices = CleanQmat_list;
            numTimePoints = length(Cell_of_Matrices);
           
           

            DefectData = BuildDefectData(Cell_of_Matrices, Phi_list, maxDistance_tracking, MaxDefect_disappearance);
            defects=DefectData.defects;
            defects = ComputeDefectLifetime(defects,numTimePoints,"False");
            defects = CalculateMSD(defects,Pixel_size,Grid_size);
            defects = ComputeAvgOppositeDistance(defects,Cell_of_Matrices);
            defects = ComputeAvgSameChargeDistance(defects,Cell_of_Matrices);
            defects = CalculateInstantVelocity(defects, Pixel_size, Grid_size,timePerPoint);

            pairs = pairPositiveDefects_phiproduct(defects);
            %stats = analyzeDefectCreationDestruction(defects);

            defects = CalculateAngularMSD(defects);
            defects = CalculateMSD_ParPerp(defects,Pixel_size, Grid_size);

            defects= CalculateMSDnotTA(defects, Pixel_size, Grid_size);

            pair_opposite = pairOppositeDefects_phiproduct(defects);
        [meanMSD_par, meanMSD_perp, meanMSAD,coeffs] = ComputeMeanMSD_MSAD_Positive(defects, 15);

            % Distances
            numdefects = length(defects);
            Same=[]; Opposite=[]; AllDist=[]; All_pos_dx=[];All_pos_dy=[]; All_neg_dx=[];All_neg_dy=[]; All_pos_phi=[]; All_neg_phi=[];
            posVel = [];   % pour charge +0.5
            negVel = [];   % pour charge -0.5
            for d=1:numdefects
                dsame = defects(d).AvgSameChargeDist;
                dopp  = defects(d).AvgOppositeDist;
                % dall  = [dsame(:)'; dopp(:)'];

                             dall=[];
                   for j=1:length(dsame)
                       dall(j)=nanmin(dsame(j),dopp(j));
            
                   end
                Same     = [Same; dsame(:)'];
                Opposite = [Opposite; dopp(:)'];
                AllDist  = [AllDist; dall];

                v = defects(d).meanVelocity_um_per_h;     % velocity vectors in µm/h

        dispx=defects(d).dx_inst_um;
        dispy=defects(d).dy_inst_um;
        phidef=defects(d).phi;

       

    if defects(d).charge == 0.5
        posVel = [posVel; v(:)'];
        All_pos_dx= [All_pos_dx; dispx(:)'];
        All_pos_dy= [All_pos_dy; dispy(:)'];
        All_pos_phi=[All_pos_phi; phidef(:)'];

    elseif defects(d).charge == -0.5
        negVel = [negVel; v(:)'];
        All_neg_dx= [All_neg_dx; dispx(:)'];
        All_neg_dy= [All_neg_dy; dispy(:)'];
         All_neg_phi=[All_neg_phi; phidef(:)'];
    end

            end
            Same     = Same * Pixel_size*Grid_size;
            Opposite = Opposite * Pixel_size*Grid_size;
            AllDist  = AllDist * Pixel_size*Grid_size;

         

            % MSD
            meanpos  = PlotMSDTime(defects, 0.5);
            meanneg  = PlotMSDTime(defects, -0.5);

           numTimePoints = length(Cell_of_Matrices);
GlobalNumDefects = zeros(1, numTimePoints);

for t = 1:numTimePoints
    NumberofDefects = nnz(Qmat_list{t});
    GlobalNumDefects(t) = NumberofDefects;
end

            ImageSize = length(OrientationMat)^2 * Grid_size^2 * Pixel_size^2 ; % area for defect density
            DefDensity = GlobalNumDefects / ImageSize;



            %Detectionof defect position and placement in lists for further
            %analysis 

 Qmat = CleanQmat_list{t};
    % ensure Qmat is square
    matrixSize = size(Qmat,1);

    tmp_pos_x = [];
    tmp_pos_y = [];
    tmp_neg_x = [];
    tmp_neg_y = [];

    for yy = 1:matrixSize
        for xx = 1:matrixSize
            charge = Qmat(xx,yy); % follow your indexing convention
            if xx <= margin || xx > matrixSize-margin || ...
               yy <= margin || yy > matrixSize-margin
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
            %


            [All_u_ori, All_v_ori, All_xBox_ori, All_yBox_ori] = analyze_defects_shape(posPOSx,posPOSy,posNEGx,posNEGy,OrientationMat_list,Phi_list,'pos'); %Extract cells which contain a box of vetors extarcted around the center of each positive defects 

            [Xq,Yq,nematicMat,gridStep,count_ori] = average_nematic_shape(All_u_ori,All_v_ori,All_xBox_ori,All_yBox_ori,1); %Averages all the rotated defect and realigns the averaged matrix on the same grid of coordinate 

            [ThetaValues, PhiValues, ThetaNeg, PhiNeg] = AngleDistributionv2_withCoords(nematicMat,Xq,Yq,0,0); %Computes the shape of defect and distribution of Theta and Phi , ThetaNeg, PhiNeg contain negative distribution of angles while ThetaValues, PhiValues contain both positive and negative PhiValues  

           
% =========================================================
% SAUVEGARDE DANS LA STRUCTURE S
% =========================================================

S.(sprintf('Same_%s_%d',     cellType, movieIdx)) = Same;
S.(sprintf('Opposite_%s_%d', cellType, movieIdx)) = Opposite;
S.(sprintf('AllDist_%s_%d',  cellType, movieIdx)) = AllDist;

S.(sprintf('meanpos_%s_%d', cellType, movieIdx)) = meanpos;
S.(sprintf('meanneg_%s_%d', cellType, movieIdx)) = meanneg;

S.(sprintf('meanparMSD_%s_%d',  cellType, movieIdx)) = meanMSD_par;
S.(sprintf('meanperpMSD_%s_%d', cellType, movieIdx)) = meanMSD_perp;
S.(sprintf('meanMSAD_%s_%d',     cellType, movieIdx)) = meanMSAD;

S.(sprintf('DefDensity_%s_%d', cellType, movieIdx)) = DefDensity;
S.(sprintf('posVel_%s_%d',  cellType, movieIdx)) = posVel;
S.(sprintf('negVel_%s_%d',  cellType, movieIdx)) = negVel;

S.(sprintf('negDx_%s_%d',  cellType, movieIdx)) = All_neg_dx;
S.(sprintf('negDy_%s_%d',  cellType, movieIdx)) = All_neg_dy;
S.(sprintf('posDx_%s_%d',  cellType, movieIdx)) = All_pos_dx;
S.(sprintf('posDy_%s_%d',  cellType, movieIdx)) = All_pos_dy;
S.(sprintf('posPhi_%s_%d', cellType, movieIdx)) = All_pos_phi;
S.(sprintf('negPhi_%s_%d', cellType, movieIdx)) = All_neg_phi;
S.(sprintf('DefStruct_%s_%d', cellType, movieIdx)) = defects;



            S.(sprintf('posPairs_%s_%d', cellType, movieIdx)) = pairs;
            %S.(sprintf('posStats_%s_%d', cellType, movieIdx)) = stats %Work In Progress

            S.(sprintf('pos_opposite_Pairs_%s_%d', cellType, movieIdx)) = pair_opposite;

  
S.(sprintf('Defect_shape_%s_%d', cellType, movieIdx)) = {ThetaValues, PhiValues};


            % Save once all done
if isfile(outputFile)
    save(outputFile,'-struct','S','-append'); % append new fields
else
    save(outputFile,'-struct','S'); % first save
end
        end
    end
end

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
