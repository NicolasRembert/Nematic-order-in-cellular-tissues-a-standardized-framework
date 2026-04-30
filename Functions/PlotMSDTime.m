function  [meanMSD,STEMSD,coeffs]= PlotMSDTime(DefectStructure,Charge)
%PLOTMSDTIME Plot mean squared displacement (MSD) over time for a given defect charge
%
% [meanMSD,STEMSD,coeffs] = PlotMSDTime(DefectStructure, Charge)
%
% Inputs
%   DefectStructure - structure array containing tracked defects with fields:
%                     .charge and .MSD
%   Charge          - defect charge to analyze (0.5 or -0.5)
%
% Outputs
%   meanMSD         - mean MSD across all selected defects
%   STEMSD          - standard error of the mean MSD
%   coeffs          - linear fit coefficients of MSD vs time for each defect
%


% Select defects with the requested charge
indices=find([DefectStructure.charge]==Charge);

DefectStructure=DefectStructure(indices);

Numberdefect=length(DefectStructure);

figure();

if Charge == 0.5
    colourlist= [linspace(1, 0.5, Numberdefect)', linspace(0.9, 0.2, Numberdefect)', linspace(0.9, 0.2, Numberdefect)'];
end

if Charge == -0.5
    colourlist= [linspace(0, 0.5, Numberdefect)', linspace(0, 0.8, Numberdefect)', linspace(1, 1, Numberdefect)'];
end

coeffs=zeros(Numberdefect,1);

all_MSD=[];
MSDvalues_list = {};

for num=1:Numberdefect
% Collect MSD values for each defect
% Generate xlist corresponding to the number of Δt used to compute MSD

MSDvalues=DefectStructure(num).MSD;
MSDvalues_list{num}=MSDvalues';

xlist=linspace(1,length(MSDvalues),length(MSDvalues));

xlist=xlist*15;
f=polyfit(xlist,MSDvalues,1);

hold on

coeffs(num)=f(1);

end

maxLength = max(cellfun(@length, MSDvalues_list));

% Concatenate MSD curves by padding with NaNs
for num = 1:length(MSDvalues_list)
    MSDval = MSDvalues_list{num};

    padded_MSD = [MSDval, NaN(1, maxLength - length(MSDval))];

    all_MSD = [padded_MSD; all_MSD];
end

meanMSD=nanmean(all_MSD);
STEMSD=nanstd(all_MSD)/sqrt(length(all_MSD));

xlist=linspace(1,length(meanMSD),length(meanMSD));
xlist=xlist*15;

f=polyfit(xlist,meanMSD,1); 

if Charge == 0.5
    color= [242, 74, 82] / 255;
    disp(Numberdefect+ "positive defects")
end

if Charge == -0.5
    color= [87, 105, 229] / 255;
    disp(Numberdefect+ "negative defects")
end

end