function VOI= CreateVOI(MRI, Spacing, LinePnts, Ratio)

%Find Center
VOI.Center = round((LinePnts(1,:) + LinePnts(2,:))/2);

%Find Diameter
Diameter = ((LinePnts(1,1)-LinePnts(2,1))*Spacing(1))^2 + ((LinePnts(1,2)-LinePnts(2,2))*Spacing(2))^2 + ((LinePnts(1,3)-LinePnts(2,3))*Spacing(3))^2;
Diameter = sqrt(Diameter);

%Calculate Size of subvolume
if (nargin == 4)
    Ratio = 1.3;
end
VOI.Size = round(Spacing .^ (-1) * Ratio * Diameter) + [2 2 2];

VOI.Start = round(VOI.Center - (VOI.Size ./ 2));
VOI.End = round(VOI.Center + (VOI.Size ./ 2));

% Crop input volumes
%VOI.MRI = MRI(VOI.Start(1):VOI.End(1), VOI.Start(2):VOI.End(2), VOI.Start(3):VOI.End(3));

%Check the boundaries
if (VOI.Start(1)<1)
    VOI.Start(1) = 1;
end
if (VOI.Start(2)<1)
    VOI.Start(2) = 1;
end
if (VOI.Start(3)<1)
    VOI.Start(3) = 1;
end
if (VOI.End(1) > size(MRI,1))
    VOI.End(1) = size(MRI,1);
end
if (VOI.End(2) > size(MRI,2))
    VOI.End(2) = size(MRI,2);
end
if (VOI.End(3) > size(MRI,3))
    VOI.End(3) = size(MRI,3);
end

VOI.MRI = MRI(VOI.Start(1):VOI.End(1), VOI.Start(2):VOI.End(2), VOI.Start(3):VOI.End(3));

%VOI.CI = CI(VOI.Start(1):VOI.End(1), VOI.Start(2):VOI.End(2), VOI.Start(3):VOI.End(3));
VOI.Spacing = Spacing;

% Set border as background seeds
VOI.BGSeeds = false(size(VOI.MRI));

VOI.BGSeeds(1,:,:) = true;
VOI.BGSeeds(size(VOI.BGSeeds,1),:,:) = true;
VOI.BGSeeds(:,1,:) = true;
VOI.BGSeeds(:,size(VOI.BGSeeds,2),:) = true;
VOI.BGSeeds(:,:,1) = true;
VOI.BGSeeds(:,:,size(VOI.BGSeeds,3)) = true;

% We will use matlabs line function to find FG Seeds along the user
% supplied line
LineCrop = 0.15;
VOI.FGSeeds = false(size(VOI.MRI));

FGSeedP1 = LinePnts(1,:)-VOI.Start+[1 1 1];
FGSeedP2 = LinePnts(2,:)-VOI.Start+[1 1 1];

FGSeedP1 = round(FGSeedP1 + (double(FGSeedP2-FGSeedP1).*LineCrop));
FGSeedP2 = round(FGSeedP2 + (double(FGSeedP1-FGSeedP2).*LineCrop));
FGSlice = FGSeedP1(3);

% hold on;axis image;
% colormap(gray);
% hline = line([LinePnts(1,1)-VOI.Start(1)+1 LinePnts(2,1)-VOI.Start(1)+1],[LinePnts(1,2)-VOI.Start(2)+1 LinePnts(2,2)-VOI.Start(2)+1]);
% VOI.XData = get(hline,'XData');
% VOI.YData = get(hline,'YData');
VOI.FGSeeds(:,:,FGSlice) = VOI.FGSeeds(:,:,FGSlice) | CreateLine(VOI.FGSeeds(:,:,FGSlice), [FGSeedP1(1) FGSeedP1(2)],[FGSeedP2(1) FGSeedP2(2)]); 
VOI.FGSeeds(:,:,FGSlice) = VOI.FGSeeds(:,:,FGSlice) | CreateLine(VOI.FGSeeds(:,:,FGSlice), [FGSeedP1(1)+1 FGSeedP1(2)],[FGSeedP2(1)+1 FGSeedP2(2)]); 
VOI.FGSeeds(:,:,FGSlice) = VOI.FGSeeds(:,:,FGSlice) | CreateLine(VOI.FGSeeds(:,:,FGSlice), [FGSeedP1(1)-1 FGSeedP1(2)],[FGSeedP2(1)-1 FGSeedP2(2)]); 
%hfig = figure, imagesc(VOI.MRI(:,:,FGSlice));
%figure,imagesc(VOI.FGSeeds(:,:,FGSlice));
%VOI.Seeds = int8(zeros(size(VOI.MRI)));
VOI.Seeds = int8(VOI.FGSeeds) + 2*int8(VOI.BGSeeds);

function Output = CreateLine(Input, FirstPoint, SecondPoint)
Output = false(size(Input));
Length = norm(SecondPoint-FirstPoint);
step = 1 / (2*Length);
for t=0:step:1
    Pt = round(FirstPoint + (double(SecondPoint-FirstPoint).*t));
    Output(Pt(1),Pt(2))=true;
end


