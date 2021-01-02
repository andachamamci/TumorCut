function [Result,CAMap] = tumorsegmentation_L1_new(Image,Seeds,Spacing,parameter,smoothingweight)

% Calculate Strengths and Labels using Cellular Automata
[FCAMap,FCAStr] = castrength(Image, int8(Seeds==1), Spacing, parameter);
[BCAMap,BCAStr] = castrength(Image, int8(Seeds==2), Spacing, single(1.0));
disp('Cellular Automata Done');


DF = -log(FCAStr);
DB = -log(BCAStr);

P = DB ./ (DB+DF);

% Shift strengths [x 1] -> [0 1-x]
%MinStr = min(min(min(CAStr)));
%CAStr = CAStr - MinStr;
%MaxStr = max(max(max(CAStr)));
%CAStr = CAStr ./ MaxStr;
% Flip Bacground strengths around 0
%CAStr(CAMap == 2) = -CAStr(CAMap == 2);
% Convert to uint16 image data to use with level-set
StrIm = uint16(P*2000);
CAMap = (P >= 0.5);
%figure,imagesc(P(:,:,34));
%SliceBrowser(StrIm);
% Use kernel method on strength map
Result = ls3D(StrIm,int8(P>=0.5),Spacing,int32(500), smoothingweight );
