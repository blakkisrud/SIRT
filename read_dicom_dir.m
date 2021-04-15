function [imref,IM, Info] = read_dicom_dir(path)
%read_dicom_dir 
%   Modified from "bildlesing_earl_func.m"

% TODO: There seem to be a problem with the info on the extent in the 
% Z-direction (got from a SIRT-PET 0.4384 while in the dcm it says 
% 

d = dir(path); % it has all the dcm files plus the two fist dots. 
info = dicominfo(fullfile(d(3).folder, d(3).name));
IM = zeros([info.Columns, info.Rows, numel(d)-2]);

for i = 3:numel(d)
    fn = fullfile(d(i).folder, d(i).name);
    info = dicominfo(fn) ;
    im = dicomread(fn);
    im = double(im) .* info.RescaleSlope + info.RescaleIntercept;
    % info.RescaleSlope=1, info.RescaleIntercept=-1024
    % output units = m*SV+b m=RescaleSlope, b=RescaleIntersept,
    % SV=Stored Value
    % figure; imagesc(im);% fainetai ligo strabo
    Info(i-2) = info;
    IM(:,:,i-2) = im;
   % IM er 512 x 512 x 125
end

%% lage referense : antar standard referensesystem.
IPP = [Info.ImagePositionPatient];  % 3x125
dZ = (diff(IPP(3,:))); % 1 x124
dZ = mean(dZ); % antar alle like.(men de er ikke!!!)
if dZ < 0,
    IM = flipdim(IM,3);
    dZ = abs(dZ);
end
PE = [Info(1).PixelSpacing', dZ];
XWorldLimits = Info(1).ImagePositionPatient(1);
XWorldLimits = XWorldLimits + [0, PE(2).*(size(IM,2)-1)] + [-PE(2)/2, PE(2)/2];

% vectorize
WorldLimits = Info(1).ImagePositionPatient';
WorldLimits = WorldLimits + [0 0 0; PE.*(size(IM)-1)] + [-PE/2; PE/2];
WorldLimits = WorldLimits';

imref = imref3d(size(IM), WorldLimits(1,:), WorldLimits(2,:), WorldLimits(3,:));

end

