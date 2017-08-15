%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLE SCRIPT (FAN-BEAM FILTERED BACK-PROJECTION)           %
%                                                             %
% INPUT PARAMETERS:                                           %
%    Log: Sinogram input.                                     %
%    start_angle: Starting angular position in degrees.       % 
%    SOD: Source Object Distance in m.                        %
%    SDD: Source detector distance.                           %
%    Fan_angle: Fan beam opening angle in degrees.            %
%    Norg: Original radial sample size.                       %
%    weighting: Sinogram weighting (Parker or Differential)   %
%    OutputSize: Output size of image in number of pixels     %
%    total_angle: Coverage angle of the source                %
%    Filt: Filter choce (ramp, shepp-logan, cosine,           %
%                        hann, hanning etc.)                  %
%                                                             %
% OUTPUT:                                                     %
%    Reconstruction: Reconstruction output.                   %
%    Bp_RotationIncrement: Angular increments of projections. %
%    Fan_sensor_spacing: Detector Sensor angular spacing.     %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR:                                                     %
%    Kunal Kumar,                                             %
%    Copyright, 2016                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read Sinogram 

clear  
clc
fileid = fopen('Sample_sinogram.sino','rb'); 
proj = fread(fileid,[140,inf],'float32'); 
fclose(fileid);
%proj=fliplr(proj);  % Apply this as per requirement to flip the sinogram

%% Reconstruction Parameters

detector_rows = 1;
SOD = 100; 
SDD = 200;
Fan_angle = 28; 
start_angle = -14; 
total_angle = 180 + Fan_angle; 
OutputSize = 256;
Filter = 'shepp-logan';
weighting = 'parker'; % Parker or differential weighting can be applied.
Norg=size(proj,1);

%% Projection data from logarithmic transformation of intensity ratio

[Max_Intensity,Index] = max(proj(:,:));
X =log(max(ind2sub(size(proj),Max_Intensity)));
Y=log(proj); 
Log=X-Y; 

% Log=proj; % For Pre-Transformed data

%% FFBP Reconstruction

[Reconstruction, Bp_RotationIncrement, Fan_sensor_spacing,Bp_spacing] = ...
                      FFBP_Weighted(Log, start_angle, SOD, SDD, Fan_angle, Norg,...
                                    weighting, OutputSize, total_angle,Filter); 
Reconstruction(isnan(Reconstruction)) = 0;
Reconstruction = Reconstruction*Fan_sensor_spacing*Bp_spacing;

%% Write to file

fileid = fopen('Reconstruction.raw','w+'); 
wrt = fwrite(fileid,Reconstruction,'float32'); 
fclose(fileid);
