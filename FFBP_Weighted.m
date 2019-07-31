function [Reconstruction, Bp_RotationIncrement, Fan_sensor_spacing,Bp_spacing] = ...
                            FFBP_Weighted(Log, start_angle, SOD, SDD, Fan_angle, Norg,...
                                          weighting,OutputSize,total_angle,Filt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAN-BEAM FILTERED BACK-PROJECTION (SHORT SCAN)              %
%                                                             %
% INPUT:                                                      %
%    Log: Sinogram input.                                     %
%    start_angle: Starting angular position in degrees.       % 
%    SOD: Source Object Distance in mm.                       %
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

Log = Log(:,end:-1:1);
out_x = OutputSize;
out_y = size(Log,2);
[in_x,in_y] = size(Log);
SIZE = padarray(Log, [floor((out_x - in_x)/2) floor((out_y - in_y)/2)], 0, 'post');
Log = padarray(SIZE, [ceil((out_x - in_x)/2) ceil((out_y - in_y)/2)], 0, 'pre');

[Npix,Nproj] = size(Log);
Fan_sensor_spacing = (degtorad(Fan_angle)/Norg);
Bp_RotationIncrement = ((degtorad(total_angle))/Nproj);
thetaI = degtorad(0 + start_angle);
thetaF = degtorad(total_angle + start_angle);
Bp_angles = linspace(thetaI, thetaF, Nproj);
Bp_park_ang = linspace(0, degtorad(total_angle), Nproj);
Bp_spacing = abs(Bp_angles(1) - Bp_angles(2));
Det_size = Norg;
Det_angles = Fan_sensor_spacing*linspace(-Det_size/2, Det_size/2, Norg)';
dif = diff(Det_angles(1:2));

if (mod(Npix - Norg, 2) == 0)
    N = (Npix - Norg) / 2;
    dt_in = (fliplr(Det_angles(1) - dif + (0:N-1)*(-dif)))';
    dt_fi = (Det_angles(end) + dif + (0:N-1)*(dif))';
else
    N1 = floor((Npix-Norg)/2);
    N2 = ceil((Npix - Norg)/2);
    dt_in = (fliplr(Det_angles(1) - dif + (0:N1-1)*(-dif)))';
    dt_fi = (Det_angles(end) + dif + (0:N2-1)*(dif))';
end
Det_angles = [dt_in; Det_angles; dt_fi;];

pad(1,1:Nproj) = 1;
DTA = Det_angles*pad;

if strcmp(weighting,'differential')
    
    % Differential weighting
    mat = zeros(Npix,Nproj);
    for i = 1:Nproj
      for j = 1:Npix
        if Bp_angles(i) >= 0 && Bp_angles(i) < 2*(degtorad(Fan_angle/2) + Det_angles(j))
            mat(j,i) = Bp_angles(i)/(2*(degtorad(Fan_angle/2) + Det_angles(j)));
        elseif Bp_angles(i) >= (2*(Det_angles(j) + degtorad(Fan_angle/2))) && Bp_angles(i) < (pi + 2*Det_angles(j))
            mat(j,i) = 1;
        elseif Bp_angles(i) >= (pi + 2*Det_angles(j)) && Bp_angles(i) <= (pi+2*degtorad(Fan_angle/2))
            mat(j,i) = (pi + 2*degtorad(Fan_angle/2) - Bp_angles(i))/(2*(degtorad(Fan_angle/2) - Det_angles(j)));
        else
            mat(j,i) = 0;
        end
      end
    end
    mat=3.*mat.^2-2.*mat.^2;
 
elseif strcmp(weighting,'parker')
    
    % Parker weighting
    mat = zeros(Npix,Nproj);
    for i = 1:Nproj
      for j = 1:Npix
        if Bp_park_ang(i) >= 0 && Bp_park_ang(i) <= 2*(degtorad(Fan_angle/2) - Det_angles(j))
            mat(j,i) = (sin((pi/4)*(Bp_park_ang(i)/((degtorad(Fan_angle/2) - Det_angles(j))))))^2;
        elseif Bp_park_ang(i) >= 2*(degtorad(Fan_angle/2) - Det_angles(j)) && Bp_park_ang(i) <= (pi - 2*Det_angles(j))
            mat(j,i) = 1;
        elseif Bp_park_ang(i) >= (pi - 2*Det_angles(j)) && Bp_park_ang(i) <= (pi + 2*degtorad(Fan_angle/2))
            mat(j,i) = (sin((pi/4)*((pi + 2*degtorad(Fan_angle/2) - Bp_park_ang(i))/((degtorad(Fan_angle/2) + Det_angles(j))))))^2;
        else
            mat(j,i) = 0;
        end
      end
    end
    
else
    disp('Incorrect weighting choice, choose either parker or differential')
   
end

mat(isnan(mat)) = 0;
mat = flipud(mat);

% Gain correction and filtering
[N_req,Filter] = filterProjections(Log, Filt, DTA, SOD, SDD, Fan_sensor_spacing);
Filter = Filter.*mat;

Reconstruction(Npix,Npix) = 0;
RC = (1 + Reconstruction(:,1))*linspace(1, Nproj, Nproj);
MRC = 1 + Reconstruction;
srad = sin(degtorad(Fan_angle/2))*SOD;
Reg_DFOV = linspace(-srad,srad,Npix);
HDFOV = zeros(length(Reg_DFOV),length(Reg_DFOV));
VDFOV = HDFOV;

for i = 1:length(Reg_DFOV)
    HDFOV(:,i) = Reg_DFOV(i);
    VDFOV(:,i) = Reg_DFOV;
end

for n = 1:Nproj
    [H,V] = pol2cart(Bp_angles(n),SOD);
    Hr = V - HDFOV;
    Vr = H - VDFOV;
    Interp = (Bp_angles(n) - (pi/2 - atan2(Vr,Hr)));
    Interp(Interp < -pi) = rem(Interp(Interp < -pi) + 2*pi, 2*pi);
    Interp(Interp >= pi) = rem(Interp(Interp >= pi) - 2*pi, 2*pi);
    Interpolate = interp2(RC, DTA, Filter, MRC*n,Interp);
    Reconstruction = Reconstruction + Interpolate./((Hr).^2 + (Vr).^2);
    figure(1)
    axis image
    colormap(gray(256))
    imagesc(Reconstruction)
    drawnow
end

end
