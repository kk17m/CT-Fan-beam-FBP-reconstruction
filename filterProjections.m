function [p,sino,Hk,hn] = filterProjections(p_in, filter, DTA,SOD,Fan_sensor_spacing)

p_in = p_in.*(cos(DTA)*SOD);
p = p_in; 

% Filtering the sinogram 
% (Jeffrey A Fessler, "Michigan Image Reconstruction Toolbox," 
%  http://web.eecs.umich.edu/~fessler/code/index.html).
[sino,Hk,hn,nn] = fbp_sino_filter(p, degtorad(Fan_sensor_spacing),...
                                  200, filter, 0);

% Following functions are taken fron iradon function
% Copyright 1993-2013 The MathWorks, Inc.
% Design the filter
len = size(p,1);
H = designFilter(filter, len, 1);

if strcmpi(filter, 'none')
return;
end

p(length(H),1)=0;                       % Zero pad projections
p = fft(p);                             % p holds fft of projections

for i = 1:size(p,2)
     p(:,i) = p(:,i).*H.*0.5*...
                (Fan_sensor_spacing/...
                  sin(Fan_sensor_spacing)^2); % frequency domain filtering
end

p = real(ifft(p));                      % p is the filtered projections
p(len+1:end,:) = [];                    % Truncate the filtered projections

end 
