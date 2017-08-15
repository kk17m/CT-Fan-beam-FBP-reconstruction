function [sino, Hk, hn, nn] = fbp_sino_filter(sino, ...
	ds, dsd, window, extra)

[nb,na] = size(sino);

	npad = 2^ceil(log2(2*nb-1)); % padded size
%	printm('nb=%d npad=%d', nb, npad)

sino = [sino; zeros(npad-nb,na)]; % padded sinogram

%[hn,nn] = fbp_ramp(type, npad, ds, dsd);

% ARC
n=npad;
nn = [-(n/2):(n/2-1)]';
h = zeros(size(nn));
h(nn==0) = 1 / (4 * ds^2);
odd = mod(nn,2) == 1;
h(odd) = -1 ./ (pi * dsd * sin(nn(odd) * ds / dsd)).^2;
hn=h;

% FLAT
% nn = [-(n/2):(n/2-1)]';
% h = zeros(size(nn));
% h(n/2+1) = 1 / 4; 
% odd = mod(nn,2) == 1;
% h(odd) = -1 ./ (pi * nn(odd)).^2;
% hn = h / ds^2;


Hk = real(fft(fftshift(hn)));

Hk = Hk .* fbp2_window(npad, window);

Hk = ds * Hk; % differential for discrete-space convolution vs integral

% linear interpolation is like blur with a triangular response,
% so we can compensate for this approximately in frequency domain
% if decon1
 %	Hk = Hk ./ fftshift(nufft_sinc(nn / npad).^2);
% end

sino = ifft( fft(sino, [], 1) .* repmat(Hk, [1 na]), [], 1, 'symmetric'); % apply filter

% trick: possibly keep extra column(s) for zeros!
sino = sino([1:(nb+extra)],:);
sino([(nb+1):(nb+extra)],:) = 0;

function y = nufft_sinc(x)
iz = find(x == 0); % indices of zero arguments
x(iz) = 1;
y = sin(pi*x) ./ (pi*x);
y(iz) = 1;