function window = fbp2_window(n, window)
%|function window = fbp2_window(n, window)
%| compute an apodizing window of length n and fft shift it

if nargin == 1 && streq(n, 'test'), fbp2_window_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

if ischar(window)
	if isempty(window) || streq(window, 'boxcar') || streq(window, 'ramp')
		window = ones(n,1);

	elseif streq(window, 'hamming,', 8)
		cut = sscanf(window, 'hamming,%g');
		window = my_hamming(n, cut);
        
        elseif streq(window, 'shepp-logan')
				window = my_shepp(n, 1.0);
                
%                 elseif streq(window, 'ramp')
% 				window = my_ramp(n, 1.0);

                elseif streq(window, 'cosine')
				window = my_cosine(n, 1.0);
                
	elseif streq(window, 'hanning,', 8)
		cut = sscanf(window, 'hanning,%g');
		window = my_hann(n, cut);

	elseif streq(window, 'hann')
		window = my_hann(n, 1.0);
	elseif streq(window, 'hann50')
		window = my_hann(n, 0.5);
	elseif streq(window, 'hann75')
		window = my_hann(n, 0.75);
	elseif streq(window, 'hann80')
		window = my_hann(n, 0.80);

	else
		fail('unknown window %s', window)
	end

elseif length(window) ~= n
	error 'bad window length'
end

window = fftshift(window);

function tf = streq(a,b,n)
%function tf = streq(a, b [,n])
% return 1 if two strings are equal (optionally only up to 1st n chars)
if nargin == 1 && strcmp(a, 'test'), streq_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

if ~ischar(a) | ~ischar(b), tf = 0; return, end
if nargin == 2
	tf = strcmp(a,b);
elseif nargin == 3
	tf = strncmp(a,b,n);
else
	error(mfilename)
end

% my_hann()
function window = my_hann(n, cutoff)
w = round(cutoff * n);
ii = [0:n-1]'-n/2;
window = 0.5 * (1 + cos(2*pi*ii/w)) .* (abs(ii) < w/2);

function window = my_shepp(n, cutoff)
% sh = round(cutoff * n);
% ii = [0:order-1]'-order/2;
% window = (sin(2*pi*ii/(2*sh))./(2*pi*ii/(2*sh))) .* (abs(ii) < sh/2);
order = n;
fi = 2*( 0:(order/2) )./order;
w = 2*pi*(0:size(fi,2)-1)/order;
shp=(sin(w(2:end)/(2*cutoff))./(w(2:end)/(2*cutoff)));
% window = [shp(end-1:-1:2)' ;shp'];
window = [shp(end:-1:1)' ;shp'];

% function window = my_ramp(n, cutoff)
% % sh = round(cutoff * n);
% % ii = [0:order-1]'-order/2;
% % window = (sin(2*pi*ii/(2*sh))./(2*pi*ii/(2*sh))) .* (abs(ii) < sh/2);
% order = n;
% fi = 2*( 0:(order/2) )./order;
% w = 2*pi*(0:size(fi,2)-1)/order;
% window = [w(end:-1:1)' ;w'];

function window = my_cosine(n, cutoff)
order = n;
fi = 2*( 0:(order/2) )./order;
w = 2*pi*(0:size(fi,2)-1)/order;
shp= cos(w(2:end)/(2*cutoff));
% window = [shp(end-1:-1:2)' ;shp'];
window = [shp(end:-1:1)' ;shp'];

% my_hamming()
function window = my_hamming(n, cutoff)
w = round(cutoff * n);
ii = [0:n-1]'-n/2;
window = (0.54 + 0.46 * cos(2*pi*ii/w)) .* (abs(ii) < w/2);
