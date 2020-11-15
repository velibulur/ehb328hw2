clear, clc, close all;

% Firstly read the audio.
materialfolder = 'hw1materials';
soundname = dir([materialfolder filesep 'polyushka.wav']);
[poly, Fs] = audioread([materialfolder filesep soundname.name]);

% Figure the original spectrum of audio.
figure
spectrogram(poly,'yaxis')
title('Original Sound Spectrum')

notesfolder = 'notes15';
portion = 1:Fs*5;

poly_spectrum = stft(poly', 2048, 256, 0, hann(2048));
poly_stft = abs(poly_spectrum);
poly_phase = poly_spectrum./(poly_stft+eps);

% Adding the 15 notes in to 'notes' array.
listname = dir([materialfolder filesep notesfolder filesep '*.wav']);
notes = [];
for i=1:length(listname)
  [s, Fs] = audioread([materialfolder filesep notesfolder filesep listname(i).name]);
  s = s(:,1);
  s = resample(s, 16000, Fs);
  spectrum = stft(s', 2048, 256, 0, hann(2048));
  middle = ceil(size(spectrum, 2) /2); 
  note = abs(spectrum(:, middle)); 
  note(find(note<max(note(:))/100)) = 0 ;
  note = note/norm(note);
  notes = [notes, note];
end

% To find the W, pseudoinverse of notes * stft of the audio poly.
W = pinv(notes)*poly_stft;

% W's negative values go to zero.
for i=1:numel(W)
    if W(i)<0
        W(i) = 0;
    end
end

polysynth_stft = notes*W;

% Figure of synthesised power spectrum of audio's stft.
figure
imagesc(polysynth_stft)
title('Synthesised Power Spectrum2')

% For figure of synthesised sound spectrum, invert the stft of audio.
polysynth_spectrum = polysynth_stft.*poly_phase;
polysynth = stft(polysynth_spectrum, 2048, 256, 0, hann(2048));
polysynth = polysynth';

% Figure the synthesised sound spectrum.
figure
spectrogram(polysynth,'yaxis')
title('Synthesised Sound Spectrum')

% Obtain the synthesised audio in format of wav.
audiowrite('poly_synth.wav', polysynth, Fs);

% Matlab file in the attached zip combines the two codes. Please inspect that one as well.
% These codes can be found at https://github.com/velibulur/ehb328hw2 and https://github.com/meserbetcioglu/ehb328hw2
function [f,fp] = stft( x, sz, hp, pd, w)
% [f,fp] = stft( x, sz, hp, pd, w)
%x = signal
%sz = fft size
%hp = hopsize between adajcent frames (in points)
%pd = 0 padding (in points)
%w = window (optional; default is boxcar)
%Returns:
%f = stft (complex)
%fp = phase
%
%To reconstruct, x must be a complex array (i.e. an stft)
%                rest stays the same
%
% This code traces its ownership to several people from Media labs, MIT
%


% Forward transform
if isreal( x)

	% Defaults
	if nargin < 5
		w = 1;
	end
	if nargin < 4
		pd = 0;
	end
	if nargin < 3
		hp = sz/2;
	end

	% Zero pad input
%	x = [x zeros( 1, ceil( length(x)/sz)*sz-length(x))];
        extra = (length(x)-sz)/hp;
        padding = ceil(extra)*hp + sz - length(x);
	x = [x zeros( 1, padding)];
%	x = [zeros( 1, sz+pd) x zeros( 1, sz+pd)];

	% Pack frames into matrix
	s = zeros( sz, (length(x)-sz)/hp);
	j = 1;
	for i = sz:hp:length( x)
		s(:,j) = w .* x((i-sz+1):i).';
		j = j + 1;
	end

	% FFT it
	f = fft( s, sz+pd);

	% Chop redundant part
	f = f(1:end/2+1,:);
	
	% Return phase component if asked to
	if nargout == 2
		fp = angle( f);
		fp = cos( fp) + sqrt(-1)*sin( fp);
	end

% Inverse transform
else

	% Defaults
	if nargin < 5
		w = 1;
	end
	if nargin < 4
		pd = 0;
	end
	if nargin < 3
		hp = sz/2;
	end

	% Ignore padded part
	if length( w) == sz
		w = [w; zeros( pd, 1)];
	end

	% Overlap add/window/replace conjugate part
	f = zeros( 1, (size(x,2)-1)*hp+sz+pd);
	v = 1:sz+pd;
	for i = 1:size( x,2)
		f((i-1)*hp+v) = f((i-1)*hp+v) + ...
			(w .* real( ifft( [x(:,i); conj( x(end-1:-1:2,i))])))';
	end

	% Norm for overlap
	f = f / (sz/hp);
	f = f(sz+pd+1:end-sz-2*pd);
end
end