clear all, clc, close all;
[smagNote, smagMusic, sphaseMusic, Fs_1] = load_data();

%% Solution for Problem 2_1 here
% Store W in a file called "problem2_1.mat"
% W will be a 15xT matrix, where T is the number of frames in the music.
W = pinv(smagNote)*smagMusic;
for i=1:numel(W)
    if W(i)<0
        W(i) = 0;
    end
end
save('results/problem2_1.mat','W');

%% Synthesize Music
% Use the 'synthesize_music' function here.
% Use 'wavwrite' function to write the synthesized music as 'problem2_1_synthesis.wav' to the 'results' folder.
smagMusicProj = smagNote*W;
synMusic = synthesize_music2048(sphaseMusic, smagMusicProj);
audiowrite('results/problem2_1_synthesis.wav', synMusic, Fs_1);


%% Solution for Problem 2_2 here

% Find and store the transformation matrix
[smagAL, smagBL, smagCL, sphaseCL, Fs_C] = load_data2();

X =  smagBL * pinv(smagAL);

save('results/problem2_2.mat','X');

smagMusicProj2 = X * smagCL;
synMusic2 = synthesize_music1024(sphaseCL, smagMusicProj2);

audiowrite('results/problem2_2_synthesis.wav', synMusic2, Fs_C);

% Apply the transformation matrix to audio C and store the created music using 'synthesize_music' function.
% Use 'wavwrite' function to write the synthesized music as 'problem2_2_synthesis.wav' to the 'results' folder.

function [smagNote, smagMusic, sphaseMusic, Fs_1] = load_data()

materialfolder = 'materials';
notesfolder = 'notes15';

[poly, Fs] = audioread('materials/polyushka.wav');

poly_spectrum = stft(poly', 2048, 256, 0, hann(2048));
poly_stft = abs(poly_spectrum);
poly_phase = poly_spectrum./(poly_stft+eps);

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

smagNote = notes;
smagMusic = poly_stft;
sphaseMusic = poly_phase;
Fs_1 = Fs;
end
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
function [synMusic] = synthesize_music2048(sphaseMusic,smagMusicProj)
%% Argument Descriptions
% Required Input Arguments:
% sphaseMusic: 1025 x K matrix containing the spectrum phases of the music after STFT.
% smagMusicProj: 1025 x K matrix, reconstructed version of smagMusic using transMatT

synMusicSpectrum = smagMusicProj.*sphaseMusic;
synMusic = stft(synMusicSpectrum, 2048, 256, 0, hann(2048));
synMusic = synMusic';

% Required Output Arguments:
% synMusic: N x 1 music signal reconstructed using STFT.


%% Music synthesis
% Fill your code here to return 'synMusic'
end
function [smagAL, smagBL, smagCL, sphaseCL, Fs_C] = load_data2()

[audioA, fsA] = audioread('materials\Audio\silentnight_piano.aif');
[audioB, fsB] = audioread('materials\Audio\silentnight_guitar.aif');
[audioC, fsC] = audioread('materials\Audio\littlestar_piano.aif');

audioAL = audioA(:,1);
spectrumA = stft(audioAL', 1024, 256, 0, hann(1024));
music_stftA = abs(spectrumA);

audioBL = audioB(:,1);
spectrumB = stft(audioBL', 1024, 256, 0, hann(1024));
music_stftB = abs(spectrumB);

audioCL = audioC(:,1);
spectrumC = stft(audioCL', 1024, 256, 0, hann(1024));
music_stftC = abs(spectrumC);
sphaseC = spectrumC ./(abs(spectrumC)+eps);

for i=1:numel(music_stftC)
    if music_stftC(i)<0
        music_stftC(i) = 0;
    end
end

smagAL = music_stftA;
smagBL = music_stftB;
smagCL = music_stftC;
sphaseCL = sphaseC;
Fs_C = fsC;

end
function [synMusic] = synthesize_music1024(sphaseMusic,smagMusicProj)
%% Argument Descriptions
% Required Input Arguments:
% sphaseMusic: 1025 x K matrix containing the spectrum phases of the music after STFT.
% smagMusicProj: 1025 x K matrix, reconstructed version of smagMusic using transMatT

synMusicSpectrum = smagMusicProj.*sphaseMusic;
synMusic = stft(synMusicSpectrum, 1024, 256, 0, hann(1024));
synMusic = synMusic';

% Required Output Arguments:
% synMusic: N x 1 music signal reconstructed using STFT.


%% Music synthesis
% Fill your code here to return 'synMusic'
end