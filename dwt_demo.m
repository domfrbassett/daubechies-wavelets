% Prompt user to select an audio file
choice = questdlg('This code demonstrates how the multilevel DWT can be a used to visualize chirp and frequency break signals. It also neatly demonstrates Heisenberg''s famous indeterminacy principle.', ...
    'DWT Demo', ...
    'Continue', 'Quit', 'Continue');  % Button options

switch choice
    case 'Continue'
    case 'Quit'
        return;  % Exit the script immediately
end
[filename, pathname] = uigetfile({'*.wav;*.mp3;*.flac', 'Audio Files (*.wav, *.mp3, *.flac)'}, 'Select an Audio File For DWT Analysis');
if isequal(filename,0)
    disp('User canceled file selection.');
    return;
end
filePath = fullfile(pathname, filename);

% Read audio file
[x, fs] = audioread(filePath);

% Check for stereo and ask user how to proceed
if size(x,2) > 1
    choice = questdlg('Stereo file detected. How do you want to proceed?', ...
        'Stereo Handling', ...
        'Sum to Mono', 'Use Left Channel', 'Use Right Channel', 'Sum to Mono');
    switch choice
        case 'Sum to Mono'
            x = mean(x, 2); % Average left and right
        case 'Use Left Channel'
            x = x(:,1);     % Left channel
        case 'Use Right Channel'
            x = x(:,2);     % Right channel
        otherwise
            disp('No option selected. Exiting.');
            return;
    end
end

len = length(x);
nLevels = floor(log2(len)); % Automatically decides appropriate number of 
                            % levels (or dyadic decompositions to be precise)
                            % based on signal length.
                            % From Ameera Patel (2015) on an internet
                            % forum: '(For DWT), J0 (the maximum number of 
                            % levels), cannot be greater than j, where N=2^j.' 
                            % Feel free to override as you please.
nbcol = 64; % Number of color levels for coefficient visualisation

% BEANS

[c, l] = discrete_wavelet_transform(x, nLevels);

% Try commenting out the previous line and uncommenting the below 4 lines 
% and running. This shows how discrete_Wavelet_transform.m and MATLAB's 
% wavedec differ in processing my debauchies coefficients.

% [aZ, bZ] = db_construction(6);
% LoD = aZ(:);  % Scaling filter (low-pass)
% HiD = bZ(:);  % Wavelet filter (high-pass)
% [c,l] = wavedec(x,nLevels,LoD, HiD);

% Then comment out all of the above up to the word BEANS and just run the 
% line below. This replaces both my discrete_Wavelet_transform.m and
% db_construction.m scripts and uses the built_in wavelet families with the built in dwt.

% [c,l] = wavedec(x,nLevels,"db6");

% This visualisation part is not my own - it can be found here: https://uk.mathworks.com/help/wavelet/ug/continuous-and-discrete-wavelet-analysis.html

cfd = zeros(nLevels,len);

for k = 1:nLevels
    d = detcoef(c,l,k);
    d = d(:)';
    d = d(ones(1,2^k),:);
    cfd(k,:) = wkeep(d(:)',len);
end

cfd =  cfd(:);
I = find(abs(cfd)<sqrt(eps));
cfd(I) = zeros(size(I));
cfd = reshape(cfd,nLevels,len);
cfd = wcodemat(cfd,nbcol,'row');

h211 = subplot(2,1,1);
h211.XTick = [];
plot(x,'r'); 
title('Analyzed signal');
ax = gca;
ax.XLim = [1 len];
subplot(2,1,2);
colormap(cool(128));
image(cfd);
tics = 1:nLevels; 
labs = int2str(tics');
ax = gca;
ax.YTickLabelMode = 'manual';
ax.YDir = 'reverse';
ax.Box = 'On';
ax.YTick = tics;
ax.YTickLabel = labs;
title('Discrete Transform, absolute coefficients');
ylabel('Level');