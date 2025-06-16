[aZ, bZ] = db_construction(6);

% 4 filters comparison
% See 'Compute 4 filters section': https://uk.mathworks.com/help/wavelet/ref/wfilters.html

% Extract the low-pass (scaling) and high-pass (wavelet) decomposition filters
LoD = aZ(:);
HiD = bZ(:);

% Construct reconstruction filters from the decomposition filters
LoR = fliplr(aZ); % Low-pass reconstruction filter (reverse of LoD)
HiR = -fliplr(bZ); % High-pass reconstruction filter (negated reverse of HiD)

% Plot custom decomposition and reconstruction filters
figure(1)
subplot(2,2,1)
stem(LoD)
title("Decomposition Lowpass Filter")
subplot(2,2,2)
stem(HiD)
title("Decomposition Highpass Filter")
subplot(2,2,3)
stem(LoR)
title("Reconstruction Lowpass Filter")
subplot(2,2,4)
stem(HiR)
title("Reconstruction Highpass Filter")

% Retrieve MATLAB's built-in Daubechies 6 filters
[LoD,HiD,LoR,HiR] = wfilters('db6');

% Plot for comparison
figure(2)
subplot(2,2,1)
stem(LoD)
title("Built-In Decomposition Lowpass Filter")
subplot(2,2,2)
stem(HiD)
title("Built-In Decomposition Highpass Filter")
subplot(2,2,3)
stem(LoR)
title("Built-In Reconstruction Lowpass Filter")
subplot(2,2,4)
stem(HiR)
title("Built-In Reconstruction Highpass Filter")