function [c, l] = discrete_wavelet_transform(x, nLevels)

% Performs multilevel DWT

% Inputs: x - Input signal as a row or column vector
%         nLevels - Number of decomposition levels
% Outputs: c - concatenated wavelet coefficients
%          l - vector of lengths

    % Get Daubechies-6 low-pass and high-pass decomposition filters
    [aZ, bZ] = db_construction(6);

    LoD = aZ(:); % LPF (column vector)
    HiD = bZ(:); % HPF (column vector)

    % Ensure input is a row vector
    x = x(:).';

    A = x; % Approximation coefficients at current level
    cD = cell(1, nLevels); % To store detail coefficients at each level

    for level = 1:nLevels
        % Symmetric padding
        lf = length(LoD); % Filter length
        padLen = floor(lf-1); % Padding length on each side
        A_padded = padarray(A, [0, padLen], 'symmetric', 'both');

        % Convolve with filters and downsample
        A_conv = conv(A_padded, LoD, 'valid'); % Approximation convolution
        D_conv = conv(A_padded, HiD, 'valid'); % Detail convolution

        % Downsample by 2
        A_next = A_conv(1:2:end);
        D_next = D_conv(1:2:end);

        % Store detail coefficients at the current level
        % Stored in reverse order (finest to coarsest)
        cD{nLevels - level + 1} = D_next;

        % Prepare approximation for the next level
        A = A_next;
    end

    % Combine coefficients into a single vector (coarse-to-fine)
    c = [A, cell2mat(cD)];

    % Track the lengths of all coefficient vectors (including original signal length)
    l = [length(A), cellfun(@length, cD), length(x)];
end