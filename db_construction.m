function [aZ, bZ, pZ, PA] = db_construction(A)

% Daubechies Wavelet Filter Construction
% 
% Define Minimal Polynomial P_A(X):
%    - This is a specific polynomial used in the Daubechies construction,
%      defined as:
%       - P_A(X) = sum_{k=0}^{A-1} binomial(A+k-1, A-1) * 2^(-k) * X^k
%      and is chosen to guarantee A vanishing moments and orthonormality.
%      It satisfies the identity:
%       - (1-X)^A * P_A(X) + X^A * P_A(1-X) = 1
%
% Spectral Factorization:
%    - Find p(Z) from P_A(X) roots using quadratic:
%      Z^2 - 2*(1-mu)*Z + 1 = 0
%    - Choose roots inside the unit circle to form p(Z).
%
% Construct Final Filter a(Z):
%    - a(Z) = 2^(1-A) * (1+Z)^A * p(Z)
%    - Normalize a(Z), compute high-pass b(Z) using QMF:
%      b(Z) = (-1)^n * reversed a(Z)
%
% References:
% Daubechies (1988), "Orthonormal bases of compactly supported wavelets."
% https://en.wikipedia.org/wiki/Daubechies_wavelet
% 3blue1brown was also very informative

% Outputs:
% aZ = LPF coefficients (scaling function)
% bZ = HPF coefficients (wavelet function)
% pZ = Spectral factor polynomial (numeric)
% PA = Minimal polynomial (symbolic)

    % Define Minimal Polynomial P_A(X)
    % P_A(X) = sum_{k=0}^{A-1} [ (A+k-1 choose A-1) * 2^{-k} * X^k ]
    % Satisfies: (1-X)^A P_A(X) + X^A P_A(1-X) = 1
    syms X;
    PA = sym(0);
    for k = 0:A-1
        coeff = nchoosek(A + k - 1, A - 1) * 2^(-k); % Binomial coefficient
        PA = PA + coeff * X^k; % Build P_A(X) term-by-term
    end

    fprintf('Minimal polynomial P_A(X):\n');
    disp(simplify(PA));

    % Spectral Factorization to Extract p(Z)
    % Factor P_A(X) into p(Z) with roots inside the unit circle (minimum phase)
    %   - Find roots mu_i of P_A(X).
    %   - For each mu_i, solve Z^2 - 2(1-mu_i)Z + 1 = 0
    %   - Choose root |Z| < 1 to ensure stability
    coeffsP = double(fliplr(coeffs(PA, 'All'))); % Numeric coefficients (highest power first)
    rootsP = roots(coeffsP); % Roots mu_i of P_A(X)

    pZroots = [];
    for i = 1:length(rootsP)
        mu = rootsP(i);
        % Quadratic equation: Z^2 - 2(1-mu)Z + 1 = 0 (from X(Z) = 1 - 0.5(Z + Z^{-1}))
        Zsol = roots([1, -2*(1 - mu), 1]);
        % Choose root inside or on unit circle (minimum phase)
        if abs(Zsol(1)) <= 1
            pZroots(end+1) = Zsol(1);
        else
            pZroots(end+1) = Zsol(2);
        end
    end

    pZ = poly(pZroots); % Convert roots to polynomial
    pZ = real(pZ); % Enforce real coefficients
    fprintf('Approximated p(Z) coefficients:\n');
    disp(pZ);

    % Construct Lowpass Filter a(Z)
    % a(Z) = 2^{1-A} (1+Z)^A p(Z)
    %   - (1+Z)^A enforces A vanishing moments
    %   - p(Z) ensures orthonormality and minimum phase
    syms Z;
    onePlusZ = (1 + Z)^A; % Vanishing moments term
    aSym = 2^(1 - A) * onePlusZ * poly2sym(pZ, Z); % Combine terms symbolically
    aZ = double(fliplr(coeffs(aSym, 'All'))); % Convert to numeric coefficients
    aZ = aZ / norm(aZ); % Normalise energy (||aZ||_2 = 1)

    fprintf('Approximate a(Z) coefficients (normalized):\n');
    disp(aZ);

    % Construct Highpass Filter b(Z) via Quadrature Mirror Filter (QMF) relation
    % b(Z) = (-1)^k * a_{L-1-k}
    bZ = (-1).^(0:length(aZ)-1) .* fliplr(aZ);
    fprintf('Detail filter b(Z) coefficients (QMF):\n');
    disp(bZ);

% Now you can: 
%   - Visualise the wavelet and scaling functions from 
%     runplot_wavelet_and_scaling.m
%   - Compare the resulting decomposition and reconstruction filters with
%     MATLAB's preprogrammed debauchies wavelets from
%     runplot_four_filters.m - db6 exhibits a relatively high level of 
%     correlation (I have not implemented reconstruction as it is outside 
%     the scope of the brief, however I am in the process of building a 
%     DWT-SPIHT algorithm for scalable-to-lossless audio compression with 
%     bandwidth constraints which I should be able to show you next year!)
%   - Run dwt_demo.m to see DWT in action (be aware, this calls
%     discrete_wavelet_transform.m, where the approximation order is set to
%     6 by default - this can be changed but there's no UI)