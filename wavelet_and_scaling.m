function wavelet_and_scaling(aZ, bZ, numIterations)

% Function to compute and plot the scaling function (phi) and wavelet function (psi)
% using the cascade algorithm

    if nargin < 3
        numIterations = 10; % Set default number of iterations if not passed
    end

    % Normalise filters
    aZ = aZ / norm(aZ);
    bZ = bZ / norm(bZ);

    % Generate scaling function (phi) using cascade
    phi = 1; % Start with an impulse
    for i = 1:numIterations
        phi = conv(phi, aZ); % Convolve repeatedly with low-pass filter
    end

    % Generate the wavelet function (psi) by convolving phi with the high-pass filter
    psi = conv(phi, bZ);

    % Time axes
    t_phi = linspace(0, 1, length(phi));
    t_psi = linspace(0, 1, length(psi));

    % Plot
    figure;
    subplot(2,1,1);
    plot(t_phi, phi, 'b');
    title('Scaling Function (\phi)');
    grid on;

    subplot(2,1,2);
    plot(t_psi, psi, 'r');
    title('Wavelet Function (\psi)');
    grid on;
end