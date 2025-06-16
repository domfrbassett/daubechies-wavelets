[aZ, bZ, pZ, PA] = db_construction(6);

% Display the low-pass analysis filter coefficients
disp('aZ (low-pass analysis filter):')
disp(aZ)

% Display the high-pass analysis filter coefficients
disp('bZ (high-pass analysis filter):')
disp(bZ)

% Display the spectral factor polynomial used in filter construction
disp('pZ (spectral factor polynomial):')
disp(pZ)

% Display the minimal polynomial used during spectral factorization
disp('PA (minimal polynomial):')
disp(PA)

% Compute and plot the scaling function (phi) and wavelet function (psi)
wavelet_and_scaling(aZ, bZ, 10); % Ten iterations