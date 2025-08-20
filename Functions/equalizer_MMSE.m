function x_hat = equalizer_MMSE(y,H,Es,N0)

% Create MMSE matrix and the corresponding output
W_mmse = (H' * H + (N0/Es) * eye(length(y))) \ H';
x_hat = W_mmse * y;