function [BER,SER,FER] = sim_fun_ODDM_v3(new_frames,parameters)

% Make parameters
fields = fieldnames(parameters);
for i = 1:numel(fields)
    eval([fields{i} ' = parameters.(fields{i});']);
end

% Cause error if settings are inconsistent with theory
Ts = T / M;
if 2*Q + floor(2510*10^(-9) / Ts) >= M
    error("Settings can not satisfy the ambiguity assumption.")
end

% Define parameters
res = 10;
Es = 1;
Eb = Es / log2(M_ary);
N_iters = 10;
syms_per_f = M*N;
N0 = Eb / (10^(EbN0 / 10)) * ((N+2)/N);

% Add redundancy for rectangular and sinc pulses
if shape == "rect"
    Q = 1;
    alpha = 1;
elseif shape == "sinc"
    alpha = 1;
end

% Data setup
if M_ary == 2
    % Define bit order
    bit_order = [0;1];

    % Define alphabet
    alphabet_set = linspace(1,M_ary,M_ary)';
    S = sqrt(Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ M_ary);
elseif M_ary == 4
    % Define bit order
    bit_order = [0,0;0,1;1,0;1,1];

    % Define alphabet
    S = zeros(4,1);
    S(1) = (sqrt(2)/2) + (1j*sqrt(2)/2);
    S(2) = (sqrt(2)/2) - (1j*sqrt(2)/2);
    S(3) = -(sqrt(2)/2) + (1j*sqrt(2)/2);
    S(4) = -(sqrt(2)/2) - (1j*sqrt(2)/2);
    S = sqrt(Es) .* S;
end

% Render ambiguity table
[Ambig_Table.vals,Ambig_Table.t_range,Ambig_Table.f_range] = gen_DD_cross_ambig_table(N,M,T,Fc,vel,shape,alpha,Q,res);

% Reset bit errors for each SNR
bit_errors = zeros(new_frames,syms_per_f*log2(M_ary));
sym_errors = zeros(new_frames,syms_per_f);
frm_errors = zeros(new_frames,1);

for frame = 1:new_frames

    % Generate data
    [TX_bit,TX_sym,xDD] = gen_data(bit_order,S,syms_per_f);

    % Generate H matrix and channel information
    [HDD,L1,L2] = gen_HDD_direct(T,N,M,Fc,vel,Q,Ambig_Table);

    % Generate noise
    zDD = sqrt(N0/2) * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));

    % Construct received signal - Discrete
    yDD = HDD * xDD + zDD;

    % Equalize received signal
    x_hat = OTFS_pulse_equalizer_AWGN(yDD,HDD,N,M,L2,-L1,Es,N0,S,N_iters);

    % Hard detection for final x_hat
    dist = abs(x_hat.' - S).^2;
    [~,min_index] = min(dist);
    RX_sym = min_index.' - 1;

    % Convert final RX_sym to RX_bit
    RX_bit = bit_order(RX_sym+1,:);

    % Error calculation
    diff_bit = TX_bit ~= RX_bit;
    diff_sym = TX_sym ~= RX_sym;
    bit_errors(frame,:) = diff_bit(:).';
    sym_errors(frame,:) = diff_sym(:).';
    if sum(sym_errors(frame,:)) > 0
        frm_errors(frame) = 1;
    end

end

% Calculate BER, SER and FER
BER = sum(bit_errors,"all") / (new_frames*syms_per_f*log2(M_ary));
SER = sum(sym_errors,"all") / (new_frames*syms_per_f);
FER = sum(frm_errors,"all") / (new_frames);