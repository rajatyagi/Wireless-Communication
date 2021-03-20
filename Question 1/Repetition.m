%% QUESTION 1

% PART (a)
% To obtain the BER vs. SNR curves with L = 1 (no repetition) and 
% L = {2,3,4,5} repetitions over a Rayleigh fading channel.

% Parameters
signal_length = 10000000;

% Generating a random signal with binary bits (0 and 1)
tx_bits = randi([0,1],signal_length,1);
disp('First 20 Transmitted Bits');
disp(tx_bits(1:20));

% L is the number of times we repeat each bit
L = 5;

% Vary SNR from -10 dB to 10 dB
SNR = -10 : 0.1 : 10;
prac_BER = zeros(numel(SNR),1); % array to store practical BER values

for i = 1:numel(SNR)

    prac_BER(i) = system(tx_bits,L,SNR(i));
    disp(SNR(i));
    
end

%% Theoretical vs. Practical BER
% The theoretical bit error rate in BPSK with Rayleigh fading channel is
% given as follows:

%       (2L-1)      1
%  P_e =     C  ----------
%             L  (4SNR)^L

theo_BER = nchoosek(2*L-1,L)./((4*(10.^(SNR/10))).^L)';

figure(1);
semilogy(SNR, prac_BER);
hold on;
semilogy(SNR, theo_BER);
hold off;

%% FUNCTIONS

function error = system(tx_bits,reps,SNR)

    % This function takes the message bits, L and SNR as input.
    % It performs the modulation, repetition coding, reception in fading
    % channel, demodulation, decoding and error calculation (output).
    
    % ##### Modulation #################
    msg_symbols = bits2bpsk(tx_bits);

    % ##### Repeat #####################
    rep_symbols = repetition_coding(msg_symbols,reps);
    
    % ##### Channel ####################
    [rx_symbols, h] = channel(rep_symbols,SNR);
    
    % ##### Maximal Ratio Combiner #####
    rx_vector = maximal_ratio_combiner(rx_symbols, h, reps);
    
    % ##### BPSK Demodulation ##########
    rx_bits = bpsk_demodulation(rx_vector);
    
    % ##### Error ######################
    error = sum(bitxor(rx_bits,tx_bits)) / numel(rx_bits);
    
end

function symbols = bits2bpsk(bits)

    % This function converts bits to BPSK symbols
    % 0 -> -1
    % 1 -> 1
    symbols = bits*2 - 1;

end

function rep_symbols = repetition_coding(symbols,L)

    % Creating a matrix with L number of columns, each columns containing
    % all the message bits.
    Symbol_matrix = repmat(symbols,1,L);  
    
    % Converting the matrix to a linear vector of repeated bits (L times)
    rep_symbols = Symbol_matrix';
    rep_symbols = rep_symbols(:); % Flattening the matrix

end

function [rx_symbols, h] = channel(rep_symbols,SNR)

    length = numel(rep_symbols);
    
    %signal to noise ratio on linear scale for this iteration
    snr_now = 10^(SNR/10); 

    % Channel gain
    
    % real part
    h_real = randn(length, 1);
    
    % imaginary part
    h_imag = randn(length, 1);
    
    h = h_real + 1i * h_imag;
    
    % AWGN Generation
    sigma = sqrt(1/2*snr_now); % standard devation of noise for BPSK
    
    w = sigma*(randn(length,1) + 1i*randn(length,1));
    
    % received signal
    rx_symbols =  h .* rep_symbols + w;

end

function rx_vector = maximal_ratio_combiner(rx_symbols, h, reps)

    length = numel(rx_symbols)/reps;
    
    rx_matrix = reshape(rx_symbols, reps, length);
    
    h_matrix = reshape(h, reps, length);
    
    % Maximal Ratio Combiner
    
    rx_vector = zeros(length,1);
    
    for i = 1:length
        
        h_i = h_matrix(:,i);
        y_i = rx_matrix(:,i);
        
        h_i_norm = norm(h_i);
        
        rx_symbol = real((h_i' * y_i) / h_i_norm);
        
        rx_vector(i) = rx_symbol;
        
    end

end

function rx_bits = bpsk_demodulation(rx_vector)

    rx_bits = (sign(rx_vector) + 1)/2;

end
