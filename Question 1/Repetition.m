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
L = [1 2 3 4 5];

% Vary signal to noise ratio (SNR) from -10 dB to 10 dB
SNR = -10 : 0.5 : 20;

% Array to store practical bit error rate (BER) values
prac_BER = zeros(numel(SNR),numel(L));

% Run the sysmtem for each value of L over all values of SNR
for j = 1:numel(L)
    disp(L(j))
    for i = 1:numel(SNR)

        prac_BER(i, j) = system(tx_bits, L(j), SNR(i));

    end
end

%% PLOTTING THE BER FOR DIFFERENT L

figure(1);
for i = 1:numel(L)
    
    semilogy(SNR, prac_BER(:,i))
    hold on
    title('BER for Different L values');
    xlabel('SNR (dB)');
    ylabel('BER');
    text(SNR(35),prac_BER(35, i),"L="+num2str(L(i)))
    
end
hold off
grid on;

%% FUNCTIONS

function error = system(tx_bits, reps, SNR)

    % This function takes the message bits, L and SNR as input.
    % It performs the modulation, repetition coding, reception in fading
    % channel, demodulation, decoding and error calculation (output).
    
    % ##### MODULATION #################
    msg_symbols = bits2bpsk(tx_bits);

    % ##### REPETITION CODING ##########
    rep_symbols = repetition_coding(msg_symbols, reps);
    
    % ##### CHANNEL ####################
    [rx_symbols, h] = channel(rep_symbols, SNR);
    
    % ##### MAXIMAL RATIO COMBINER #####
    rx_vector = maximal_ratio_combiner(rx_symbols, h, reps);
    
    % ##### BPSK DEMODULATION ##########
    rx_bits = bpsk_demodulation(rx_vector);
    
    % ##### ERROR CALCULATION ##########
    error = sum(bitxor(rx_bits, tx_bits)) / numel(rx_bits);
    
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
    rep_symbols = Symbol_matrix.'; % transpose
    rep_symbols = rep_symbols(:); % flattening the matrix

end

function [rx_symbols, h] = channel(rep_symbols, SNR)

    % This function emulates a Rayleigh fading channel with
    % Channel gain => h ~ CN(0,1)
    % AWGN => w ~ CN(0, sigma^2)
    % where, sigma = sqrt(1/2*SNR)
    
    length = numel(rep_symbols);
    
    % Signal to noise ratio on linear scale for the iteration
    snr_now = 10^(SNR/10); 

    % CHANNEL GAIN
    % Rayleigh fading channel gain is given by a 
    % complex normal distribution.
    
    h_real = sqrt(0.5)*randn(length, 1); % real part
    h_imag = sqrt(0.5)*randn(length, 1); % imaginary part
    
    h = h_real + 1i*h_imag; % complex number
    
    % ADDITIVE WHITE GAUSSIAN NOISE
    sigma = sqrt(1/(2*snr_now)); % standard devation of noise for BPSK
    
    w = sigma*(randn(length,1) + 1i*randn(length,1)); % complex noise
    
    % Received signal
    rx_symbols =  h .* rep_symbols + w;

end

function rx_vector = maximal_ratio_combiner(rx_symbols, h, reps)

    % This function performs the maximal ratio combining:
    % Weigh the received signal in each branch in proportion to the signal
    % strength and also align the phases of the signals in the summation to
    % maximize the output SNR. (coherent combining)
    
    %   h*                 h*
    % ----- y = ||h||x + ----- w
    % ||h||              ||h||
    
    % h => L length vector of channel gains (h* => Hermitian of h)
    % y => L length vector of repeated bits
    % w => L length vector of AWGN
    
    length = numel(rx_symbols)/reps;
    
    % Reshaping the received symbol sequence and the channel gain sequence
    % according to repetitions encoded.
    rx_matrix = reshape(rx_symbols, reps, length);
    h_matrix = reshape(h, reps, length);
    
    % MAXIMAL RATIO COMBINER (Matched Filter)
    
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

    % This function demodulates the received symbols to BPSK symbols by
    % using Maximum Likelihood (ML) detection technique.
    rx_bits = (sign(rx_vector) + 1)/2;

end
