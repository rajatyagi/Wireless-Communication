%% QUESTION 1

% PART (b)
% Let us say that each coherence period corresponds to 8 BPSK samples.
% Create your own random channel which becomes “bad” once in every four
% coherence time periods (you can define “bad” based on intuition).
% Show through BER vs. SNR simulations that "interleaving + repetition"
% has better performance than "only interleaving" or "only repetition".

% Parameters
signal_length = 10000000;
SNR = -10 : 0.5 : 20;

%Generating a random signal with binary bits (0 and 1)
tx_bits = randi([0,1],signal_length,1);
disp('First 20 Transmitted Bits');
disp(tx_bits(1:20));

% Only Interleaving
packet_size = 5;
reps = 1;

interleave_BER = zeros(numel(SNR), 1);

for i = 1:numel(SNR)

    interleave_BER(i) = system(tx_bits, reps, SNR(i), packet_size);

end

% Only Repetition Coding
packet_size = 1;
reps = 5;

repeat_BER = zeros(numel(SNR), 1);

for i = 1:numel(SNR)

    repeat_BER(i) = system(tx_bits, reps, SNR(i), packet_size);

end

% Interleaving + Repetition Coding
packet_size = 5;
reps = 5;

inter_rep_BER = zeros(numel(SNR), 1);

for i = 1:numel(SNR)

    inter_rep_BER(i) = system(tx_bits, reps, SNR(i), packet_size);

end

%% PLOTTING THE BER FOR ALL CASES

figure(1);
semilogy(SNR, interleave_BER)
hold on
semilogy(SNR, repeat_BER)
semilogy(SNR, inter_rep_BER)
title('BER for Different Cases');
xlabel('SNR (dB)');
ylabel('BER');
legend('Interleaving', 'Repetition Coding', 'Interleaving + Repetition Coding');
legend
hold off
grid on;

%% FUNCTIONS

function error = system(tx_bits, reps, SNR, packet_size)

    % This function takes the message bits, repetitions, SNR and packet
    % size as input. It performs the below defined tasks to give error.

    % ##### Modulation #####
    msg_symbols = bits2bpsk(tx_bits);

    % ##### Repeat #####
    rep_symbols = repetition_coding(msg_symbols,reps);
    
    % ##### Interleaving #####
    if(packet_size ~= 1)
        interleaved_symbols = interleaver(rep_symbols,packet_size);
    else
        interleaved_symbols = rep_symbols;
    end
    
    % ##### Channel ##### 
    [rx_symbols, h] = channel(interleaved_symbols,SNR);
    
    % ##### De-Interleaving #####
    if(packet_size ~= 1)
        de_interleaved_symbols = interleaver(rx_symbols,packet_size);
        de_interleaved_h = interleaver(h,packet_size);
    else
        de_interleaved_symbols = rx_symbols;
        de_interleaved_h = h;
    end
    
    % ##### Maximal Ratio Combiner #####
    rx_vector = maximal_ratio_combiner(de_interleaved_symbols, de_interleaved_h, reps);
        
    % ##### BPSK Demodulation #####
    rx_bits = bpsk_demodulation(rx_vector);
    
    % ##### Error #####
    error = sum(bitxor(rx_bits,tx_bits)) / numel(rx_bits);
    
end

function symbols = bits2bpsk(bits)

    symbols = bits*2 - 1;

end

function rep_symbols = repetition_coding(symbols,L)

    % creating a matrix with reps number of columns and each containing the message bits
    Symbol_matrix = repmat(symbols,1,L);  
    
    %converting the matrix to a linear vector of repeated bits
    rep_symbols = Symbol_matrix';
    rep_symbols = rep_symbols(:);

end

function interleaved_symbols = interleaver(symbols, packet_size)

    % dividing the symbol stream into blocks of packet_size
    buffer_length = mod(numel(symbols),packet_size^2);
    
    % removing the buffer symbols (not a multiple of packet_size)
    buffer_symbols = symbols(end - buffer_length + 1 : end);
    symbols = symbols(1 : end - buffer_length);
            
    rx_matrix = reshape(symbols, packet_size, packet_size, []);
    
    sz = size(rx_matrix);
    
    buf = zeros(sz(3),packet_size^2);
    
    % Interleaving all the symbols of a packet with consecutive
    % (packet_size - 1) packets.
    for i = 1 : sz(3)
       
        buf_1 = rx_matrix(:,:,i);
        buf_1 = buf_1.';
        buf(i,:) = buf_1(:);
        
    end
    
    interleaved_symbols = buf.';
    interleaved_symbols = interleaved_symbols(:);
    
    % Adding the removed buffer symbols
    interleaved_symbols = vertcat(interleaved_symbols, buffer_symbols);

end

function [rx_symbols, h] = channel(rep_symbols,SNR)

    length = numel(rep_symbols);
    
    %signal to noise ratio on linear scale for this iteration
    snr_now = 10^(SNR/10); 
    
    % Channel gain
    coherence_frames = ceil(length/8);
    
    h_real = sqrt(0.5)*randn(coherence_frames, 1); % real part
    h_imag = sqrt(0.5)*randn(coherence_frames, 1); % imaginary part
    
    h = h_real + 1i*h_imag;
    
    for i = 1 : numel(h)
       
        if(mod(i,4) == 0)            
            h(i) = h(i)*(1/(10*snr_now));            
        end
                
    end
    
    h = repetition_coding(h,8);
    h = h(1:length);
    
    % AWGN Generation
    sigma = sqrt(1/(2*snr_now)); % standard devation of noise for BPSK
        
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
