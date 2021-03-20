
% Parameters
signal_length = 1000000;
packet_size = 10;
reps = 5;

%Generating a random signal with binary bits (0 and 1)
tx_bits = randi([0,1],signal_length,1);
disp('First 20 Transmitted Bits');
disp(tx_bits(1:20));

function error = system(tx_bits, reps, SNR, packet_size)

    % ##### Modulation #####
    msg_symbols = bits2bpsk(tx_bits);

    % ##### Repeat #####
    rep_symbols = repetition_coding(msg_symbols,reps);
    
    % ##### Interleaving #####
    interleaved_symbols = interleaver(rep_symbols,packet_size);
    
    % ##### Channel ##### 
    [rx_symbols, h] = channel(interleaved_symbols,SNR);
    
    % ##### De-Interleaving #####
    de_interleaved_symbols = interleaver(rx_symbols,packet_size);
    de_interleaved_h = interleaver(h,packet_size);
    
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

function interleaved_symbols = interleaver(symbols,packet_size)

    buffer_length = mod(numel(symbols),packet_size^2);
    buffer_symbols = symbols(end - buffer_length + 1 : end);
    symbols = symbols(1 : end - buffer_length);
            
    rx_matrix = reshape(symbols, packet_size, packet_size, []);
    
    sz = size(rx_matrix);
    
    buf = zeros(sz(3),packet_size^2);
    
    for i = 1 : sz(3)
       
        buf_1 = rx_matrix(:,:,i);
        buf_1 = buf_1.';
        buf(i,:) = buf_1(:);
        
    end
    
    interleaved_symbols = buf.';
    interleaved_symbols = interleaved_symbols(:);
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
