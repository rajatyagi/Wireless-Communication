% Parameters
signal_length = 100000;
SNR = -10 : 0.5 : 20;

%Generating a random signal with binary bits (0 and 1)
tx_bits = randi([0,1],signal_length,1);

BER = zeros(numel(SNR),1);

for i = 1:numel(SNR)

    BER(i) = system(tx_bits, SNR(i));

end

%% PLOTTING THE BER FOR ALL CASES

figure(1);
semilogy(SNR, BER);
title('BER for Different Cases');
xlabel('SNR (dB)');
ylabel('BER');
grid on;

%% FUNCTIONS

function error = system(tx_bits, SNR)

    % ##### Modulation #####
    msg_symbols = bits2qpsk(tx_bits);

    % ##### Alamouti Encoding #####
    [tx_1, tx_2] = alamouti_encoder(msg_symbols);
    
    % ##### Channel ##### 
    [rx_symbols, h1, h2] = channel(tx_1, tx_2, SNR);
    
    % ##### Alamouti Decoding #####
    decoded_symbols = alamouti_decoder(rx_symbols, h1, h2);
    
    % ##### QPSK Demodulation #####
    detected_symbols = qpsk_detector(decoded_symbols);
    rx_bits = qpsk_demod(detected_symbols);
    
    % ##### Error #####
    error = sum(bitxor(rx_bits,tx_bits)) / numel(rx_bits);
    
end

function symbols = bits2qpsk(bits)

    symbols = zeros(length(bits)/2, 1);
    %With Gray Labelling    
    for i = 1 : 2 : length(bits)
        % 00 --> 1+j
        if bits(i) == 0 && bits(i+1) == 0
            symbols((i + 1)/2) = 0.5 + 0.5*1i;
        % 01 --> -1+j
        elseif bits(i) == 0 && bits(i+1) == 1
            symbols((i + 1)/2) = -0.5 + 0.5*1i;
        % 10 --> -1-j    
        elseif bits(i) == 1 && bits(i+1) == 0
            symbols((i + 1)/2) = 0.5 - 0.5*1i;
        % 11 --> 1-j    
        elseif bits(i) == 1 && bits(i+1) == 1
            symbols((i + 1)/2) = -0.5 - 0.5*1i; 
        end

        %          With Gray Labelling
        %                  |
        %          01 *    |    * 00
        %                  |
        %         _________|_________
        %                  |
        %                  |
        %          11 *    |    * 10
        %                  |
    end
    
end

function [tx_1, tx_2] = alamouti_encoder(symbols)

    length = numel(symbols);
    
    tx_1 = symbols;
    tx_2 = zeros(length,1);
        
    for i = 1 : length
    
        if(mod(i,2) == 0)
            tx_1(i) = -1 * symbols(i)';
            tx_2(i) = symbols(i - 1)';
        else
            tx_2(i) = symbols(i+1);
        end
        
    end
    
end

function [rx_symbols, h1, h2] = channel(tx_1, tx_2, SNR)

    length = numel(tx_1);
    
    %signal to noise ratio on linear scale for this iteration
    snr_now = 10^(SNR/10); 
    
    % Channel gain
    coherence_frames = ceil(length/2);
    
    h1_real = sqrt(0.5)*randn(coherence_frames, 1); % real part
    h1_imag = sqrt(0.5)*randn(coherence_frames, 1); % imaginary part
    
    h2_real = sqrt(0.5)*randn(coherence_frames, 1); % real part
    h2_imag = sqrt(0.5)*randn(coherence_frames, 1); % imaginary part
    
    h1 = h1_real + 1i*h1_imag;
    h2 = h2_real + 1i*h2_imag;
    
    h1 = repetition_coding(h1,2);
    h2 = repetition_coding(h2,2);
    
    % AWGN Generation
    sigma = sqrt(1/(2*snr_now)); % standard devation of noise for BPSK
        
    w = sigma*(randn(length,1) + 1i*randn(length,1));
    
    % received signal
    rx_symbols =  h1 .* tx_1 +  h2 .* tx_2 + w;

end

function rep_symbols = repetition_coding(symbols,L)

    % creating a matrix with reps number of columns and each containing the message bits
    Symbol_matrix = repmat(symbols,1,L);  
    
    %converting the matrix to a linear vector of repeated bits
    rep_symbols = Symbol_matrix';
    rep_symbols = rep_symbols(:);

end

function decoded_symbols = alamouti_decoder(symbols, h1, h2)

    length = numel(symbols);
    decoded_symbols = zeros(length,1);
    
    for i = 1 : 2 : length
    
        y = [symbols(i) symbols(i+1)'].';
        h_orth_1 = [h1(i) h2(i)'].';
        h_orth_2 = [h2(i) (-1 * h1(i)')].';
        
        decoded_symbols(i) = ((h_orth_1') * y) / norm(h_orth_1);
        decoded_symbols(i+1) = ((h_orth_2') * y) / norm(h_orth_2);
        
    end

end

function detected_symbol = qpsk_detector(symbol) 

    %Phase of the current symbol in degrees
    phase = angle(symbol)*180/pi;
    length = numel(phase);
    detected_symbol = zeros(length,1);
    
    %Symbol Detection using ML Rule:
    for i = 1:numel(phase)
        %if phase --> [0,90] symbol = 1+j     I Quadrant
        if phase(i) >= 0 && phase(i) < 90
            detected_symbol(i) = 0.5+0.5*1i;

        %if phase --> [90,180] symbol = -1+j  II Quadrant
        elseif phase(i) >= 90 && phase(i) < 180
            detected_symbol(i) = -0.5+0.5*1i;

        %if phase --> [0,-90] symbol = 1-j    IV Quadrant
        elseif phase(i) < 0 && phase(i) >= -90
            detected_symbol(i) = 0.5-0.5*1i;

        %if phase --> [-90,-180] symbol = -1-j  III Quadrant  
        elseif phase(i) < -90 && phase(i) >= -180
            detected_symbol(i) = -0.5-0.5*1i;
        end

    end
        
end

function recvd_bits = qpsk_demod(recvd_symbols)

    recvd_bits = zeros(2*length(recvd_symbols),1);
        
    %bits assignment for every symbol
    for i = 1 : 2 : length(recvd_bits)

        symbol = recvd_symbols((i+1)/2);

        % 1+j --> 00
        if symbol == 0.5 + 0.5*1i
            recvd_bits(i) = 0;
            recvd_bits(i+1) = 0;               
        % -1+j --> 01
        elseif symbol == -0.5 + 0.5*1i
            recvd_bits(i) = 0;
            recvd_bits(i+1) = 1;            
        % -1-j --> 11    
        elseif symbol == -0.5 - 0.5*1i
            recvd_bits(i) = 1;
            recvd_bits(i+1) = 1;
        % 1-j --> 10    
        elseif symbol == 0.5 - 0.5*1i
            recvd_bits(i) = 1;
            recvd_bits(i+1) = 0;
        end

    end

    %          With Gray Labelling
    %                  |
    %          01 *    |    * 00
    %                  |
    %         _________|_________
    %                  |
    %                  |
    %          11 *    |    * 10
    %                  |

end










