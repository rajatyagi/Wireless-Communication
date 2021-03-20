% Parameters
signal_length = 10000000;
SNR = -10 : 0.5 : 25;

%Generating a random signal with binary bits (0 and 1)
tx_bits = randi([0,1],signal_length,1);

alamouti_BER = zeros(numel(SNR),1);
repetition_BER = zeros(numel(SNR),1);
theoretical_BER = zeros(numel(SNR),1);

for i = 1:numel(SNR)

    snr_now = 10^(SNR(i)/10); 
    alamouti_BER(i) = alamouti_system(tx_bits, SNR(i));
    repetition_BER(i) = repetition_system(tx_bits, SNR(i));
    theoretical_BER(i) = 1/(snr_now^2);
    
end
%% PLOTTING THE BER FOR ALL CASES

figure(1);
semilogy(SNR, alamouti_BER);
hold on;
semilogy(SNR, repetition_BER);
semilogy(SNR, theoretical_BER);
hold off;
title('BER vs. SNR for Different Cases');
legend('Alamouti Coding','Repetition Coding', 'Theoretical Upper Bound');
xlabel('SNR (dB)');
ylabel('BER');
grid on;

%% FUNCTIONS

function error = alamouti_system(tx_bits, SNR)

    % Bit Energy
    Eb = 0.25;

    % ##### Modulation #####
    msg_symbols = bits2qpsk(tx_bits);

    % ##### Alamouti Encoding #####
    [tx_1, tx_2] = alamouti_encoder(msg_symbols);
    
    % ##### Channel ##### 
    [rx_symbols, h1, h2] = channel(tx_1, tx_2, SNR, Eb);
    
    % ##### Alamouti Decoding #####
    decoded_symbols = alamouti_decoder(rx_symbols, h1, h2);
    
    % ##### QPSK Demodulation #####
    detected_symbols = qpsk_detector(decoded_symbols);
    rx_bits = qpsk_demod(detected_symbols);
    
    % ##### Error #####
    error = sum(bitxor(rx_bits,tx_bits)) / numel(rx_bits);
    
end

function error = repetition_system(tx_bits, SNR)

    % Bit Energy
    Eb = 0.5;
    
    % ##### Modulation #####
    msg_symbols = bits2qpsk(tx_bits);
    
    msg_symbols = sqrt(2)*msg_symbols;
    
%     disp(length(msg_symbols));
%     disp(msg_symbols(1:10));

    % ##### Repetition Coding #####
    [tx_1, tx_2] = repetition_encoder(msg_symbols);
    
%     disp(length(tx_1));
%     disp(length(tx_2));
%     disp(tx_1(1:10));
%     disp(tx_2(1:10));
    
    % ##### Channel ##### 
    [rx_symbols, h1, h2] = channel(tx_1, tx_2, SNR, Eb);
    
%     disp(length(h1));
%     disp(length(h2));
%     disp(length(rx_symbols));
%     disp(h1(1:10));
%     disp(h2(1:10));
%     disp(rx_symbols(1:10));
    
    % ##### getting h vector #####  
    h1 = reshape(h1, 2, numel(h1)/2);
    h1 = h1(1,:);
    h2 = reshape(h2, 2, numel(h2)/2);
    h2 = h2(1,:);
    h = [h1' h2'];
    h = h';
    h = h(:);
    
%     disp(length(h));
%     disp(h(1:10));
        
    % ##### Maximal Ratio Combiner #####
    rx_vector = maximal_ratio_combiner(rx_symbols, h, 2);
    
%     disp(rx_vector(1:10));
    
    % ##### QPSK Demodulation #####
    detected_symbols = qpsk_detector(rx_vector);
    
%     disp(detected_symbols(1:10));
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

function [tx_1, tx_2] = repetition_encoder(msg_symbols)

    tx_1 = [msg_symbols zeros(numel(msg_symbols),1)];
    tx_1 = tx_1.';
    tx_1 = tx_1(:);
    
    tx_2 = [zeros(numel(msg_symbols),1) msg_symbols];
    tx_2 = tx_2.';
    tx_2 = tx_2(:);
    
end

function [rx_symbols, h1, h2] = channel(tx_1, tx_2, SNR, Eb)

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
    sigma = sqrt(Eb/(2*snr_now)); % standard devation of noise for BPSK
        
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
        
        rx_symbol = (h_i' * y_i) / (h_i_norm);
        
        rx_vector(i) = rx_symbol;
        
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
        elseif phase(i) >= 90 && phase(i) <= 180
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
