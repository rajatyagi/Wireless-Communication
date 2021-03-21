% Parameters
num_subcarriers = 64;
taps = 4;
signal_length = num_subcarriers * 10000;
SNR = -10 : 0.5 : 30;

h_real = sqrt(0.5)*randn(taps, 1); % real part
h_imag = sqrt(0.5)*randn(taps, 1); % imaginary part

h = h_real + 1i*h_imag; % complex number

%Generating a random signal with binary bits (0 and 1)
tx_bits = randi([0,1],signal_length,1);

ofdm_BER = zeros(numel(SNR),1);

for i = 1:numel(SNR)

    ofdm_BER(i) = ofdm_system(tx_bits, SNR(i), taps, 64, 8, h);
    
end

%% PLOTTING THE BER FOR ALL CASES

figure(1);
semilogy(SNR, ofdm_BER);
title('BER for Different Cases');
xlabel('SNR (dB)');
ylabel('BER');
grid on;

%% FUNCTIONS

function error = ofdm_system(tx_bits, SNR, taps, num_subcarriers, cp_len, h)

    % Bit Energy
    Eb = 0.5;

    % ##### QPSK Modulation #####
    msg_symbols = sqrt(2)*bits2qpsk(tx_bits);
    
    % ##### Serial to Parallel #####
    tx_ofdm_symbols = reshape(msg_symbols, num_subcarriers, numel(msg_symbols)/num_subcarriers);
    
    % ##### IFFT #####
    ifft_symbols = ofdm_ifft(tx_ofdm_symbols);
    
    % ##### Add Cyclic Prefix #####
    cyclic_prefix = ifft_symbols(end - cp_len + 1 : end,:);
    cp_ifft_symbols = vertcat(cyclic_prefix,ifft_symbols);
    
    % ##### Parallel to Serial #####
    tx_symbols = cp_ifft_symbols(:);    
    
    % ##### Channel #####
    [rx_symbols, h] = channel(tx_symbols, SNR, Eb, taps, num_subcarriers, cp_len, h);
    
    % ##### h matrix #####
    sz = size(h);
    H = zeros(numel(rx_symbols)/(num_subcarriers + cp_len),num_subcarriers);
    for i = 1 : sz(1)
        if(mod(i,(num_subcarriers + cp_len)) == 1)
           
            h_buf = h(i,:);
            H(1 + (i-1)/(num_subcarriers + cp_len),:) = h_buf(1:num_subcarriers);
            
        end
    end
    
    H = H.';
    
    % ##### Serial to Parallel #####
    rx_ofdm_symbols = reshape(rx_symbols, num_subcarriers + cp_len, numel(rx_symbols)/(num_subcarriers + cp_len));
    
    % ##### Remove Cyclic Prefix #####
    rx_ofdm_symbols = rx_ofdm_symbols(cp_len + 1 : end, :);
    
    % ##### FFT #####
    fft_symbols = fft(rx_ofdm_symbols);
    fft_H = fft(H);
    
    % ##### Maximal Ratio Combiner #####
    rx_vector = maximal_ratio_combiner(fft_symbols, fft_H);
    
    % ##### QPSK Demodulation ##### 
    detected_symbols = qpsk_detector(rx_vector);
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

function ifft_symbols = ofdm_ifft(ofdm_symbols)

    sz = size(ofdm_symbols);

    ifft_symbols = zeros(sz(1),sz(2));
    for i = 1 : sz(2)
       
        ifft_symbols(:,i) = ifft(ofdm_symbols(:,i));
        
    end

end

function fft_symbols = ofdm_fft(ofdm_symbols)

    sz = size(ofdm_symbols);

    fft_symbols = zeros(sz(1),sz(2));
    for i = 1 : sz(2)
       
        fft_symbols(:,i) = fft(ofdm_symbols(:,i));
        
    end

end

function [rx_symbols, h] = channel(symbols, SNR, Eb, taps, num_subcarriers, cp_len, h1)

    y = zeros(numel(symbols),1);
    symbols_padded = vertcat(zeros(num_subcarriers + cp_len - 1, 1),symbols);
    h = [];
    for i = 1 : (numel(symbols) / (num_subcarriers + cp_len))
       
        h = [h channel_gain(num_subcarriers, taps, cp_len, h1)];
        
    end
    
    h = h.';
    
    for i = num_subcarriers + cp_len : length(symbols_padded)
    
        buf_symbols = flipud(symbols_padded(i - num_subcarriers - cp_len + 1 : i));
        
        y(i - num_subcarriers - cp_len + 1) = h(i - num_subcarriers - cp_len + 1,:) * buf_symbols;
        
    end
    
    % This function emulates a Rayleigh fading channel with
    % Channel gain => h ~ CN(0,1)
    % AWGN => w ~ CN(0, sigma^2)
    % where, sigma = sqrt(1/2*SNR)
        
    % Signal to noise ratio on linear scale for the iteration
    snr_now = 10^(SNR/10);
    
    % ADDITIVE WHITE GAUSSIAN NOISE
    sigma = sqrt(Eb/(2*snr_now)); % standard devation of noise for BPSK
    
    w = sigma*(randn(numel(symbols),1) + 1i*randn(numel(symbols),1)); % complex noise
    
    % Received signal
    rx_symbols =  y + w;

end

function h = channel_gain(num_subcarriers, taps, cp_len, h1)

    % CHANNEL GAIN
    % Rayleigh fading channel gain is given by a 
    % complex normal distribution.
    
%     h_real = sqrt(0.5)*randn(taps, 1); % real part
%     h_imag = sqrt(0.5)*randn(taps, 1); % imaginary part
%     
%     h = h_real + 1i*h_imag; % complex number
    
    h = vertcat(h1,zeros(num_subcarriers + cp_len -taps,1));
    
    h = repmat(h, 1, num_subcarriers + cp_len);

end

function rx_vector = maximal_ratio_combiner(rx_symbols, H)

    H_flat = H(:);
    rx_symbols_flat = rx_symbols(:);
    
    % Maximal Ratio Combiner
    rx_vector = (conj(H_flat) .* rx_symbols_flat) ./ abs(H_flat);

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
