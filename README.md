# EE 404: Wireless Communication

**1. Implement a BPSK modulated signal transmission adopting coherent
reception.**

    a. Obtain the BER vs. SNR. curves with L = 1 (no repetition) and L =
{2,3,4,5} repetitions over a Rayleigh fading channel. Determine the
range of SNR values you want to simulate using your own intuition.

    b. Let us say that each coherence period corresponds to 8 BPSK
samples. Create your own random channel which becomes “bad” once
in every four coherence time periods (you can define “bad” based on
intuition). Show through BER vs. SNR simulations that "interleaving +
repetition" has better performance than "only interleaving" or "only
repetition".

**2. Implement the 2 Tx x 1 Rx Alamouti scheme using QPSK modulated symbols.
Plot the BER vs. SNR for a Rayleigh fading channel with and without Alamouti
scheme. Compare the simulation results with the theoretical upper bound for
the probability of error vs. SNR**

**3. Build your own point-to-point OFDM system (transmitter and receiver) with 64-
point IFFT/FFT and 8 samples of CP. The bandwidth is 1 MHz. Each
subcarrier should carry QPSK modulated signals.**

    a. Plot BER vs. SNR for a single tap Rayleigh fading channel assuming
that the channel is perfectly known.

    b. Come up with your own pilot transmission scheme. Perform your own
channel estimation at the receiver such that it works for frequency
selective channels. Plot BER vs. SNR for a 4-tap Rayleigh fading
channel. 
