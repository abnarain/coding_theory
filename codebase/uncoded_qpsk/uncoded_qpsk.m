% Model parameters
%close all 
%clear all
Es = 1; % energy per QPSK symbol
n=10000;
k=n;
Eb = (Es*n)/(2*k); % energy per coded bit
Eb_by_N0_db = linspace(-2,10,13);
N0= Eb./(10.^(Eb_by_N0_db/10))
[row,col]=size(N0);

global_ber=[]
[row,col]=size(N0);
for n0=1:col
    col
sigma2 = N0(n0)/2; % noise variance
% Simulation of single transmission
msg = (rand(1,k)>0.5); % generates k bits uniformly at random
cdwrd = msg;%encoded(msg); % encodes message msg into codeword cdwrd of n bits
col_=n;
row_=0;
if mod(n,2)==1
    cdwrd = [cdwrd 0]; % padding if odd length
    [row_,col_]=size(cdwrd);
end
symb = sqrt(Es/2)*(1-2*cdwrd(1:2:end))+j*sqrt(Es/2)*(1-2*cdwrd(2:2:end)); % QPSK modulation
rcvd = symb + sqrt(sigma2)*(randn(size(symb))+j*randn(size(symb))); % channel noise
rcvd_dem = zeros(1,2*length(symb)); % demodulation with optimal decision
rcvd_dem(1:2:end) = (real(rcvd)<0); % hard decoding of in-phase (I) component
rcvd_dem(2:2:end) = (imag(rcvd)<0); % hard decoding of out-of-phase (Q) component
rcvd_dem = rcvd_dem(1:n); % removes the padded bits
cdwrd_est = rcvd_dem ; % decodes the demodulated bit sequence to form an estimate cdwrd_est
ber = sum(msg~=cdwrd_est)/n % computes the BER over the codeword
global_ber=[global_ber,ber];
end
semilogy(Eb_by_N0_db,global_ber,'b+');
set(gca,'fontsize',12);
h=xlabel('Eb/N0 (dB)');set(h,'fontsize',12);
h=ylabel('Bit Error Rate');set(h,'fontsize',12);
title('BER vs Power Spectral Efficiency of Uncoded QPSK')
