% For (7,3) hamming codes
Es=1
k=3
n=7
format long
Eb = (n*Es)/(2*k); % energy per coded bit
Eb_by_N0_db = linspace(-2,10,13);
N0= Eb./(10.^(Eb_by_N0_db/10))
[row,col]=size(N0);

global_ber=[];
global_fer=[];
CI_ber=[];
CI_fer=[];

d_alpha=1.96;
messages=[[0,0,0];
    [0,0,1];
    [0,1,0];
    [0,1,1];
    [1,0,0];
    [1,0,1];
    [1,1,0];
    [1,1,1]];


for n0=1:col
g = [ [0,1,1,1,1,0,0];
      [1,0,1,1,0,1,0];
      [1,1,0,1,0,0,1]];

h = [[1,0,0,0,0,1,1];
     [0,1,0,0,1,0,1];
     [0,0,1,0,1,1,0];
     [0,0,0,1,1,1,1]]; 

codewords=mod(messages*g,2);
sigma2=N0(n0)/2;
n0
frame_err = 1; %serialize the number of frame error s
bit_err = 0; % initialize the number of bit errors
frame_nb = 1; % initializethe number of codeword transmissions
pl_f=0; %initialize lower bound on Wilson interval;
pl_b=0; %initialize lower bound on Wilson interval;
pu_f=1; %initialize upper bound on Wilson interval;
pu_b=1; %initialize upper bound on Wilson interval;

m=0;

while(1) %abs(pu_f-pl_f)>0.1*frame_err/frame_nb)
    frame_nb = frame_nb +1; % update count of frame simulations
    m=frame_nb
    msg = (rand(1,k)>0.5); % generates k bits uniformly at random
    
    cdwrd = encode_simplex_bd(msg,g); % encodes message m into codeword c of n bits
    if mod(n,2)==1
        cdwrd = [cdwrd 0]; % padding if odd length
    end
    %symb = modulate(cdwrd); % modulates binary digits of c onto complex symbols in constellation
    symb = sqrt(Es/2)*(1-2*cdwrd(1:2:end))+j*sqrt(Es/2)*(1-2*cdwrd(2:2:end)); % QPSK modulation
    %done with modulation
    rcvd = symb + sqrt(sigma2)*(randn(size(symb))+j*randn(size(symb))); % channel noise
    %rcvd_dem = demodulate(rcvd); % demodulates received signal r (hard decoding); done above
    rcvd_dem = zeros(1,2*length(symb)); % demodulation with optimal decision
    rcvd_dem(1:2:end) = (real(rcvd)<0); % hard decoding of in-phase (I) component
    rcvd_dem(2:2:end) = (imag(rcvd)<0); % hard decoding of out-of-phase (Q) component
    %done with demodulation
    rcvd_dem = rcvd_dem(1:n); % removes the padded bits
    cdwrd_est = decode_bounded_distance(rcvd_dem,h,codewords);%decoding_table);% decodes
    if mod(n,2)==1 % removes the extra padded bit to get original cdwrd
        cdwrd=cdwrd(1:n);
    end
    bit_err = bit_err + sum(cdwrd~=cdwrd_est);% updates the bit errors
    frame_err = frame_err + (sum(cdwrd~=cdwrd_est)>0); % updates the frame errors
    fer=frame_err/frame_nb; % updates frame error rate
    p_hat_fer=fer; % updates p_hat_fer
    ber=bit_err/(frame_nb*n); % updates bit error rate
    p_hat_ber=ber; % updates p_hat_ber
    constant=d_alpha/(2*m) ;%constants used below
    constant_2 = (d_alpha^2)/m;
    determ_f=sqrt((p_hat_fer*(1-p_hat_fer)/m) +(constant^2));
    determ_b=sqrt((p_hat_ber*(1-p_hat_ber)/m) +(constant^2));
    pu_f=((p_hat_fer +(constant_2/2)+(d_alpha*determ_f))/ (1+constant_2));
    pl_f=((p_hat_fer +(constant_2/2)-(d_alpha*determ_f))/ (1+constant_2));
    
    pu_b=((p_hat_ber +(constant_2/2)+(d_alpha*determ_b))/ (1+constant_2));
    pl_b=((p_hat_ber +(constant_2/2)-(d_alpha*determ_b))/ (1+constant_2));
    if (frame_nb>1000000)
        break
    end
end
global_ber=[global_ber,ber];
global_fer=[global_fer,fer];
CI_ber=[CI_ber, abs(pu_b-pl_b)];
CI_fer=[CI_fer, abs(pu_f-pl_f)];

end
a=errorbar(Eb_by_N0_db,global_ber,CI_ber)
%hold on;
%errorbar(Eb_by_N0_db,global_fer,CI_fer);
%hold off;
set(gca,'fontsize',12);
h=xlabel('Power Efficiency (dB)');set(h,'fontsize',12);
h=ylabel('Bit Error Rate');set(h,'fontsize',12);
set(gca,'YScale','log');
%legend();
title('Simplex Code (bounded distance)')
