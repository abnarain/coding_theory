% For Simplex Codes
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

g = [ [0,1,1,1,1,0,0];
      [1,0,1,1,0,1,0];
      [1,1,0,1,0,0,1]];

h = [[1,0,0,0,0,1,1];
     [0,1,0,0,1,0,1];
     [0,0,1,0,1,1,0];
     [0,0,0,1,1,1,1]]; 

max_wt=7;
errors=[];
for i=1:max_wt
	combinations=combnk(1:n,i)
	[row_c,col_c]=size(combinations)
	for jk=1:row_c
		arr=zeros(1,n);
        for jj=1:col_c
    	 arr(n-combinations(jk,jj)+1)=1;
        end
	    errors = [errors; arr];
     end
end

dd= bi2de(errors,'left-msb');
b= dec2bin(dd);
str2num(b);
errors=arrayfun(@str2double,b);
actual_syndromes={};
index=1
for r=1:size(errors,1)
    err_elem=errors(r,:);
    synd= mod(err_elem*h',2);
    flag =0;
    size(actual_syndromes);
    for i=1:size(actual_syndromes,1)
        if synd== cell2mat(actual_syndromes(i,1))
            flag=1;
             break;
        end
    end
    if flag ==0
        actual_syndromes(index,1)={synd};
        actual_syndromes(index,2)={err_elem};
        index=index+1
        r
    end
end

for n0=1:col
sigma2=N0(n0)/2;
n0
frame_err = 0; %serialize the number of frame error s
bit_err = 0; % initialize the number of bit errors
frame_nb = 0; % initializethe number of codeword transmissions
pl_f=0; %initialize lower bound on Wilson interval;
pl_b=0; %initialize lower bound on Wilson interval;
pu_f=1; %initialize upper bound on Wilson interval;
pu_b=1; %initialize upper bound on Wilson interval;
m=0;
while(1) %(frame_nb<10000)||(abs(pu_f-pl_f)>0.1*frame_err/frame_nb))
    frame_nb = frame_nb +1; % update count of frame simulations
    m=frame_nb
    msg = (rand(1,k)>0.5); % generates k bits uniformly at random
    cdwrd = encode_simplex_syndrome(msg,g) ;% encodes message m into codeword c of n bits
    if mod(n,2)==1
        cdwrd = [cdwrd 0]; % padding if odd length
    end
    symb = sqrt(Es/2)*(1-2*cdwrd(1:2:end))+j*sqrt(Es/2)*(1-2*cdwrd(2:2:end)); % QPSK modulation
    rcvd = symb + sqrt(sigma2)*(randn(size(symb))+j*randn(size(symb))); % channel noise
    rcvd_dem = zeros(1,2*length(symb)); % demodulation with optimal decision
    rcvd_dem(1:2:end) = (real(rcvd)<0); % hard decoding of in-phase (I) component
    rcvd_dem(2:2:end) = (imag(rcvd)<0); % hard decoding of out-of-phase (Q) component
    %done with demodulation
    rcvd_dem = rcvd_dem(1:n); % removes the padded bits
    cdwrd_est = decode_simplex_syndrome(rcvd_dem,h,actual_syndromes); % decodes
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
    pu_f=(p_hat_fer +(constant_2/2)+(d_alpha*determ_f))/ (1+constant_2);
    pl_f=(p_hat_fer +(constant_2/2)-(d_alpha*determ_f))/ (1+constant_2);
    
    pu_b=(p_hat_ber +(constant_2/2)+(d_alpha*determ_b))/ (1+constant_2);
    pl_b=(p_hat_ber +(constant_2/2)-(d_alpha*determ_b))/ (1+constant_2);
    if (frame_nb >100000)
        break
    end    
end
global_ber=[global_ber,ber];
global_fer=[global_fer,fer];
CI_ber=[CI_ber,pu_b-pl_b];
CI_fer=[CI_fer,pu_f-pl_f];

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
title('Simplex Code (Sydrome decoding)')
