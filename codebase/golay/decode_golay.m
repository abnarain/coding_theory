function [cdwrd_est,flag] = decode_hamming(m,h,actual_syndromes)
%encodes hamming code
s= mod(m * h',2);
flag=0;
cdwrd_est= m;
for i=1:size(actual_syndromes,1)
    if s==cell2mat(actual_syndromes(i,1))
        cdwrd_est=mod(mm+cell2mat(actual_syndromes(i,2)),2);
        break
    end
end 
