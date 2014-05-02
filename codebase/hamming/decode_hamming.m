function [cdwrd_est] = decode_hamming(m,h)%,decoding_table)
%encodes hamming code
syndrome = mod(m * h',2);
% c =sum(syndrome);
% flag=0;
a = [ 4 2 1];
pos = sum(syndrome.*a);
if pos~=0
m(pos) = ~m(pos);
end
cdwrd_est = m;
% if sum(syndrome) ==1
%     flag=0;
%     m(pos) = ~m(pos);
%     m=mod(m,2);
%     cdwrd_est = m;
% elseif sum(syndrome)==2
%     flag=0;
%     m(sum(syndrome)) =m(sum(syndrome))+1;
%     m=mod(m,2);
%     cdwrd_est = m;
% elseif sum(syndrome)==4
%     flag=0;
%     m(sum(syndrome)) =m(sum(syndrome))+1;
%     m=mod(m,2);
%     cdwrd_est = m;
% else
%     flag=1;
%     cdwrd_est = [5];
% end

% pause