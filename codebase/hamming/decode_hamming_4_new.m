function [cdwrd_est] = decode_hamming_3_new(m,h,a)%,decoding_table)
%encodes hamming code
syndrome = mod(m * h',2);
% c =sum(syndrome);
% flag=0;
pos = sum(syndrome.*a,2);
%size(pos,1)
%pause
%if pos~=0
%m(pos) = ~m(pos);
%end

for jj=1:size(pos,1)
    if pos(jj) ~=0
    m(jj,pos(jj)) = ~m(jj,pos(jj));
    end
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
