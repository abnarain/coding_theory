function cdwrd_est = decode_bounded_distance(m,h,codewords)
m_l = repmat(m,size(codewords,1),1);
new_list = mod (m_l + codewords,2);
sb=sum(new_list,2);
[m_val,index]=min(sb);
%index=0;
%for i=1:size(sb)
%    if sb[i]==m_val
%        index=i
%        disp('index');
%        break
%    end
%end        
%index = find(min(sum(new_list,2)))%,2)) %<=1
if m_val>3
    cdwrd_est = m;
else
    cdwrd_est = codewords(index,:);
end
