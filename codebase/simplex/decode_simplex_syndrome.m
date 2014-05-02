function cdwrd_est = decode_simplex_syndrome(mm,h,actual_syndromes)
%encodes simplex code using syndrome table

s=mod(mm*h',2);
%cdwrd_est=mm;
if (s == [0 0 0 0])
   cdwrd_est=  mod(mm+  [0 0 0 0 0 0 0],2);
elseif (s== [1 0 0 0])
   cdwrd_est= mod(mm+ [1 0 0 0 0 0 0],2);
elseif (s== [0 1 0 0])
   cdwrd_est= mod(mm+ [0 1 0 0 0 0 0],2);
elseif (s==  [0 0 1 0])
    cdwrd_est= mod(mm+ [0 0 1 0 0 0 0],2);
elseif (s== [0 0 0 1])
     cdwrd_est= mod(mm+ [0 0 0 1 0 0 0],2);
elseif (s== [0 1 1 1])
     cdwrd_est= mod(mm+ [0 0 0 0 1 0 0],2);
elseif (s== [1 0 1 1 ])
     cdwrd_est= mod(mm+ [0 0 0 0 0 1 0],2);
elseif (s== [1 1 0 1])
     cdwrd_est= mod(mm+ [0 0 0 0 0 0 1],2);
elseif (s== [1 1 0 0])
     cdwrd_est=  mod(mm+ [1 1 0 0 0 0 0],2);
elseif (s== [1 0 1 0])
     cdwrd_est= mod(mm+ [1 0 1 0 0 0 0],2);
elseif (s== [0 1 1 0])
     cdwrd_est= mod(mm+ [0 1 1 0 0 0 0],2);
elseif (s==  [1 0 0 1])
     cdwrd_est= mod(mm+ [1 0 0 1 0 0 0],2);
elseif (s==  [0 1 0 1])
     cdwrd_est= mod(mm+ [0 1 0 1 0 0 0],2);
elseif (s==  [0 0 1 1])
    cdwrd_est= mod(mm+ [0 0 1 1 0 0 0],2);
elseif (s==  [1 1 1 1])
    cdwrd_est=  mod(mm+ [1 0 0 0 1 0 0],2);
elseif (s==  [1 1 1 0])
    cdwrd_est=  mod(mm+ [1 1 1 0 0 0 0],2);
end
%cdwrd_est
%for i=1:size(actual_syndromes,1)
%    if s==cell2mat(actual_syndromes(i,1))
        %cdwrd_est=mod(mm+cell2mat(actual_syndromes(i,2)),2);
%        cdwrd_est=cell2mat(actual_syndromes(i,2));
%        break
%    end
%end       
