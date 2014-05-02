function[cdwrd] = encode_syndrome(mm,g)
cdwrd = mod(mm*g,2);
