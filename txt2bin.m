function bin = txt2bin(txt)
S = dec2bin(txt,8);
txt_bin = [];
l = size(S);
for i = 1:l(1)
    S_i = S(i,:);
    txt_bin = [txt_bin S_i];
end
bin = [];
for j = 1:length(txt_bin)
    bin_j = str2num(txt_bin(j));
    bin = [bin bin_j];
end