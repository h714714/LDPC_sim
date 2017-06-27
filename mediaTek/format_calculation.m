% calculation of the padding and puncture
Z = 64;
CBS = 528;

cr = [0.89, 0.83 0.75 0.67 0.5 0.4 0.33];
N = [4 6 8 10 18 26 34];

info_data = CBS - 2*Z
info_padding = 16*Z - CBS
parity_data = round(CBS./cr - (CBS-2*Z))
parity_punc = N.*Z - parity_data