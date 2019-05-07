function bank=test_bank

wid=32;
% npos=5;
npos=6;
% nori=4;
nori=12;
pos=wid/(npos+1):wid/(npos+1):wid-wid/(npos+1);
oris=0:pi/nori:pi-pi/nori;
% freqs=[1/8 1/6 1/4];
% freqs=[1/6];
freqs=[1/8];
% phases=0:pi/2:2*pi-pi/2;
phases=[0 pi/2];
bank=gabor_bank(wid,0.4,0.4,pos,freqs,oris,phases);
bank=bank./repmat(sqrt(sum(bank.^2,1)),[size(bank,1),1,1,1,1,1]);

end
