function psgset=create_spectrogram_from_cmu_arctic(indir)

files=dir([indir '/*.wav']);
nfiles=length(files);
lens=zeros(nfiles,1);

for I=1:nfiles
    fprintf('processing %s...\n',files(I).name);
    aud=audioread([indir '/' files(I).name]);
    aud=downsample(aud,4);
    [S,F,T,P]=spectrogram(aud,256,250,256,4E3);
    sgset{I}.power=P;
    lens(I)=size(P,2);
end;

len=min(lens);
psgset=zeros(size(P,1),len,nfiles);

for I=1:nfiles
    psgset(:,:,I)=10*log10(sgset{I}.power(:,1:len));
end;


end


