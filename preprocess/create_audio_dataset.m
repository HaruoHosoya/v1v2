function dataset=create_audio_dataset(psgset,wid,dur,N)

[flen,tlen,naud]=size(psgset);

dataset=zeros(wid,wid,N);
fidx=floor((1:wid)*flen/wid);
tidx=floor((1:wid)*dur/wid);

for I=1:N
    t=randi(tlen-dur);
    n=randi(naud);
    d=psgset(:,t:t+dur-1,n);
    d=d(fidx,tidx);
    dataset(:,:,I)=d;
end;

dataset=reshape(dataset,wid^2,N);

end

    