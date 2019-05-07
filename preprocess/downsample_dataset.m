function dataset=downsample_dataset(dirname,wid,len)

dataset=zeros(wid,wid,len);

files=dir([dirname '/*.h5']);
nfiles=length(files);
nsample_per_file=ceil(len/nfiles);

pos=1;
for I=1:nfiles
    fprintf('loading %s...\n',files(I).name);
    data=h5read([dirname '/' files(I).name],'/data');
    datalen=size(data,3);
    dataset(:,:,pos:pos+nsample_per_file-1)=data(:,:,randi(datalen,nsample_per_file,1));
    pos=pos+nsample_per_file;
end;

end


    