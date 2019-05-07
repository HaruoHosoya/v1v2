function dataset=create_image_dataset(imset,width,ndata,varargin)

pr=inputParser;
pr.addParamValue('normalize',false,@islogical);
pr.addParamValue('minvar',0,@isnumeric);
pr.parse(varargin{:});
pr=pr.Results;


% IMAGES=load('~/workspace/bv/matlab/sparsenet/IMAGES.mat','IMAGES');
% IMAGES=IMAGES.IMAGES;
% dataset=load('~/workspace/bv/matlab/sparsenet/IMAGES_RAW.mat','IMAGESr');
% dataset=dataset.IMAGESr;

if ndims(imset)==3
    [imheight,imwidth,imnum]=size(imset);
    nch=1;
elseif ndims(imset)==4
    [imheight,imwidth,nch,imnum]=size(imset);
end;
imset=reshape(imset,imheight,imwidth,nch,imnum);

if pr.normalize
    imset=double(imset)/255;
end;


% ndata=10000;
datadim=width^2*nch;
dataset=zeros(datadim,ndata);

for I=1:ndata
    while(true)
        x=randi(imwidth-width+1); y=randi(imheight-width+1); t=randi(imnum);
        im=imset(y:y+width-1,x:x+width-1,:,t);
        if var(im(:))>pr.minvar break; end;
    end;
    dataset(:,I)=im(:);
end;

end
