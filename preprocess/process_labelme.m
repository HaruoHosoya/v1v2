function images=process_labelme(indir,varargin)

pr=inputParser;
pr.parse(varargin{:});
pr=pr.Results;

subdirs=dir([indir '/*']);
files={}; K=1;
for I=1:length(subdirs)
    if subdirs(I).isdir && subdirs(I).name(1)~='.' 
        files1=dir([indir '/' subdirs(I).name '/*.jpg']);
        for J=1:length(files1)
            if files1(J).name(1)~='.'
                files{K}=[indir '/' subdirs(I).name '/' files1(J).name];
                K=K+1;
            end;
        end;
    end;
end;

nfiles=length(files);
img1=imread(files{1});
[height,width,nch]=size(img1);
images=zeros(height,width,nch,nfiles,'uint8');

for I=1:nfiles
    img=imread(files{I});
    images(:,:,:,I)=img;
end;


end

    
    


