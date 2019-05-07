function joint_dist_table(r,xunits,yunits,varargin)

pr=inputParser;
pr.addParamValue('style','histogram',@isstr);
pr.addParamValue('maxval',10,@isnumeric);
pr.parse(varargin{:});
options=pr.Results;

nx=length(xunits);
ny=length(yunits);

m=options.maxval;
P=1;
for y=1:ny
    for x=1:nx
        subplot(ny,nx,P);P=P+1;
        switch(options.style)
            case 'histogram'
                cnt=hist3([r(:,xunits(x)) r(:,yunits(y))],{-m:0.5:m,-m:0.5:m});
                imagesc(log10(cnt)');
                axis xy;axis square;set(gca,'xtick',[]);set(gca,'ytick',[]);
            case 'cond-histogram'
                cnt=hist3([r(:,xunits(x)) r(:,yunits(y))],{-m:0.25:m,-m:0.25:m});
                c=bsxfun(@rdivide,cnt,sum(cnt,2)); 
                imagesc(-m:0.25:m,-m:0.25:m,c',[0 0.2]);
                axis xy;axis square;set(gca,'xtick',[]);set(gca,'ytick',[]);
            case 'scatter'
                scatter(r(:,xunits(x)),r(:,yunits(y)));
                axis xy;axis square;set(gca,'xtick',[]);set(gca,'ytick',[]);
            case 'contour'
                cnt=hist3([r(:,xunits(x)) r(:,yunits(y))],{-m:0.5:m,-m:0.5:m});
                s=log10(cnt)'; s(isinf(s))=0;
                contour(s,5);
                axis xy;axis square;set(gca,'xtick',[]);set(gca,'ytick',[]);
        end;                
        if y==ny xlabel(num2str(xunits(x))); end;
        if x==1 ylabel(num2str(yunits(y))); end;
    end;
end;

colormap(jet);

end