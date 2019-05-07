function v1_f1f0ratios(net,varargin)
% v1_f1f0ratios(data)
%  data : network data 
% shows histograms for f1/f0 ratios for each layer
% assumes data to contain unit props "f1f0ratio" and "optorient"

pr=inputParser;
pr.addParamValue('min_resp',0.01,@isnumeric);
pr.addParamValue('layers',2:length(net.structure.layers),@isnumeric);
pr.parse(varargin{:});
pr=pr.Results;

num_layers=length(net.structure.layers);


nlayer=length(pr.layers);
W=nlayer;

for I=1:nlayer
    L=pr.layers(I);
    resps=net.content.layers{L}.unitProperties3.grating;
    resps=max(reshape(resps,size(resps,1)*size(resps,2)*size(resps,3),size(resps,4)*size(resps,5)));
    min_resp=pr.min_resp;
    valid=find(resps(:)>min_resp);
    f1f0ratio=net.content.layers{L}.unitProperties.f1f0ratio;
    optorient=net.content.layers{L}.unitProperties.optorient;
    orientbw=net.content.layers{L}.unitProperties.orientbw;
    optfreq=net.content.layers{L}.unitProperties.optfreq;
    freqbw=net.content.layers{L}.unitProperties.freqbw;
    orients=net.content.layers{L}.layerProperties.orients(:);
    freqs=net.content.layers{L}.layerProperties.freqs(end:-1:1)';
    subplot(W,3,I*3-2);
    histe(f1f0ratio(valid),0:0.1:2); 
    set(gca,'XTickLabel',{0,'','','','',0.5,'','','','',1,'','','','',1.5,'','','','',2});
    set(gca,'FontName','Times','FontSize',12)
%     title(['layer' int2str(L-1) ': f1/f0']);
    subplot(W,3,I*3-1);
    histe(orientbw(valid)*180/pi,0:10:180); 
    set(gca,'XTickLabel',{0,'','',30,'','',60,'','',90,'','',120,'','',150,'','',180});
    set(gca,'FontName','Times','FontSize',12)
%     title('orientation bandwidth');
    subplot(W,3,I*3);
    histe(freqbw(valid),0.2:0.2:2.6); 
    set(gca,'XTickLabel',{0.2,'','',0.8,'','',1.2,'','',2.0,'','',2.6});
    set(gca,'FontName','Times','FontSize',12)
%     title('frequency bandwidth');
    fprintf('layer %d: simple %d, complex %d, meanfreq=%1.2f\n', L-1, sum(f1f0ratio(valid)>=1), sum(f1f0ratio(valid)<1), nanmean(optfreq(valid)));

%     subplot(5,W,W*3+I);
%     histe(optorient(valid)*180/pi,0:180/24:180); title('orientation');
%     set(gca,'XTick',(0:180/4:180)-180/48,'XTickLabel',0:180/4:180);
%     subplot(5,W,W*4+I);
%     histloge(optfreq(valid)*freqscale,[0 0.5 0.7 1.0 1.4 2.0 2.8 4.0 5.6 8.0 11.2 16.0]); title(sprintf('frequency [scale=%1.3f]',freqscale));
end;


end

function histe(x,edges,varargin)
    N=histc(x,edges);
    tick=edges(2)-edges(1);
    bar(edges,N/sum(N)*100,'stacked',varargin{:});
    xlim([edges(1)-tick/2,edges(end)+tick/2]);
    set(gca,'XTick',edges-tick/2,'XTickLabel',edges);
end

function histloge(x,edges,varargin)
    N=histc(x,edges);
    tick=edges(2)-edges(1);
    bar(1:length(edges),N/sum(N)*100,'stacked',varargin{:});
    xlim([0.5,length(edges)+0.5]);
    set(gca,'XTickLabel',num2str(edges(:),'%5.1f'));

end

function histe2(x1,x2,edges,varargin)
    N1=histc(x1,edges);
    N2=histc(x2,edges);
    tick=edges(2)-edges(1);
    bar(edges,[N1(:) N2(:)],'stacked',varargin{:});
    xlim([edges(1)-tick/2,edges(end)+tick/2]);
    set(gca,'XTick',edges-tick/2,'XTickLabel',edges);
    colormap([0 0 0.5625; 1 1 1]);
end
