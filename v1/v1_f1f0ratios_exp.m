function v1_f1f0ratios(net,varargin)
% v1_f1f0ratios(data)
%  data : network data 
% shows histograms for f1/f0 ratios for each layer
% assumes data to contain unit props "f1f0ratio" and "optorient"

pr=inputParser;
pr.addParamValue('min_resp',0.01,@isnumeric);
pr.parse(varargin{:});
pr=pr.Results;

v1data=v1_expdata;

W=2; nlayer=2;

subplot(W,3,1);
bar(1.0:0.1:2.1,v1data.f1f0.y(11:end-1)); 
xlim([-0.05 2.15]); set(gca,'XTick',-0.05:0.1:2.15);
set(gca,'XTickLabel',{0,'','','','',0.5,'','','','',1,'','','','',1.5,'','','','',2,''});
set(gca,'FontName','Times','FontSize',12);
box on;

subplot(W,3,2); hold on;
bar(0:10:80,v1data.ori_bw.simple_foveal.y(1:9)/v1data.ori_bw.simple_foveal.N*100);
xlim([-5 185]); set(gca,'XTick',-5:10:175);
set(gca,'XTickLabel',{0,'','',30,'','',60,'','',90,'','',120,'','',150,'','',180});
set(gca,'FontName','Times','FontSize',12);
box on;

subplot(W,3,5); hold on;
bar(0:10:80,v1data.ori_bw.complex_foveal.y(1:9)/v1data.ori_bw.complex_foveal.N*100);
xlim([-5 185]); set(gca,'XTick',-5:10:175);
set(gca,'XTickLabel',{0,'','',30,'','',60,'','',90,'','',120,'','',150,'','',180});
set(gca,'FontName','Times','FontSize',12);
box on;
%     subplot(W,3,nlayer*3+3); hold on;
%     plot(1:12,pr.v1data.peak_freq.simple_foveal.y(1:12)/pr.v1data.peak_freq.simple_foveal.N*100,'r');
%     subplot(W,3,nlayer*3+4); hold on;
%     plot(1:12,pr.v1data.peak_freq.complex_foveal.y(1:12)/pr.v1data.peak_freq.complex_foveal.N*100,'r');

subplot(W,3,4);
bar(0:0.1:0.9,v1data.f1f0.y(1:10)); 
xlim([-0.05 2.15]); set(gca,'XTick',-0.05:0.1:2.15);
set(gca,'XTickLabel',{0,'','','','',0.5,'','','','',1,'','','','',1.5,'','','','',2,''});
set(gca,'FontName','Times','FontSize',12)
box on;

subplot(W,3,3); hold on;
bar(0.4:0.2:2.6,v1data.freq_bw.simple_foveal.y(1:12)/v1data.freq_bw.simple_foveal.N*100);
xlim([0.1 2.7]); set(gca,'XTick',0.3:0.2:2.5);
set(gca,'XTickLabel',{0.2,'','',0.8,'','',1.2,'','',2.0,'','',2.6});
set(gca,'FontName','Times','FontSize',12);
box on;

subplot(W,3,6); hold on;
bar(0.4:0.2:2.6,v1data.freq_bw.complex_foveal.y(1:12)/v1data.freq_bw.complex_foveal.N*100);
xlim([0.1 2.7]); set(gca,'XTick',0.3:0.2:2.5);
set(gca,'XTickLabel',{0.2,'','',0.8,'','',1.2,'','',2.0,'','',2.6});
set(gca,'FontName','Times','FontSize',12);
box on;

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
