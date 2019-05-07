function net=v2_anzai_stat(net,varargin)

pr=inputParser;
pr.addParamValue('layer',4,@isnumeric);
pr.addParamValue('minresp',0.5,@isnumeric);
pr.addParamValue('minrsquare1',0.5,@isnumeric);
pr.addParamValue('minrsquare2',0.5,@isnumeric);
pr.addParamValue('addSingle',false,@islogical);
pr.addParamValue('onlyExperim',false,@islogical);
pr.addParamValue('showTypes',false,@islogical);
pr.addParamValue('display',true,@islogical);
pr.parse(varargin{:});
pr=pr.Results;

L=pr.layer;
verbose=false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vmparams=net.content.layers{L}.unitProperties2.vmparams;
[~,numLocation,numUnits,numNodes]=size(vmparams);
numPairs=(numLocation*2)*(numLocation*2-1)/2;

oridiff=ones(numPairs,numUnits,numNodes)*NaN;
maxori=zeros(numUnits,numNodes);

numSingle=0;

for t = 1:numNodes
  for u = 1:numUnits
    peak_oris=ones(2,numLocation)*NaN;
    for v = 1:numLocation
      ps=num2cell(vmparams(:,v,u,t));
      [rs1,a1,m1,s1, rs2,a21,m21,s21,a22,m22,s22]=ps{:};      
      if (rs1 > pr.minrsquare1)
          a=vonmises(a1,m1,s1,s1);
          if a >= pr.minresp peak_oris(1,v)=s1; end; 
      elseif rs2 > pr.minrsquare2
          a1=vonmises(a21,m21,s21,s21)+vonmises(a22,m22,s22,s21);
          a2=vonmises(a21,m21,s21,s22)+vonmises(a22,m22,s22,s22);
          if a1 >= pr.minresp peak_oris(1,v)=s21; end;
          if a2 >= pr.minresp peak_oris(2,v)=s22; end;
      end;
    end;
    peak_oris=peak_oris(:);
    c=combnk(1:length(peak_oris),2);
    diffs=orientationDiff(peak_oris(c(:,1)),peak_oris(c(:,2)));
    [~,ix]=max(abs(diffs));
    maxori(u,t)=diffs(ix);
    oridiff(:,u,t)=diffs;
    if all(isnan(peak_oris)) && verbose fprintf('empty (%d,%d)\n',t,u); end;
    if sum(~isnan(peak_oris))==1 
        if pr.addSingle maxori(u,t)=0; end;
        if verbose fprintf('singleton (%d,%d)\n',t,u); end;
        numSingle=numSingle+1; 
    end;
  end;
end;

maxori=abs(maxori);
oridiff=abs(oridiff);
oridiff(:,maxori<pi/6)=NaN;
maxori_orig=maxori;

net.content.layers{L}.unitProperties.maxori=reshape(maxori,1,numUnits,numNodes);
net.content.layers{L}.unitProperties.oridiff=reshape(oridiff,numPairs,numUnits,numNodes);


if pr.showTypes
    inh_types=net.content.layers{L}.unitProperties.inh_types;
    desc_valid=net.content.layers{L}.unitProperties.desc_valid;
    exc_types=net.content.layers{L}.unitProperties.exc_types;
    maxori2=ones(numUnits,numNodes,5)*NaN;
%     maxori2(inh_types<=3,1)=maxori(inh_types<=3);
%     maxori2(inh_types==4,2)=maxori(inh_types==4);
%     maxori2(exc_types<=2,1)=maxori(exc_types<=2);
%     maxori2(exc_types==3,2)=maxori(exc_types==3);
    class5=exc_types==3&inh_types==4;
    class4=exc_types<=2&inh_types==4;
    class3=inh_types==3;
    class2=inh_types==2;
    class1=inh_types==1;
    maxori(desc_valid~=1)=NaN;
    oridiff(desc_valid~=1)=NaN;
    maxori2(class1,1)=maxori(class1);
    maxori2(class2,2)=maxori(class2);
    maxori2(class3,3)=maxori(class3);
    maxori2(class4,4)=maxori(class4);
    maxori2(class5,5)=maxori(class5);
    
%     oridiff2=ones(numPairs,numUnits,numNodes,2)*NaN;
%     oridiff2(:,exc_types<=2,1)=oridiff(:,exc_types<=2);
%     oridiff2(:,exc_types==3,2)=oridiff(:,exc_types==3);

%     maxori2(desc_valid~=1,:)=NaN;
%     oridiff2(:,desc_valid~=1,:)=NaN;

    net.content.layers{L}.unitProperties.maxori=reshape(maxori,1,numUnits,numNodes);
    net.content.layers{L}.unitProperties.oridiff=reshape(oridiff,numPairs,numUnits,numNodes);

    maxori=maxori2;
%     oridiff=oridiff2;

end;

maxori=reshape(maxori,numNodes*numUnits,size(maxori,3));
oridiff=reshape(oridiff,numPairs*numNodes*numUnits,size(oridiff,4));

if ~pr.display return; end;

fprintf('singletons %d\n', numSingle);

% display

% figure('Position',[0 0 800 600]);

subplot(2,2,1);
h=hist(isnan(nansum(maxori,2)),[0,1]);
pie(h);
title(['used/discarded N=' int2str(sum(h))]);

subplot(2,2,2);
h=histc(nansum(maxori,2),[0,pi/6,pi]);
pie(h(1:2));
title(['Unit groups (uniform/non-uniform) N=' int2str(sum(h))]);

X=0:15:90;

subplot(2,2,3);
Y1=histc(maxori/pi*180,X);
Y2=[121 0 51 0 189 414 345 261 19 51 51 225 0]';
Y2=floor(Y2/sum(Y2(:))*100);  %N=100
Y2=[Y2(6:-1:1)+Y2(7:12); 0];
if pr.onlyExperim
    bar(X,[Y1/sum(Y1(:))*0,Y2/sum(Y2(:))*100],'histc');
else
%     bar(X,[Y1/sum(Y1(:))*100,Y2/sum(Y2(:))*100],'histc');
    b=bar(X+15/2,Y1/sum(Y1(:))*100,'stacked'); hold on;
    plot(X(1:end-1)+15/2,Y2(1:end-1)/sum(Y2(:))*100,'mx-','LineWidth',2);
    if pr.showTypes
        set(b,{'FaceColor'},{'b';'c';'y';'r';[0.5 0 0]}); 
        legend('broad','side','cross','end (aligned)','end (converging)');
    end;
end;
set(gca,'XLim',[0,90]);
set(gca,'XTick',X);
gof=hellingerdist(Y1/sum(Y1),Y2/sum(Y2));
% title(sprintf('(a) Maximum orientation differences [h=%1.3f]',gof));
xlabel('degree');
ylabel('# of units (%)');
set(gca,'FontName','Times','FontSize',12);

subplot(2,2,4);
Y1=histc(oridiff/pi*180,X);
Y2=[141 44 29 66 78 263 410 122 46 46 53 220 0]';
Y2=floor(Y2/sum(Y2(:))*424);  %N=424
Y2=[Y2(6:-1:1)+Y2(7:12); 0];
if pr.onlyExperim
    bar(X,[Y1/sum(Y1(:))*0,Y2/sum(Y2(:))*100],'histc');
else
%     bar(X,[Y1/sum(Y1(:))*100,Y2/sum(Y2(:))*100],'histc');
    bar(X+15/2,Y1/sum(Y1(:))*100,'stacked'); hold on;
    plot(X(1:end-1)+15/2,Y2(1:end-1)/sum(Y2(:))*100,'mx-','LineWidth',2);
end;
set(gca, 'XLim', [0,90]);
set(gca,'XTick',X);
gof=hellingerdist(Y1/sum(Y1),Y2/sum(Y2));
% title(sprintf('(b) All orientation differences for non-uniforms [h=%1.3f]',gof));
xlabel('degree');
% ylabel('occurrence (%)');
set(gca,'FontName','Times','FontSize',12);



% if pr.showTypes
%     subplot(3,2,5);
%     exc_types=net.content.layers{L}.unitProperties.exc_types;
%     inh_types=net.content.layers{L}.unitProperties.inh_types;
%     params=net.content.layers{L}.unitProperties2.desc_params;
%     rsq=net.content.layers{L}.unitProperties.rsq;
% %     e=exc_types==3 & rsq>0.5;
%     e=inh_types==4 & rsq>0.5;
%     ps=params(:,4,e);
%     mo=maxori_orig(e)*180/pi;
% %     divg=atan2(ps(5,:),ps(7,:)-ps(6,:))*2*180/pi;
%     divg=ps(7,:);
%     scatter(mo(:),divg(:));
%     [co,p]=corr(mo(:),divg(:));
%     xlabel('max orientation diff');
%     ylabel('divergence');
%     title(sprintf('[corr=%1.3f p=%1.3e]',co,p));
% end;

end

