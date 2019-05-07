function v1s_stat( net, db )
% v1s_stat( net, db )
%   net : network data
%   db : experimental data provided by Ringach

width=net.structure.layers{1}.patchWidth;
[x0 y0 amp sigmax sigmay theta freq phi] = getParams(net);
res = getResiduals(net);
resTol = 0.1;

margin = 1;
well_within = x0 >= margin+1 & x0 <= width-margin & y0 >= margin+1 & y0 <= width-margin & ...
 sigmax <= width/2 & sigmay <= width/2  & freq > 1/width/3;
valid = well_within & res <= resTol;


x0(~valid)=NaN;
y0(~valid)=NaN;
amp(~valid)=NaN;
sigmax(~valid)=NaN;
sigmay(~valid)=NaN;
theta(~valid)=NaN;
freq(~valid)=NaN;
phi(~valid)=NaN;


bandwidth=log2((sigmax.*freq*pi + sqrt(log(2)/2))./(sigmax.*freq*pi - sqrt(log(2)/2)));
bandwidth(~imag(bandwidth)==0)=NaN;

aspect=sigmay ./ sigmax;
aspect2=sigmay*5 ./ (1./freq/2);
aspect3=aspect;



disp(sprintf('discarded %d units for bad fiting',sum(res>resTol)));
disp(sprintf('discarded %d units for boundary excess / too-low frequency',sum(~valid)-sum(res>resTol)));
disp(sprintf('use %d units for normal analysis',sum(valid)));
disp(sprintf('discarded %d units for invalid bandwidth',sum(isnan(bandwidth))-sum(~valid)));
disp(sprintf('use %d units for bandwidth analysis',sum(~isnan(bandwidth))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[0 0 800 800]); 

% Spatial frequency

subplot(3,2,1);
X=[0.5 0.7 1.0 1.4 2.0 2.8 4.0 5.6 8.0 11.2 16.0];
Y2=[36 36 69 228 298 238 250 110 46 18 11 0]';    % DeValois et al. 1982 "X cells foveal"
Y2=floor(Y2/sum(Y2)*148);  % N=148
[freq_scale,freq_gof]=scaleFit([X,Inf],freq,Y2);
Y1=histc(freq * freq_scale,[X,Inf]);
bar(1:(length(X')+1),[Y1/sum(Y1)*100,Y2/sum(Y2)*100],'histc');
set(gca, 'XLim', [0 length(X')+2], 'XTickLabel',{num2str(X'); '>16'}); 
title(sprintf('(a) Spatial frequency (cy/deg) [scale=%3.1f h=%1.3f]',freq_scale, freq_gof));
ylabel('occurrence (%)');

% Spatial Frequency Bandwidth

subplot(3,2,3);
X=(0.0:0.5:6)';
Y1=histc(bandwidth,X);
bandwidth2=log2((db.nx'.*pi + sqrt(log(2)/2))./(db.nx'*pi - sqrt(log(2)/2)));
bandwidth2(~imag(bandwidth2)==0)=NaN;
Y2=histc(bandwidth2,X);
bandwidth_gof=goodnessOfFit(X,bandwidth,Y2,1);
bar(1:length(X'),[Y1/sum(Y1)*100,Y2/sum(Y2)*100],'histc');
set(gca, 'XLim', [0 length(X')+2], 'XTickLabel',{num2str(X)}); 
title(sprintf('(c) Spatial frequency bandwidth (octave) [h=%1.3f]',bandwidth_gof));
ylabel('occurrence (%)');

% X=0.4:0.2:2.6;
% Y1=histc(bandwidth,[X,Inf]);
% Y2=[12 77 96 161 193 166 96 128 65 65 65 56 0]';    % DeValois et al. 1982 "X cells foveal"
% Y2=floor(Y2/sum(Y2)*145);  % N=148
% bandwidth_gof=goodnessOfFit(X,bandwidth,Y2,1);
% bar(1:(length(X')+1),[Y1/sum(Y1)*100,Y2/sum(Y2)*100],'histc');
% set(gca, 'XLim', [0 length(X')+2], 'XTickLabel',{num2str(X'); '>2.6'}); 
% title(sprintf('(c) Spatial frequency bandwidth (octave) [h=%1.3f]',bandwidth_gof));
% ylabel('occurrence (%)');


% Bandwidth vs Frequency

subplot(3,2,5);
X=[3.0 3.5 4.2 5.0 6.0 7.0];
freq2=freq*freq_scale;
for I=1:length(X)-1
    f=bandwidth(~isnan(bandwidth) & freq2>=X(I) & freq2<X(I+1));
    m(I)=mean(f);
    s(I)=std(f);
end;
errorbar(1:length(X)-1,m,s,'k');
set(gca, 'XLim', [0 length(X')], 'XTick', (0:length(X)), 'XTickLabel', {''; X(1:length(X)-1)'; ''});
title('(e) Bandwidth vs spatial frequency');
ylabel('Bandwidth (octave)');
xlabel('Frequency (cy/deg)');
% semilogx(freq*freq_scale,bandwidth,'bo');
% title('Bandwidth vs spatial frequency');
% ylabel('Bandwidth (octave)');
% xlabel('Frequency (cy/deg)');
% set(gca, 'XTick', [2.0 3.0 4.0 5.0 6.0]);
% xlim([3.0 7.0]);

% Apsect ratio

subplot(3,2,4);
if exist('db')==1
    X=0:0.25:4;
    exper_aspect=(db.ny./db.nx)';
else
    X=0:0.1:1;
    exper_aspect=[72 123 127 127 144 190 141 171 195 149 162 198 219 154 166 177 205 224 252 262 195 226 238 253 288]' / 293;
        % Palmer 1987
    aspect3(aspect>1)=1./(aspect(aspect>1));
end;
Y1=histc(aspect3,X);
Y2=histc(exper_aspect,X);
aspect_gof=goodnessOfFit([X,Inf],aspect3,[Y2;0],1);
bar(X,[Y1/sum(Y1)*100,Y2/sum(Y2)*100]);
set(gca, 'XLim', [0 max(X)]);
title(sprintf('(d) Apsect ratio [h=%1.3f]',aspect_gof));
ylabel('occurrence (%)');

% subplot(3,3,4);
% aspect_scale = 3.0;
% X=0:1:16;
% %Y1=histc(aspect2,X);
% Y1=histc(aspect * aspect_scale,X);
% Y2=[0 3 14 11 12 7 2 2 2 0 1 1 0 0 0 0 1]'; % Parker&Hawken 1988
% bar(X,[Y1/sum(Y1)*100,Y2/sum(Y2)*100],'histc');
% set(gca, 'XLim', [-0.5 16.5]);
% title('Aspect ratio');
% ylabel('occurrence (%)');

% Length

subplot(3,2,2);
X=0:5:40;
Y2=[3 13 11 11 6 6 3 2 1 0]';   % Parker&Hawken 1988
length_scale=60/freq_scale;
length_gof=goodnessOfFit([X,Inf],2*sigmay,Y2,length_scale);
Y1=histc(2*sigmay * length_scale,[X,Inf]);
bar(1:(length(X')+1),[Y1/sum(Y1)*100,Y2/sum(Y2)*100],'histc');
set(gca, 'XLim', [0 length(X')+2], 'XTickLabel',{num2str(X'); '>40'}); 
title(sprintf('(b) Length (arcmin) [scale=%3.1f h=%1.3f]',length_scale, length_gof));
ylabel('occurrence (%)');
[length_scale2,length_gof2]=scaleFit([X,Inf],2*sigmay,Y2);
fprintf('If lengths data are fit, scale=%3.1f h=%1.3f\n', length_scale2, length_gof2);

% Orientation

subplot(3,2,6);
hist(theta,32);
set(gca, 'XLim', [0 pi]);
title('Orientation');
ylabel('occurrence (%)');

    function [scale,gof] = scaleFit(X,data,Y2)
        obj=@(s)(goodnessOfFit(X,data,Y2,s));
        opt=optimset('Algorithm','interior-point');
        [scale,gof]=fminbnd(obj,min(X)/max(data),max(find(X<Inf))/min(data));
    end;
    
    function g = goodnessOfFit(X,data,Y2,s)
        Y=histc(data*s,X);
        N=sum(Y);
        N2=sum(Y2);
        g=hellingerdist(Y/N,Y2/N2);
        %g=rms(Y/N - Y2/N2);
        %g=1-rsquare(Y2/N2,Y/N);
    end;
    
end


