function net=v1_gratinganalyze(net,varargin)
% data=v1_gratinganalyze(data)
%  data : network data (in/out)
%  min_resp : minimum responses of units to consider
% analyzes f1/f0 ratios and optimal orientation/wavelength
% assumes unit props3 "grating"
% writes out unit props "f1f0ratio", "optorient", and "optwlen"

pr=inputParser;
pr.addParamValue('min_resp',0.33,@isnumeric);
pr.addParamValue('layers',2:length(net.structure.layers),@isnumeric);
pr.parse(varargin{:});
options=pr.Results;

for L=options.layers
    fprintf('layer %d\n',L);

    resp_all=net.content.layers{L}.unitProperties3.grating;
    [num_phases,num_freqs,num_orients,num_units,num_nodes]=size(resp_all);
    orients=net.content.layers{L}.layerProperties.orients(:);
    phases=net.content.layers{L}.layerProperties.phases(:);
    freqs=net.content.layers{L}.layerProperties.freqs(:);
    
    f1f0ratio=zeros(1,num_units,num_nodes)*NaN;
    optorient=zeros(1,num_units,num_nodes)*NaN;
    orientbw=zeros(1,num_units,num_nodes)*NaN;
    optfreq=zeros(1,num_units,num_nodes)*NaN;
    freqbw=zeros(1,num_units,num_nodes)*NaN;
    for N=1:num_nodes
        fprintf('%d ',N);
        parfor U=1:num_units
            resp=resp_all(:,:,:,U,N);
            resp=resp/max(resp(:));
            resp1=permute(mean(resp),[2 3 1]);
            [~,mi]=max(resp1(:));
            [optla,optt]=ind2sub([num_freqs num_orients],mi);
            resp_phase=resp(:,optla,optt);
            resp_freq=resp1(:,optt);
            resp_orient=permute(resp1(optla,:),[2 1]);
            amp=abs(fft(resp_phase));
            f1f0=amp(2)*2/amp(1);
            [optori,oribw,gof_orient]=orientation_tuning(orients,resp_orient);
            [optfr,frbw,gof_freq]=frequency_tuning(freqs,resp_freq);
            if max(resp_phase)>=options.min_resp && gof_orient>0.5 && gof_freq>0.5
                f1f0ratio(1,U,N)=f1f0;
                optorient(1,U,N)=optori;
                orientbw(1,U,N)=oribw;
                optfreq(1,U,N)=optfr; 
                freqbw(1,U,N)=frbw;
            end;            
        end;
        fprintf('\n');
    end;
    net.content.layers{L}.unitProperties.f1f0ratio=f1f0ratio;
    net.content.layers{L}.unitProperties.optorient=optorient;
    net.content.layers{L}.unitProperties.orientbw=orientbw;
    net.content.layers{L}.unitProperties.optfreq=optfreq;
    net.content.layers{L}.unitProperties.freqbw=freqbw;
    net.content.layers{L}.unitProperties.optwlen=1./optfreq;
end;

end

function [optorient,orientbw,gof]=orientation_tuning(orients,resp)
% half-maximal full-width

% resp=resp-min(resp);
[copt,gof]=vonmisesFit(resp);
optorient=copt.s1;
orientbw=acos(1-log(2)/copt.m1);
gof=gof.rsquare;

end

function [optfreq,freqbw,gof]=frequency_tuning(freqs,resp)
% half-maximal full-width

[copt,gof]=gaussFit(freqs,resp);
optfreq=copt.c;
radius=copt.sig*sqrt(2*log(2));
freqbw=log2((optfreq+radius)/(optfreq-radius));
gof=gof.rsquare;

end
