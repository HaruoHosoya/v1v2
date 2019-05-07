function compute_denoising_variances

lams=0.1:0.1:3;

sigma2=0.2;

T=50;
lamc=zeros(T,length(lams));

for I=1:length(lams)
    l=10;
    for J=1:T
        l=l^2/(l+sigma2)^2*lams(I);
        lamc(J,I)=l;
    end;
end;

figure;
plot(lamc);

end
