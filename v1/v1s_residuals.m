function v1s_residuals(net)
% shows histogam of Gabor fitting residuals

res = getResiduals(net);
figure;

X=0:0.05:0.5;
Y=histc(res,[X,Inf]);
bar(1:(length(X)+1),Y/length(res)*100,'histc');
set(gca, 'XLim', [0 length(X)+1], 'XTick', 1:(length(X)+1), 'XTickLabel',{num2str(X');''}); 
title('residuals');
ylabel('occurrence (%)');

fprintf('bad fitting (residual>0.1): %d\n', sum(res>0.1));

end

