surrogatePermuted = [3.8797    3.4820    3.5412    3.4967    3.6263    3.4914    3.4754    3.7460    3.5838    3.4548];
surrogateNonPermuted = [8.3963    9.5479    8.5000    7.9043    8.2620    8.0751    8.4548    8.4608    8.8065    8.5346];
calciumPermutePosition = [8.9176    8.5898    8.4096    8.7094    8.8012    8.4761    8.6988    8.6935    8.5499    8.3517];
calciumNonPermute = [0.4608    0.3524    0.3783    0.4348    0.6037    1.5997    1.0199    0.6882    0.5313    0.7227];
calciumPermute= [0.2094    0.1855    0.2114    0.2154    0.2094    0.1948    0.2068    0.2081    0.2061    0.1955];


figure;
boxplot([surrogatePermuted;surrogateNonPermuted;calciumPermutePosition;calciumNonPermute;calciumPermute]',{'Surrogate permuted','Surrogate non permuted','Position permuted','Non permuted','Permuted'})
ylabel('Mean absolute decoding error')
set(gcf,'color','white')
set(gcf,'units','normalized','outerposition',[0 0 1 1])
%axis square;
grid on;
%export_fig asdf.jpg -m3
%annotation('textbox',[0,0,1,1],'string',['Decoding by kfold proximity perm (',param.traceField,')'],'fontSize',20,'fontWeight','bold')
%save2pdf(['C:\Users\csc\Desktop\caImagg\Graphics\closeKDecodingPerm'])