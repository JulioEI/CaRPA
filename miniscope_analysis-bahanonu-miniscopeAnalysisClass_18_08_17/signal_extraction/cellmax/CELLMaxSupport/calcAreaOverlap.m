function areaOverlap = calcAreaOverlap(image1, image2, varargin)

options.fractionMaxThresh=0.15;
options=getOptions(options, varargin);

image1(image1<options.fractionMaxThresh*max(image1(:)))=0;
image1(image1>0)=1;
image2(image2<options.fractionMaxThresh*max(image2(:)))=0;
image2(image2>0)=1;
areaOverlap=sum(and(image1(:),image2(:)))/sum(sum(image1(:))+sum(image2(:)));