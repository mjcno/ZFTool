% Description: The program computes the minimum threshold in order that, 
% when computing cancer cells area and intensity, the autofluorescence of 
% the fish is not computed
% Final output: images from each zebrafish at 0h and 24-48-72h, with contour
% surrounding cells.
% The criterium in order to fix the threshold is explained in the paper.
%
% Data to adjust are:
% -Name of images (tiff images)
% -correc : Default value is 1 (it must be between 0.5 and 2). This is only
% used when fluorescence of cancer cells is too low or too high, which is
% not the normal case.
% -threshold (initialized to 0 or 5, until 40, in order to compute the
% evolution of fluorescence area with threshold)


%Clean all variables and all images
clear all
close all

fprintf('Image format at  0hpi for fish 2 is: 2 - 0h bn.tif and 2 - 0h gfp.tif\n');
fprintf('Image format at 48hpi for fish 2 is: 2 - 48h bn.tif and 2 - 48h gfp.tif\n');
% Reading of the input images (all in same folder)
inputfolder=input('Folder of input images: ', 's');
fish=input('Fish number: ');
%outputfolder=input('Introduzca la carpeta donde guardará los resultados: ', 's');

%Initial and final measurement hour
hour_in=input('Initial measurement hour (0hpi): ');
hour_end=input('Final measurement hour (24hpi, 48hpi or 72hpi): ');

%Reading of the 4 input images, two greys and two fluorescein GFP
%Transformation of input images to 0-255 range
image_bw_0=sprintf('%s/%d - %dh bn.tif',inputfolder,fish,hour_in);
image_gfp_0=sprintf('%s/%d - %dh gfp.tif',inputfolder,fish,hour_in);
image_bw_1=sprintf('%s/%d - %dh bn.tif',inputfolder,fish,hour_end);
image_gfp_1=sprintf('%s/%d - %dh gfp.tif',inputfolder,fish,hour_end);
%Transformation of input images to 0-255 range
f_bw_0=imread(image_bw_0); f_bw_0=uint8(double(f_bw_0)/65535*255);
f_gfp_0=imread(image_gfp_0); f_gfp_0=uint8(double(f_gfp_0)/65535*255);
f_bw_1=imread(image_bw_1);f_bw_1=uint8(double(f_bw_1)/65535*255);
f_gfp_1=imread(image_gfp_1);f_gfp_1=uint8(double(f_gfp_1)/65535*255);

% Correction factor between 0.5 and 2 which will allow to correct GFP image
% when it comes too weak or too strong and it can be confused with fish
% fluorescence. If 1, no changes are made to image, which is the normal
% case
correc=1;
f_gfp_0_c=correc*f_gfp_0;
f_gfp_1_c=correc*f_gfp_1;

% In this matrix, thresholds for each hour will be stored
thresholds=zeros(1,2);

nGFP_0=zeros(1,50/5+1);%number of pixels (areas) for 0hpi
nGFP_1=zeros(1,50/5+1);%number of pixels (areas) for 24-48-72hpi
meanGFP_0=zeros(1,50/5+1);%mean intensities for 0hpi
meanGFP_1=zeros(1,50/5+1);%mean intensities for 24-48-72 hpi
ndata=0; %number of data

%Measurement of area evolution with different thresholds
for threshold=0:5:50 
    %Computing number of points (areas) over each threshold, for images at
    %hour_in and hour_end, thresholding them in two vectors
    maskG_0=im2bw(f_gfp_0_c,threshold/255);
    maskG_1=im2bw(f_gfp_1_c,threshold/255);
    ndata=ndata+1;
    nGFP_0(ndata)=sum(maskG_0(:)==1);
    nGFP_1(ndata)=sum(maskG_1(:)==1);
    fGmascara_0=uint8(maskG_0).*f_gfp_0_c;
    meanGFP_0(ndata)=sum(fGmascara_0(:))/nGFP_0(ndata);
    fGmascara_1=uint8(maskG_1).*f_gfp_1_c;
    meanGFP_1(ndata)=sum(fGmascara_1(:))/nGFP_1(ndata);
    %Contours for image 0hpi over threshold
    perim=bwperim(fGmascara_0,8);
    %Contours overlaid on GFP image for 0hpi
    fcolor_0=cat(3,255*uint8(perim),f_gfp_0_c,0*f_gfp_0_c);
    title_0=sprintf('Perimeter GFP>%d at 0hpi',threshold);
    figure('Name',title_0), imshow(fcolor_0,'InitialMagnification','fit')
    %Contours for image 24-48-72hpi over threshold
    perim=bwperim(fGmascara_1,8);
    %Contours overlaid on GFP image for 24-58-72hpi
    fcolor_1=cat(3,255*uint8(perim),f_gfp_1_c,0*f_gfp_1_c);
    title_1=sprintf('Perimeter GFP>%d at %dhpi',threshold,hour_end);
    figure('Name',title_1), imshow(fcolor_1,'InitialMagnification','fit')
end

%Graphical representation of evolution of cancer mass area and intensity
%when increasing GFP threshold
figure, plot(0:5:50,nGFP_0,'b*-',0:5:50,nGFP_1,'ro-'),grid on
ylabel('#pixels with GFP > threshold')
xlabel('GFP thresholds')
title_area=sprintf('Evolution of #pixels with GFP threshold (B=0hpi; R=%dhpi)',hour_end);
title(title_area)
figure, plot(0:5:50,meanGFP_0,'b*-',0:5:50,meanGFP_1,'ro-'), grid on
ylabel('Mean intensity of pixels with GFP > threshold')
xlabel('GFP thresholds')
title_intensity=sprintf('Evolution of mean intensity with GFP threshold (B=0hpi; R=%dhpi)',hour_end);
title(title_intensity)

%Computation of the common automatic threshold for the fish being analyzed
%Analysis of nGFP_0 and nGFP_1 vectors in order to determine a common threshold

%Detection of a high decay in data from one component to the next in twho 
%vectors at a time, so we compare each value with its precedent
comparation_0=zeros(1,ndata-1);
comparation_1=zeros(1,ndata-1);
for k=1:ndata-1
    comparation_0(k)=nGFP_0(k+1)/nGFP_0(k);
    comparation_1(k)=nGFP_1(k+1)/nGFP_1(k);
end

%Conditions: threshold must be the same through 3 iterations
%Parameters: ratio of change between one threshold and next
difthresh=1; %Initialization.
%It means that thresholds are exactly the same 
%It will decrease in 0.005 intervals until algorithm converges

%Initialization: threshold must be the same (difumb) through 3 iterations. 
%If not, decrease difumb by 0.005 and try again.
found=zeros(1,2);
while ~all(found)
    found=zeros(1,2);
    difthresh=difthresh-0.005; %decrease threshold
    counter_0=0;counter_1=0;
    for k=1:ndata-1
        if comparation_0(k)>difthresh
            counter_0=counter_0+1;
            if counter_0==3
                %threshold for 0hpi is in first position
                thresholds(1)=(k-2)*5;
                found(1)=1;
            end
        else
            counter_0=0;
        end
        if comparation_1(k)>difthresh
            counter_1=counter_1+1;
            if counter_1==3
                %threshold for 24-48-72hpi is in first position
                thresholds(2)=(k-2)*5; 
                found(2)=1;
            end
        else
            counter_1=0;
        end
        
    end
end
%Output
%Final value of ratio of thresholds
fprintf('%% of change between thresholds: %.0f%%\n',difthresh*100);
fprintf('Final common threshold = %d (maximum of [%d,%d])\n',max(thresholds),thresholds)

%Output composed images
%Threshold 0 and final threshold with GFP over grey image
%Para paper: composición de imagen de pez 14 con umbral 0 y con umbral
%final de 35
finalthreshold=max(thresholds);
%Threshold 0
threshold=0;
maskG_0h_t0=im2bw(f_gfp_0_c,threshold/255);
maskG_1h_t0=im2bw(f_gfp_1_c,threshold/255);
perimR_0h=bwperim(maskG_0h_t0,8); %threshold 0 perimeter
%Final threhold
threshold=finalthreshold;
maskG_0h_tf=im2bw(f_gfp_0_c,threshold/255);
maskG_1h_tf=im2bw(f_gfp_1_c,threshold/255);
perimR_1h=bwperim(maskG_1h_t0,8); %threshold 1 perimeter

perimB_0h=bwperim(maskG_0h_t0-maskG_0h_tf);
perimB_1h=bwperim(maskG_1h_t0-maskG_1h_tf);

%Composition of final images
%Image at 0hpi showing threshold=0 and threshold=final
fcolor_0h=cat(3,f_bw_0+255*uint8(perimR_0h),f_bw_0+f_gfp_0_c,f_bw_0+255*uint8(perimB_0h));
figure('Name','0hpi: threshold=0(pink) - final threshold(blue)'),imshow(fcolor_0h);
%Image at 24-48-72hpi showing threshold=0 and threshold=final
fcolor_1h=cat(3,f_bw_1+255*uint8(perimR_1h),f_bw_1+f_gfp_1_c,f_bw_1+255*uint8(perimB_1h));
figure('Name','24-48-72hpi: threshold=0(pink) - final threshold(blue)'),imshow(fcolor_1h);
