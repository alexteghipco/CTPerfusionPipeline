function fnm = convert_ctp(fnm, intensityMax, isTTD, mm, clip)
%Convert Siemens RGB CT-perfusion images to linear grayscale images
% fnm : name of image[s] to convert
% intensityMax : brightest value in image
% isTTD : Siemens uses different RGB schemes 
%              TRUE for TTD, MTT or TTP images (time to drain, transit time, time to peak)
%              FALSE for CBF or CBV images (cerebral blow flow, volume)
% mm : Row, column and slice spacing in mm
% clip : voxels to clip on left and right side of image (remove colorbar)
%Examples
% convert_ctp(); %use grphical interface
% convert_ctp('208.nii', 100, false, [0.5 0.5 4.0], [64 80])
% convert_ctp('206.nii', 15, true, [0.5 0.5 4.0], [64 80])

if ~exist('fnm','var')
    fnm = spm_select(1,'image','Select RGB image'); 
    def = {'100','1','0.5 0.5 4.0', '64 80'};
    prompt = {'Intensity Maximum','TimeToDrain (0) or CBF/CBV (1)','Spatial Resolution (mm LR/AP/HF)', 'Clip Voxels L/R'};
    answer = inputdlg(prompt, 'rgb2scalar settings', 1, def);
    intensityMax = str2double(answer(1));
    isTTD = (str2double(answer(2)) == 0);
    mm = str2num(answer{3});
    clip = str2num(answer{4});
end

if ~exist('fnm','var') || isempty(fnm), return; end;
%if ~exist('fnm','var'), fnm = spm_select(1,'image','Select image to convert from RGB to scalar'); end
if ~exist('isTTD','var'),isTTD = false; end;
%fprintf('rgb->scalar %s\n',fnm);
fnm = deblank(fnm);
hdr = spm_vol(fnm);
if hdr.dt(1) ~= 128
    fprintf('%s error: input file must be RGB format %s\n',mfilename,fnm);
    return;
end
hdrM = hdr;
hdrM.dt(1)    =2; %2=8-bit char, 4= 16-bit integer; 16 =32-bit real datatype
hdrM.pinfo(3) = 352;
isPlanar = false; %is RGB voxels saved as triplets (RGBRGB...RGB or planar RR..RGG..GBB..B
if isPlanar
    hdrM.dim(3) = hdrM.dim(3) * 3; %red,green,blue each saved on separate slice
    imgRGB = spm_read_vols(hdrM);
    imgR = imgRGB(:,:,(1 : 3: hdrM.dim(3))); %red plane
    imgG = imgRGB(:,:,(2 : 3: hdrM.dim(3))); %green plane
    imgB = imgRGB(:,:,(3 : 3: hdrM.dim(3))); %blue plane
else
    hdrM.dim(1) = hdrM.dim(1) * 3; %red,green,blue each saved on separate slice
    imgRGB = spm_read_vols(hdrM);
    imgR = imgRGB((1 : 3: hdrM.dim(1)),:,:); %red plane
    imgG = imgRGB((2 : 3: hdrM.dim(1)),:,:); %red plane
    imgB = imgRGB((3 : 3: hdrM.dim(1)),:,:); %red plane
end
if isTTD
    img = ttd_rgb2scalarSub(imgR(:), imgG(:), imgB(:));
    %img = ttp_rgb2scalarSub(imgR(:)', imgG(:)', imgB(:)');
else
    img = cbfd_rgb2scalarSub(imgR(:)', imgG(:)', imgB(:)');
    %img = img * 10; %we are saving as integers, so preserve precision
end
mx = max(img);
img = img * intensityMax / mx;
img = reshape(img,hdr.dim);
img = flip(img,1);
hdrOut = hdr;
%clip Left/Right
img = img(1+clip(2):size(img,1)-clip(1), :, :, :); %clip image dimensions
hdrOut.dim = size(img);
if ~isempty(mm) 
    % set origin and scale
    hdrOut.mat = [-mm(1) 0 0 0; 0 mm(2) 0 0; 0 0 mm(3) 0; 0 0 0 1];
    %origin is center of volume
    vx = [-hdrOut.dim(1:3), 2] * 0.5;
    vx = vx * hdrOut.mat;
    hdrOut.mat(:,4) = vx;
end
[pth, nam, ext] = spm_fileparts(fnm);
hdrOut.fname = fullfile(pth, ['q' nam ext]);
hdrOut.dt(1)    = 16; %2=8-bit char, 4= 16-bit integer; 16 =32-bit real datatype
hdrOut.pinfo(3) = 352;
spm_write_vol(hdrOut,img);
fnm = hdrOut.fname;
%end convert_ctpSub()

function s = cbfd_rgb2scalarSub (R,G,B)
%converts Siemens RGB color scheme for CABFD to scalar intensity
%fprintf('CBFD %g..%g\n', min(R(:)), max(R(:)) );
s = zeros(numel(R),1);
%segment 2: indices 23..42
idx = intersect (find (R > G), find(B > G ) );
s(idx) =  ((-R(idx)+B(idx)+G(idx)-4) * (19/122)) + 23;
%segment 1: indices 1..22
idx = find (G <= 1); %1..23 red ramps 0..255
s(idx) = (R(idx)+B(idx)-122) * (22/68);
%segment 3: indices 42..78
idx = find(G == R) ;
s(idx) = ((B(idx) - 130) * (35/123)) + 43;
%segment 4 79..158
idx = intersect (find (G > R), find(B > 0 ) );
s(idx) = ((R(idx) + G(idx) - B(idx) +124) * (79/503)) + 79;
%segment 5 159..229
idx = find (B < 1);
s(idx) = ((R(idx) - 128) * (70/126)) + 159;
%segment 6 230..254
idx = intersect (find (R > G), find(R > B ) );
s(idx) = ((-R(idx) -G(idx) +B(idx) + 495) * (24/270)) + 230;
%set black
s(((R+G+B) < 1)) = 0;
s(s < 0) = 0;
%fprintf('CBFD %g..%g\n', min(s(:)), max(s(:)) );
%end cbfd_rgb2scalarSub()

function s = ttd_rgb2scalarSub (R,G,B)
%converts Siemens RGB color scheme for MTT/TTP to scalar intensity
%fprintf('%g..%g\n', min(R(:)), max(R(:)) );
s = zeros(numel(R),1);
%segment 1: blue > 0 indices 1..63
idx = intersect (find (G <= 64), find(B > 0 ) ); %0:64 gren ramps 0..64
s(idx) = (-R(idx)+G(idx)+B(idx)+4) * 0.245;
%segment 2: blue > 0 indices 64..127
idx = intersect (find (G > 64), find(B > 0 ) ); 
s(idx) = ((-R(idx)+G(idx)-B(idx)) * 0.125786164) + 95.3;
%segment 3: blue = 0 indices 128..191
idx = intersect (find (G > 252), find(B == 0 ) ); 
s(idx) = (R(idx)* 0.252) + 127.75;
%segment 4: blue = 0 indices 192..254 (Green Decay)
idx = intersect (find (R > 252), find(B == 0 ) ); 
s(idx) = ((255-G(idx))* 0.247) + 191.75;
%set black
s(((R+G+B) < 1)) = 0;
s(s < 0) = 0;
%end ttd_rgb2scalarSub()
