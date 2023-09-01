function ct_rgb(fnms,varargin)
%Process CT perfusions scans
% fnms : file name[s] of CT image[s] (optional)
% 
%Optional arguments (added by A.T.) 
% 1) add 'prebet' after fnms to stop the pipeline at brain extraction (the
% z image will be thresholded, this was originally performed during the bet
% step).
%
% 2) add 'betted' and ensure first filepath in fnms is a betted image to
% resume pipleine after bet.
%
% 3) add 'normalize' and optionally provide a smoothing value (mm) for the
% reference, then source image (e.g., 'normalize',0,8). The values in the
% example are defaults that are used by the pipeline. You can supply only a
% smoothing value for the reference (e.g., 'normalize',0) if you want. This
% will start pipeline at normalization. 
%
% Notes: Thresholding step was made into its own function to threshold brain
% without betting - A.T.
%
%Examples
% ct_rgb({'CTA_Head_Neck_4.nii', 'RAPID_rCBF_332.nii', 'RAPID_rCBV_331.nii', 'RAPID_MTT_333.nii'});
% ct_rgb({'mean.nii', 'RAPID_rCBF_332.nii', 'RAPID_rCBV_331.nii', 'RAPID_MTT_333.nii'});
% ct_rgb({'mean.nii', 'RAPID_rCBF_332.nii'});
%
%Examples including optional arguments
% ct_rgb({'CTA_Head_Neck_4.nii', 'RAPID_rCBF_332.nii', 'RAPID_rCBV_331.nii', 'RAPID_MTT_333.nii'},'prebet');
% ct_rgb({'bzCTA_Head_Neck_4.nii', 'RAPID_rCBF_332.nii', 'RAPID_rCBV_331.nii', 'RAPID_MTT_333.nii'},'betted');
% ct_rgb({'bzCTA_Head_Neck_4.nii', 'RAPID_rCBF_332.nii', 'RAPID_rCBV_331.nii', 'RAPID_MTT_333.nii'},'normalize');
% ct_rgb({'bzCTA_Head_Neck_4.nii', 'RAPID_rCBF_332.nii', 'RAPID_rCBV_331.nii', 'RAPID_MTT_333.nii'},'normalize',0);
% ct_rgb({'bzCTA_Head_Neck_4.nii', 'RAPID_rCBF_332.nii', 'RAPID_rCBV_331.nii', 'RAPID_MTT_333.nii'},'normalize',0,8);

if ~exist('fnms','var')
   %return;
   [files,pth] = uigetfile({'*.nii;';'*.*'},'Select the CT images', 'MultiSelect', 'on'); 
   fnms = cellstr(files); %make cellstr regardless of whether user selects single or multiple images 
end
%check inputs:
if exist('spm','file') ~= 2, error('Please install SPM12 or later'); end;
% if numel(fnms) < 2
%     error('At least 2 images required: structural scan and RGB scan(s)');
% end

% AT--changes to take in arguments that stop pipeline before bet, start
% pipeline after bet, or start pipeline at normalization
betOpt = 0; % we assume you want default pipeline unless you passed in some extra arguments
normOpt = 0;
smor = 4; % smoothing in mm for reference
smos = 8; % smoothing in mm for source
if length(varargin) ~= 0
    if strcmpi(varargin{1},'prebet') % check that this is the argument passed in
        betOpt = 1; % this will stop pipeline prior to bet step
        disp('Running pipeline up until brain extraction...')
    elseif strcmpi(varargin{1},'betted')
        disp('You have specified a brain extracted image to use...searching for the image')
        if isfile(fnms{1}) % check that file being passed in exists
            disp('Image was found!')
            betOpt = 2; % this starts pipeline after bet
            % gunzip if necessary...
            if strcmpi(fnms{1}(end-2:end),'.gz')
                gunzip([fnms{1} '.gz'])
                fnms{1} = fnms{1}(1:end-3);
                disp('You passed in an .nii.gz image so it was decompressed')
            end
        elseif isfile([fnms{1} '.gz']) % maybe you forgot to unzip but assumed you did
            disp('Found a gzipped version of your image...decompressing now')
            gunzip([fnms{1} '.gz'])
            betOpt = 2;
        else
            error(['The image you provided does not exist (' fnms{1} '). Maybe you did not provide the full path?'])
        end
    elseif strcmpi(varargin{1},'normalize')
        normOpt = 1;
        if length(varargin) == 1
            smor = 4;
            smos = 8;
        elseif length(varargin) == 2
            smor = varargin{2};
            smos = 8;
        elseif length(varargin) == 3
            smor = varargin{2};
            smos = varargin{3};
        else
            error('Too many arguments have been passed in')
        end
    else
        error([varargin{1} ' is not a valid option'])
    end
end

% We get mean struct if we are running entire pipeline or pipeline up to
% bet
if (betOpt == 0 && normOpt == 0) || betOpt == 1
    %convert 4D -> 3D structural
    struct = meanStruct(fnms{1});
elseif betOpt == 2 || normOpt == 1
    struct = fnms{1};
end

% If we run pipeline up to bet we need to manually get z image which you
% can then bet yourself. If we are running pipeline after bet, we need to
% set origin then can resume pipeline as normal
if betOpt == 0 && normOpt == 0
    %brain extract CT structural scan
    struct = ct_bet(struct);
elseif betOpt == 1
    [struct,~] = thresh(struct);
elseif betOpt == 2
    set_origin(struct);    
end

if betOpt ~= 1
    %convert RGB data to scalar:
    if normOpt == 0
        scalars = rgb2scalar(fnms(2:end), struct);
    end
    atlases = normTemplate2CT(struct,smor,smos);
    scalar_atlas(fileparts(fnms{1}));
end

%scalars
%end ct_rgb()

%------- LOCAL functions follow  -----

function atlases = normTemplate2CT(struct,smor,smos)
%Scalp stripped stroke template https://github.com/neurolabusc/Clinical
% template = fullfile(spm('Dir'),'toolbox','Clinical','scct_stripped.nii');
template = fullfile(fileparts(which(mfilename)),'template','ct.nii');
if ~exist(template, 'file')
    error('Unable to find %s', template);
end
pth = fullfile(fileparts(which(mfilename)), 'atlas');
fnms = dir(fullfile(pth, '*.nii'));
if isempty(fnms)
    error('Unable to find files in %s', pth);
end
atlases = {};
struct_pth = fileparts(struct);
for i = 1:numel(fnms)
    src = fullfile(pth, fnms(i).name);
    dest = fullfile(struct_pth, fnms(i).name);
    copyfile(src,dest);
    atlases{end+1} = dest;
end
atlases = atlases';
%clinical_h2c(struct);
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = {template};
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = atlases;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {struct};
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = smos; % source image smoothing; default = 8; 16
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = smor; %4; % 0
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25; 
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16; 
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;  
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [nan nan nan; nan nan nan];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [nan nan nan];
%matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [-78 -112 -70; 78 76 85];
%matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [2 2 2];
%Use nearest neighbor interpolation: boundary of area 18 and 17 is NOT 17.5
%matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'w';
spm_jobman('run',matlabbatch);
%rename 'w' prefix to 'atlas':
for i = 1:numel(atlases)
    delete(atlases{i}); %unwarped
    src = fullfile(struct_pth, ['w', fnms(i).name]);
    dest = fullfile(struct_pth, ['atlas_', fnms(i).name]);
    movefile(src, dest);
    atlases{i} = dest;
end
%end meanStuct()

function fnm = meanStruct(fnm)
hdr = spm_vol(fnm);
if (hdr(1).dt(1) == 128)
	error('Scalar (gray-scale) image required: %s', hdr.fname);
end

if (numel(hdr) == 1)
    return; %3D dataset
end
matlabbatch{1}.spm.spatial.realign.estwrite.data = {{fnm}};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
[p,n,x] = fileparts(fnm);
fnm = fullfile(p, ['mean', n, x]);
%end meanStuct()

function fnms = rgb2scalar(fnmsRGB, struct)
%convert RGB images to scalar
fnms = {};
for i = 1:numel(fnmsRGB)
    fnm = fnmsRGB{i};
    hdr = spm_vol(fnm);
    if ((numel(hdr) > 1) || (hdr(1).dt(1) ~= 128))
        error('Not a single-volume RGB image: %s', hdr.fname);
    end
    %mm = ([0,0,0,1] * hdr.private.mat')-([1,1,1,1] * hdr.private.mat' );
    %we can handle TTD, MTT or TTP  CBF or CBV
    [p,n,x] = fileparts(fnm);
    %if ~isempty(strfind(upper(n),upper('MTT')))
    if contains(n,'MTT', 'IgnoreCase',true)
        fnms{end+1} = convert_ctp(fnm, 20, false, [], [32, 32]);
    elseif contains(n,'Tmax', 'IgnoreCase',true)
        fnms{end+1} = convert_ctp(fnm, 20, false, [], [32, 32]);
    elseif contains(n,'CBF', 'IgnoreCase',true)
        fnms{end+1} = convert_ctp(fnm, 20, false, [], [32, 32]);
    elseif contains(n,'CBV', 'IgnoreCase',true)
        fnms{end+1} = convert_ctp(fnm, 20, false, [], [32, 32]);
    else
        error('RGB filenames must contain TTD,Tmax,MTT,TTP,CBF or CBV: %s', fnm);
    end
    set_origin(fnms{end});
    fnms{end} = coregSub(fnms{end}, struct, 'scalar_');
end
%end rgb2scalar()

function fnm = coregSub(src, ref, prefix)
if ~exist('prefix','var')
    prefix = 'r';
end
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {ref};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {src};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 2;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = prefix;
spm_jobman('run',matlabbatch);
[p,n,x] = fileparts(src);
fnm = fullfile(p, [prefix, n, x]);
%end coreg()

function set_origin(fnm)
%CT images use table center (near abdomen) as origin
% SPM expects this near anterior commissure
hdr = spm_vol(fnm); %load header 
img = spm_read_vols(hdr); %load image data
img = img - min(img(:));
img(isnan(img)) = 0;
%find center of mass in each dimension (total mass divided by weighted location of mass
% img = [1 2 1; 3 4 3];
sumTotal = sum(img(:));
coivox = ones(4,1);
coivox(1) = sum(sum(sum(img,3),2)'.*(1:size(img,1)))/sumTotal; %dimension 1
coivox(2) = sum(sum(sum(img,3),1).*(1:size(img,2)))/sumTotal; %dimension 2
coivox(3) = sum(squeeze(sum(sum(img,2),1))'.*(1:size(img,3)))/sumTotal; %dimension 3
XYZ_mm = hdr.mat * coivox; %convert from voxels to millimeters
hdr.mat(1,4) =  hdr.mat(1,4) - XYZ_mm(1);
hdr.mat(2,4) =  hdr.mat(2,4) - XYZ_mm(2);
hdr.mat(3,4) =  hdr.mat(3,4) - XYZ_mm(3);
spm_create_vol(hdr);
%end set_origin()

function [fnm,hdr] = thresh(fnm)
hdr = spm_vol(fnm);
img = spm_read_vols(hdr);
offset = 40;
img = img + offset;
img(img < 0) = 0;
img(img > (100+offset)) = 0;
% img(img < 0) = 0;
% img(img > 100) = 0;
[pth, nm, ext] = spm_fileparts(fnm);
hdr.fname = fullfile(pth, ['z' nm ext]);  
spm_write_vol(hdr,img);

function fnm = ct_bet(fnm)
%https://pubmed.ncbi.nlm.nih.gov/25862260/
% Threshold 0-100HU; in the other, data were not smoothed. 
% Blur Gaussian kernel (σ=1mm(3)) 
% BET (FI) thresholds: 0.01 or 0.1
%for alternative https://git.fmrib.ox.ac.uk/thanayik/ct_bet 
fsldir= '/usr/local/fsl';
if ~exist(fsldir,'dir')
    fsldir='/Users/chrisrorden/fsl';
end
if ~exist(fsldir,'dir')
    error('Unable to find %s', fsldir);
end
[fnm,hdr] = thresh(fnm); % pulled this code out into a function so we can have a prebet pipeline
betFI = 0.1;
fnm = betSub(fsldir,hdr.fname, betFI);
set_origin(fnm);
%end betSub()

function maskNam = betSub(fsldir,imgNam, betFI) %brain extract
setenv('FSLDIR', fsldir);
setenv('FSLOUTPUTTYPE', 'NIFTI');
[pth,nam,ext] = fileparts(imgNam);
if strcmpi(ext,'.gz'), [~,nam] = fileparts(nam); end;
inNam = imgNam;
maskNam = fullfile(pth, ['b', nam, ext] ); %will generate image "dti_mask.nii.gz"
if (contains(maskNam, ' ')) || (contains(inNam, ' '))
    error('FSL does not allow spaces in file paths: "%s"', inNam)
end
%command=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/bet %s %s -f %g -g 0"\n',inNam,maskNam, betFI);
command=sprintf('sh -c ". ${FSLDIR}/bin/bet "%s" "%s" -f %g -g 0"\n',inNam,maskNam, betFI);
system(command);
fprintf(command);
%end betSub()

