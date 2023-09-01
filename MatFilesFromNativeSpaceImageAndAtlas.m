addpath(genpath('/Users/alex/Documents/CSTAR/scripts/AphasiaResearchCohortQuery/AcuteStroke/ProcessingFiles'))
cd('/Volumes/Quattro/ct_alex/acute2/rerun3')
filepaths = dir(pwd)
filepaths = filepaths(1:length(filepaths))
roipath = '/Users/alex/Documents/CSTAR/scripts/NiiStat-master/roi';
roipath_extra = '/Users/alex/Documents/CSTAR/scripts/NiiStat-master/roi_uiowa';
outputPath = '/Volumes/Quattro/ct_alex/acute2/matFiles';

mkdir(outputPath);
MTTdir = [outputPath filesep 'MTT']; mkdir(MTTdir);
rCBFdir = [outputPath filesep 'rCBF']; mkdir(rCBFdir);
rCBVdir = [outputPath filesep 'rCBV']; mkdir(rCBVdir);
Tmaxdir = [outputPath filesep 'Tmax']; mkdir(Tmaxdir);

for i = 1:length(filepaths) %1:7%length(filepaths)
    disp(['working on sub ' num2str(i) ' of ' num2str(length(filepaths))])
    try
        cd ([filepaths(i).folder '/' filepaths(i).name])
       
        MTT = dir('q*MTT*'); MTTFile = [[MTT(1).folder '/' MTT(1).name]];
        rCBF = dir('q*rCBF*');rCBFFile = [[rCBF(1).folder '/' rCBF(1).name]];
        rCBV = dir('q*rCBV*');rCBVFile = [[rCBV(1).folder '/' rCBV(1).name]];
        Tmax = dir('q*Tmax*');TmaxFile = [[Tmax(1).folder '/' Tmax(1).name]];
        nii_reslice_target('T1_LM1001.nii','','jhu.nii',false);
        aalInfo = dir('atlas_aal.nii'); aalFile = [[aalInfo(1).folder '/' aalInfo(1).name]];
        jhuInfo = dir('atlas_jhu.nii'); jhuFile = [[jhuInfo(1).folder '/' jhuInfo(1).name]];

        nii_reslice_target(MTTFile,'', aalFile,false)
        nii_reslice_target(rCBFFile,'', aalFile,false)
        nii_reslice_target(rCBVFile,'', aalFile,false)
        nii_reslice_target(TmaxFile,'', aalFile,false)

        %COPY aal and jhu atlases to NiiStat roi directories to ensure they use
        %subject specific native space image!!! Need to replace old files in
        %normal space afterwards!!
        copyfile(aalFile, [roipath '/' 'aal.nii']);
        copyfile(jhuFile, [roipath '/' 'jhu.nii']);
        copyfile(aalFile, [roipath_extra '/' 'aal.nii']);
        copyfile(jhuFile, [roipath_extra '/' 'jhu.nii']);

        %FINAL RESLICED FILES TO WORK WITH!!!
        MTTfnm = strrep(MTTFile,'qR','rqR');
        rCBFfnm = strrep(rCBFFile,'qR','rqR');
        rCBVfnm = strrep(rCBVFile,'qR','rqR');
        Tmaxfnm = strrep(TmaxFile,'qR','rqR');
        
        nii_roi2stats('aal',MTTfnm,'','cbf_',[MTTdir filesep filepaths(i).name '_MTT.mat']);
        nii_roi2stats('jhu',MTTfnm,'','cbf_',[MTTdir filesep filepaths(i).name '_MTT.mat']);
        
        nii_roi2stats('aal',rCBFfnm,'','cbf_',[rCBFdir filesep filepaths(i).name '_rCBF.mat']);
        nii_roi2stats('jhu',rCBFfnm,'','cbf_',[rCBFdir filesep filepaths(i).name '_rCBF.mat']);
        
        nii_roi2stats('aal',rCBVfnm,'','cbf_',[rCBVdir filesep filepaths(i).name '_rCBV.mat']);
        nii_roi2stats('jhu',rCBVfnm,'','cbf_',[rCBVdir filesep filepaths(i).name '_rCBV.mat']);
        
        nii_roi2stats('aal',Tmaxfnm,'','cbf_',[Tmaxdir filesep filepaths(i).name '_Tmax.mat']);
        nii_roi2stats('jhu',Tmaxfnm,'','cbf_',[Tmaxdir filesep filepaths(i).name '_Tmax.mat']);
    catch
        fileID = fopen([outputPath '/errorLog.txt'],'a');
        fprintf(fileID,'error on particitpant %s \n',filepaths(i).name);
        fclose(fileID);
    end
end
    
% create global mean scaled and global scaled images to mat file. You will need load_nifti.m and save_nifti.m from freesurfer for this code but replace with whatever you use to load/save these files in matlab (e.g., readnifti.m writenifti.m now come with the image processing toolbox in matlab)
cd('/Volumes/Quattro/ct_alex/acute2/rerun')
subs = dir(pwd);
types = {'globalScaled','globalScaledRH', 'globalMean', 'globalMeanRH'};
an = readtable('/Users/alex/Documents/CSTAR/scripts/NiiStat-master/roi/backup/jhu.txt');
anrid = find(contains(an.Var2,'_R'));
for i = 1:length(subs)
    disp(['working on sub ' num2str(i) ' of ' num2str(length(subs))])
    cd([subs(i).folder '/' subs(i).name])
    fs = dir('rq*.nii*');
    a = load_nifti('atlas_jhu.nii');
    aid = find(a.vol ~= 0); % non-zero voxels in image
    rid = find(ismember(a.vol(:),anrid)); % voxels that are in RH
    for j = 1:length(fs)
        f = load_nifti(fs(j).name);
        fo = f;
        fo.vol = zeros(size(fo.vol));
        [pth,nm,ext] = fileparts(fs(j).name);
        for k = 1:length(types)
            switch types{k}
                case 'globalScaled'
                    sF = 100/mean(f.vol(aid));
                    b = f.vol(:).*sF;
                case 'globalScaledRH'
                    sF = 100/mean(f.vol(rid));
                    b = f.vol(:).*sF;
                case 'globalMean'
                    b = f.vol./mean(f.vol(aid));
                case 'globalMeanRH'
                    b = f.vol./mean(f.vol(rid));
            end
            fo.vol(:) = b;
            save_nifti(fo,[nm '_' types{k} ext]);
        end
    end
end

outputPath = '/Volumes/Quattro/ct_alex/acute2/matFiles_scaled';
types = {'globalScaled','globalScaledRH', 'globalMean', 'globalMeanRH'};
mkdir(outputPath);
MTTdir = [outputPath filesep 'MTT']; mkdir(MTTdir);
rCBFdir = [outputPath filesep 'rCBF']; mkdir(rCBFdir);
rCBVdir = [outputPath filesep 'rCBV']; mkdir(rCBVdir);
Tmaxdir = [outputPath filesep 'Tmax']; mkdir(Tmaxdir);
cd('/Volumes/Quattro/ct_alex/acute2/snrStragglers')
filepaths = dir(pwd);
filepaths = filepaths(1:length(filepaths));
for i = 1:length(filepaths) %1:7%length(filepaths)
    parfor j = 1:length(types)
        type = types{j}; %'globalScaledRH' 'globalMean' 'globalMeanRH'
        disp(['working on sub ' num2str(i) ' of ' num2str(length(filepaths))])
        disp(['working on type ' num2str(j) ' of ' num2str(length(types))])
        try
            cd ([filepaths(i).folder '/' filepaths(i).name])
            
            MTT = dir(['rq*MTT*' type '*']); MTTFile = [[MTT(1).folder '/' MTT(1).name]];
            rCBF = dir(['rq*rCBF*' type '*']);rCBFFile = [[rCBF(1).folder '/' rCBF(1).name]];
            rCBV = dir(['rq*rCBV*' type '*']);rCBVFile = [[rCBV(1).folder '/' rCBV(1).name]];
            Tmax = dir(['rq*Tmax*' type '*']);TmaxFile = [[Tmax(1).folder '/' Tmax(1).name]];
            nii_reslice_target('T1_LM1001.nii','','jhu.nii',false);
            aalInfo = dir('atlas_aal.nii'); aalFile = [[aalInfo(1).folder '/' aalInfo(1).name]];
            jhuInfo = dir('atlas_jhu.nii'); jhuFile = [[jhuInfo(1).folder '/' jhuInfo(1).name]];
            
            %COPY aal and jhu atlases to NiiStat roi directories to ensure they use
            %subject specific native space image!!! Need to replace old files in
            %normal space afterwards Roger!!
            copyfile(aalFile, [roipath '/' 'aal.nii']);
            copyfile(jhuFile, [roipath '/' 'jhu.nii']);
            copyfile(aalFile, [roipath_extra '/' 'aal.nii']);
            copyfile(jhuFile, [roipath_extra '/' 'jhu.nii']);
            
            %FINAL RESLICED FILES TO WORK WITH!!!
            MTTfnm = MTTFile; %strrep(MTTFile,'qR','rqR');
            rCBFfnm = rCBFFile; %strrep(rCBFFile,'qR','rqR');
            rCBVfnm = rCBVFile; %strrep(rCBVFile,'qR','rqR');
            Tmaxfnm = TmaxFile; %strrep(TmaxFile,'qR','rqR');
            
            nii_roi2stats('aal',MTTfnm,'','cbf_',[MTTdir filesep filepaths(i).name '_MTT_' type '.mat']);
            nii_roi2stats('jhu',MTTfnm,'','cbf_',[MTTdir filesep filepaths(i).name '_MTT_' type '.mat']);
            
            nii_roi2stats('aal',rCBFfnm,'','cbf_',[rCBFdir filesep filepaths(i).name '_rCBF_' type '.mat']);
            nii_roi2stats('jhu',rCBFfnm,'','cbf_',[rCBFdir filesep filepaths(i).name '_rCBF_' type '.mat']);
            
            nii_roi2stats('aal',rCBVfnm,'','cbf_',[rCBVdir filesep filepaths(i).name '_rCBV_' type '.mat']);
            nii_roi2stats('jhu',rCBVfnm,'','cbf_',[rCBVdir filesep filepaths(i).name '_rCBV_' type '.mat']);
            
            nii_roi2stats('aal',Tmaxfnm,'','cbf_',[Tmaxdir filesep filepaths(i).name '_Tmax_' type '.mat']);
            nii_roi2stats('jhu',Tmaxfnm,'','cbf_',[Tmaxdir filesep filepaths(i).name '_Tmax_' type '.mat']);
        catch
            fileID = fopen([outputPath '/errorLog.txt'],'a');
            fprintf(fileID,'error on particitpant %s \n',filepaths(i).name);
            fclose(fileID);
        end
    end
end
