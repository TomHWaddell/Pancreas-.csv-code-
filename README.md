# Pancreas-.csv-code-
Code to pull ROI data from .mat files to .csv using MATLAB.

searchDir='[INSERT DIRECTORY]'
% MAT.deleteUnnecessaryMatFiles(searchDir);
matFilePrefix = '/*.mat';
exportReports = false;
searchPathMat = strcat(searchDir, matFilePrefix);

mat_files = dir(searchPathMat);
numMatFiles = length(mat_files);
output = {};

labels = {'directory', 'patient_ID', 'patient_name', 'study_date', 'ROI_T2s', 'ROI_iron', 'ROI_T1', ...
    'ROI_cT1', 'ROI_fat', 'cT1_mode', 'cT1_median', 'cT1_mean', 'cT1_Q1', 'cT1_Q3', 'cT1_IQR', 'cT1_COV', ...
    't1_comment', 't2s_comment', 'fat_comment', 'additional_comment', 'analyst_ID', 'lmsD_version', 'analysis_Date', ...
    'cT1ROI1mean', 'cT1ROI2mean', 'cT1ROI3mean', 'fatROI1mean', 'fatROI2mean', 'fatROI3mean', ...
    'T1ROI1mean', 'T1ROI2mean', 'T1ROI3mean', 'patient_age', 'patient_BMI', 'patient_weight', 'patient_height', 'patient_sex'};
disp('Beginning...');

for i = 1:numMatFiles
%     try
        
        currMatFile = mat_files(i).name;
        if regexp(currMatFile,'^\.\_') % file begins with ._
            continue;
        end
        
        if(exist(fullfile(searchDir,currMatFile), 'file') == 2)
            disp(['Reading file: ',  currMatFile]);
            matFullFile = fullfile(searchDir,currMatFile);
            load(matFullFile); 
            A = cell(1, 36);
            
            
            if isprop(lms_data, 'study') && isfield(lms_data.study, 'directory')
                dirName = lms_data.study.directory;
                [a,b,~]=fileparts(dirName);
                %eid=b(1:7);
                A{1, 1} = a;
            end


            if isprop(lms_data, 'study') && isfield(lms_data.study, 'info') && isfield(lms_data.study.info, 'patientID')
                A{1, 2} = lms_data.study.info.patientID;
            end
            if isprop(lms_data, 'study') && isfield(lms_data.study, 'info') && isfield(lms_data.study.info, 'patientName')
                A{1, 3} = lms_data.study.info.patientName;
            end
            if isprop(lms_data, 'study') && isfield(lms_data.study, 'info') && isfield(lms_data.study.info, 'studyDate')
                A{1, 4} = lms_data.study.info.studyDate;
            end
            
            if isprop(lms_data, 'roi') && isfield(lms_data.roi, 't2s') && isfield(lms_data.roi.t2s, 'mean_medians')
                t2s = lms_data.roi.t2s.mean_medians/1000;
                A{1, 5} = lms_data.roi.t2s.mean_medians;
                if lms_data.study.info.MagneticFieldStrength == 3
                    
                    aLiverModel = ECF_FeModel.LiverModel_2_Comp_3T;
                    A{1, 6} = aLiverModel.ironConc(t2s);
                else
                    aLiverModel = ECF_FeModel.LiverModel_2_Comp_15T;
                    A{1, 6} = aLiverModel.ironConc(t2s);
                end
            end
                
            
            mat=[];
            if isprop(lms_data, 'images') && isfield(lms_data.images, 't1')
                slices=length(lms_data.images.t1);
                for l = 1:slices
                    eval('y=1;lms_data.images.t1{l}.roi.data;','y=0;')
                    if y ==1
                        for k = 1:length(lms_data.images.t1{l}.roi.data)
                           mat(l,k)=[lms_data.images.t1{l}.roi.data{k}.mean];
                        end
                        [ii,~,v] = find(mat);
                        A{1, 7} = mean2(v);
                    end
                end
            end
             
            mat2=[];
            if isprop(lms_data, 'images') && isfield(lms_data.images, 'ct1')
                slices=length(lms_data.images.ct1);
                for l = 1:slices
                    eval('y=1;lms_data.images.ct1{l}.roi.data;','y=0;')
                    if y ==1
                        for k = 1:length(lms_data.images.ct1{l}.roi.data)
                            mat2(l,k)=[lms_data.images.ct1{l}.roi.data{k}.mean];
                        end
                        [iii,~,u] = find(mat2);
                        A{1, 8} = mean2(u);
                        
                        
                    end
                end
            end                  
            
            mat3=[];
            if isprop(lms_data, 'images') && isfield(lms_data.images, 'pdff')
            slices=length(lms_data.images.pdff);
                for l = 1:slices
                    eval('y=1;lms_data.images.pdff{l}.roi.mean;','y=0;')
                    if y ==1
                        for k = 1:length(lms_data.images.pdff{l}.roi.mean)
                            mat3(l,k)=[lms_data.images.pdff{l}.roi.mean];
                        end
                        [iiii,~,w] = find(mat3);
                        A{1, 9} = mean2(w);
                    end
                end
            end
                        

            
            if isprop(lms_data, 'images') && isfield(lms_data.images, 'ct1')
                
                if exist('segStats','var')
                    segStats = segStats(1);
                end

                cT1segmented = any(cellfun(@(x) haveMask(x), lms_data.images.ct1));
                if 0%T1segmented
                    internal_counter = 1;
                    for s = 1:length(lms_data.images.ct1)

                        if lms_data.images.ct1{s}.haveMask
                             segStats(internal_counter) = seg.do_ct1_Stats(lms_data.images.ct1{s}.pixels(lms_data.images.ct1{s}.mask.pixels));
                             internal_counter = internal_counter + 1;
                        end
                    end
                    if internal_counter == 1
                        disp('PANIC - Andrea was wrong...');
                    end
                    A{1, 10} = mean([segStats.mode]);
                    A{1, 11} = mean([segStats.median]);
                    A{1, 12} = mean([segStats.mean]);
                    A{1, 13} = mean([segStats.Q1]);
                    A{1, 14} = mean([segStats.Q3]);
                    A{1, 15} = mean([segStats.IQR]);
                    cov = [segStats.CoV];
                    A{1, 16} = sqrt(sum(cov.^2)/numel(cov)); 
%                     if isfield(segStats, 'perc')
%                         perc = mean([segStats.perc], 2)';
%                         A{1, 18} = perc(1);
%                         A{1, 19} = perc(2);
%                         A{1, 20} = perc(3);
%                         A{1, 21} = perc(4);
%                         A{1, 22} = perc(5);
%                         A{1, 23} = perc(6);
%                     end
                end
            end
            
             if isprop(lms_data, 'images') && isfield(lms_data.images, 'ct1')
                empty1=zeros(3,1);
                
                slices=length(lms_data.images.ct1);
                for l = 1:slices
                    eval('y=1;lms_data.images.ct1{l}.roi.data;','y=0;')
                    if y ==1
                        
                        for k = 1:length(lms_data.images.ct1{l}.roi.data)
                            if isempty(lms_data.images.ct1{l}.roi.data{k}.mean)
                                empty1(k)=[];
                            else
                                empty1(k)=lms_data.images.ct1{l}.roi.data{k}.mean;
                            end
                        end
                        
                        ct1ROI=empty1;
                        
                        A{1, 24} = ct1ROI(1);
                        A{1, 25} = ct1ROI(2);
                        A{1, 26} = ct1ROI(3);
                        
                    end
                end
             end                   
            
             if isprop(lms_data, 'roi') && isfield(lms_data.roi, 'pdff')
                
                 empty2=zeros(3,1);
                
                slices=length(lms_data.images.pdff);
                for l = 1:slices
                    eval('y=1;lms_data.images.pdff{l}.roi.data;','y=0;')
                    if y ==1
                        
                        for k = 1:length(lms_data.images.pdff{l}.roi.data)
                            if isempty(lms_data.images.pdff{l}.roi.data{k}.mean)
                                empty2(k)=[];
                            else
                                empty2(k)=lms_data.images.pdff{l}.roi.data{k}.mean;
                            end
                        end
                        fatROI=empty2;
                        
                        A{1, 27} = fatROI(1);
                        A{1, 28} = fatROI(2);
                        A{1, 29} = fatROI(3);

                    end
                end
             end
             
             if isprop(lms_data, 'images') && isfield(lms_data.images, 't1')
                empty1=zeros(3,1);
                
                slices=length(lms_data.images.t1);
                for l = 1:slices
                    eval('y=1;lms_data.images.t1{l}.roi.data;','y=0;')
                    if y ==1
                        
                        for k = 1:length(lms_data.images.t1{l}.roi.data)
                            if isempty(lms_data.images.t1{l}.roi.data{k}.mean)
                                empty1(k)=[];
                            else
                                empty1(k)=lms_data.images.t1{l}.roi.data{k}.mean;
                            end
                        end
                        
                        t1ROI=empty1;
                        
                        A{1, 30} = t1ROI(1);
                        A{1, 31} = t1ROI(2);
                        A{1, 32} = t1ROI(3);
                        
                    end
                end
             end                   
             
             if isprop(lms_data, 'quality') && isfield(lms_data.quality, 't1Comment')
                A{1, 17} = lms_data.quality.t1Comment;
             end
            
             if isprop(lms_data, 'quality') && isfield(lms_data.quality, 't2sComment')
                A{1, 18} = lms_data.quality.t2sComment;
             end
            
             if isprop(lms_data, 'quality') && isfield(lms_data.quality, 'fatComment')
                A{1, 19} = lms_data.quality.fatComment;
             end
            
             if isprop(lms_data, 'quality') && isfield(lms_data.quality, 'additionalComment')
                A{1, 20} = lms_data.quality.additionalComment;
             end
            
             if isprop(lms_data, 'quality') && isfield(lms_data.quality, 'analystID')
                A{1, 21} = lms_data.quality.analystID;
             end
            
             if isprop(lms_data, 'versionInfo') && isfield(lms_data.versionInfo, 'lmsappdevTag')
                A{1, 22} = lms_data.versionInfo.lmsappdevTag;
             end
            
             if isprop(lms_data, 'quality') && isfield(lms_data.quality, 'analysisDate')
                A{1, 23} = lms_data.quality.analysisDate;
             end
             
             key1 = lms_data.images.t2s{1}.scanKeys{1};
             tags = lms_data.study.fileMap.scanMap(key1).dicom{1}.tags;
             
             A{1, 33} = tags.PatientAge;
             weight = tags.PatientWeight;
             height = tags.PatientSize;
             BMI = weight/((height)^2);
             A{1, 34} = BMI;
             A{1, 35} = tags.PatientWeight;
             A{1, 36} = tags.PatientSize;
             A{1, 37} = tags.PatientSex;

            output(:, end+1) = A;
        end
            
%     catch ME
%         err=lasterror;
%         disp(err.message);
%         error=1;
%         if error==1
%             filename = 'errorlog.txt';
%             fid4 = fopen(filename, 'a');
%             fprintf(fid4, '%s,%s,%s \n', datetime('now'),b,err.message);
%             fclose(fid4)
%         end
% 
%         continue
%     end
end

%AA = cell2table(output');
AA.Properties.VariableNames = labels;
timeStamp=char(datetime('now'));
outputFile=strcat('AnalysisResults_v4_',timeStamp(1:11),'.csv');
writetable(AA, fullfile(searchDir, outputFile));
disp('Finished');
