clear all;
fftsize = 4096;
fstep=1;
refractive_index_in_tissue=1.40;
resized_z = 512;

% folder of OCT scans of ~8 regions of a patient's brain
cd("D:\Brain_Cancer_classification\Exvivo_Data\Grey Matter_or_Cancer\2018-04-23\ExVivo\S1_Tumor")

data_folder = "C:\Users\Nathan\Deep Learning\brain_cancer_oct\invivo_cancer_classifier\training_data\";
data_label = "high_grade_cancer\"; % "high_grade_cancer", "low_grade_cancer", or "non_cancer"
parent_folder = append(data_folder, data_label);
patient_id = 'sample15_s1';
ofolder2=append(parent_folder, patient_id); % where the processed B-frames are stored

fheader=dir('*.header.tiff');
fheader=fheader(1).name;
finfo=imfinfo(fheader);

% OCT image properties from header.tiff (assume all header files contain same metadata)
DepthSize=finfo.Width;
LateralSize=finfo.Height;
BytePerSample=finfo.BitsPerSample/8;
FrameSize=DepthSize*LateralSize*BytePerSample;
disp(['Each frame of raw data has ' num2str(DepthSize) ' (D) x ' num2str(LateralSize) ' (L) pixels.']);

raw=dir('*.oct.raw');

% loop over each oct.raw file
for k=1:length(raw)
    fraw=raw(k).name;
    fileInfo = dir(fraw); % struct with fields: name, folder, date, bytes, isdir, datenum
    fileSize = fileInfo.bytes;
    FrameNumber=fileSize/FrameSize;
    % disp(['The raw file ' fraw ' contains ' num2str(FrameNumber) ' frames.']);

    filename_frag=regexp(fraw,'[_.]','split');
    fileprefix='';

    for n=1:length(filename_frag)-3
        fileprefix=[fileprefix cell2mat(filename_frag(n)) '_'];
    end

    startFrameNum=str2double(cell2mat(filename_frag(end-2)));
    fftw('planner','patient');
    str=fftw('wisdom');
    w = repmat(single(blackman(DepthSize,'periodic')),1,LateralSize);
    fid=fopen(fraw,'rb'); 
    t=startFrameNum+1:fstep:startFrameNum+FrameNumber;
    Phantomdata = zeros(2048,2048);
    
    % iteraete over each frame
    for n=1:length(t)
        if exist('display_str','var')==0
            display_str=sprintf('%05d',n);
        else
            display_str=sprintf('\b');
            display_str=repmat(display_str,1,6);
            display_str=sprintf('%s%05d',display_str,n);
        end

        disp(display_str);
        fseek(fid,FrameSize*(fstep-1),'cof');

        if BytePerSample==1
            A=fread(fid,DepthSize*LateralSize,'uint8=>single');
        else
            A=fread(fid,DepthSize*LateralSize,'uint16=>single');
        end

        if fstep~=1
            fseek(fid,FrameSize*(fstep-1),'cof');
        end

        A=A-single(2^15);
        RawData=reshape(A,DepthSize,LateralSize);
        
        % process B-frame
        I_signal2=fft(RawData.*w,fftsize,1);
        image=abs(I_signal2(1:fftsize/2,:));
        %imm = imresize(log(image),[256 256]);
        imm = log(image);
   
        [si,ei]=regexpi(fraw,'^.+[0-9]{5}');
        ftemp=fraw(si:ei-5);
        ftemp2=fraw(si:ei-16);
        
        maxP_u16=max(max(imm));
        minP_u16=min(min(imm));

        imm=min(imm,maxP_u16);
        imm=max(imm,minP_u16);
        imm=(imm-minP_u16)/(maxP_u16-minP_u16);

        fname2=sprintf('%s\\%s%04d.oct.csv',ofolder2,ftemp,n+startFrameNum);
        writematrix(1-imm, fname2);     
    end

    fclose(fid);
end
