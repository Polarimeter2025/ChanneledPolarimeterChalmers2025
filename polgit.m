%Polarimeter for github, basic methodology for retrieval. Note that some
%pre-processing has been made to make the files available as csv.
clear all
clc 
close all

%% Loading data from folder described by "path1" and spectrometer interval from "Spectrometer_interval"
Wavelength_interval=readmatrix("path1\Spectrometer_interval.csv");
Wavelength_interval=Wavelength_interval(:,1)+Wavelength_interval(:,2)/1000; %Spectrometers interval [nm]
alldatapath=dir(fullfile("path1\data.csv"));
refpath=dir(fullfile("path1\refmeas.csv"));
alldatalist=natsortfiles(alldatapath);%Sorts data in order of increasing file number
reflists=natsortfiles(refpath);
alldata=fullfile({alldatalist.folder},{alldatalist.name});
refdata=fullfile({reflists.folder},{reflists.name});

%% Data sorting

%The following loops separates the files into measurement data named "measlist"
%and measurement attributes ("attributes"). The attribute files are then
%further separated into lists of attribute names ("attribname") and
%corresponding numerical values ("attribval").
datalength=length(alldata);

for k=2:datalength-1
    workitem=alldata{k};
    measlist{k}=readmatrix(workitem,'FileType','text'); %Extracts data from file
    measlist{k}=measlist{k}(1:2050);
    reflist{k}=readmatrix(refdata{k},'FileType','text');
    reflist{k}=reflist{k}(1:2050);
    attrib{k}=readlines(workitem); %Reads file as list of strings
    attributes{k}=attrib{k}(2052:end);% Attributes are located at the end of the STD-files.

    for i=1:length(attributes{k}) %Adapting attribute-extraction to separator present in line
        if contains(attributes{k}(i),'=')
            attribname{k}(i)=extractBefore(attributes{k}(i),'=');
            attribvalue{k}(i)=extractAfter(attributes{k}(i),'=');
        elseif  contains(attributes{k}(i),' ')
            attribname{k}(i)=extractBefore(attributes{k}(i),' ');
            attribvalue{k}(i)=extractAfter(attributes{k}(i),' ');
        else
            attribname{k}(i)=attributes{k}(i);
            attribvalue{k}(i)=attributes{k}(i);
        end
    end
end
%For practical purposes "attribval" is still a list of strings.
% ref=reflist{20};
%% Cleaning raw data

Offset=measlist{2};
Dark=measlist{3};
int_time_dark=str2double(attribvalue{3}(22));

for k=4:datalength-1
    Corr_I{k-3}=measlist{k}-str2double(attribvalue{k}(22))*(Dark-Offset)/int_time_dark; %Corrected intensity
end
figure(10)

for k=1:47
    % for i=0:15
    % ref{k}(1+i*128:i*128+128)=rescale(reflist{21}(1+i*128:i*128+128),min(Corr_I{k}((1+i*128:i*128+128))),max(Corr_I{k}((1+i*128:i*128+128))));
    % end
    plot(Corr_I{31})
   max(Corr_I{k})
    ref{k}=rescale(reflist{21}(1:2048),min(Corr_I{k}),max(Corr_I{k}));
  A{k}=Corr_I{k};
Corr_I{k}=smooth(rescale(Corr_I{k}((1:2048)),0,1))./smooth(rescale(ref{k}(1:2048),0,1));
end
%OBS correct for exposure time (integration time)
%% Finding an approximation of the total intensity
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
%The following loops find the peaks of the oscillating intensity and fits a
%line between each peak. Thereafter the phase at each point in the spectra
%is approximated by calculating the fraction of the period the point is
%situated at for the local period.

for k=1:datalength-4
    [minval{k},minloc{k},~,~]=findpeaks(-smooth(smooth(Corr_I{k})),'MinPeakProminence',0.001,'NPeaks',40,'MinPeakDistance',30,'MinPeakWidth',5);
    [peakval{k},peakloc{k},~,~]=findpeaks(smooth(smooth(Corr_I{k})),'MinPeakProminence',0.001,'NPeaks',40,'MinPeakDistance',30,'MinPeakWidth',5);
    
    for j=2:length(minloc{k})
        line{k}(j*(1:2))=polyfit([minloc{k}(j-1),minloc{k}(j)],-[minval{k}(j-1),minval{k}(j)],1);
        fittedline{k}(minloc{k}(j-1):minloc{k}(j))=polyval(line{k}(j*(1:2)),minloc{k}(j-1):minloc{k}(j));
    end
    for j=2:length(peakloc{k})
        line2{k}(j*(1:2))=polyfit([peakloc{k}(j-1),peakloc{k}(j)],[peakval{k}(j-1),peakval{k}(j)],1);
        fittedline2{k}(peakloc{k}(j-1):peakloc{k}(j))=polyval(line2{k}(j*(1:2)),peakloc{k}(j-1):peakloc{k}(j));
    end
end

%% Calculating the polarization at each wavelength point "i", and instrument direction "k".

for k=1:datalength-4
if length(fittedline2{k})>length(fittedline{k})
    for i=1:length(fittedline{k})
        P{k}(i)=(fittedline2{k}(i)-fittedline{k}(i))/2;
    end
else
    for i=1:length(fittedline2{k})
        P{k}(i)=(fittedline2{k}(i)-fittedline{k}(i))/2;
    end
end
end

% After this of course also a lot of graphics and reformatting and stuff. 
% A code for relative sun position was used but can also be retrieved from online resources%