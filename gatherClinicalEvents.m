function [outS,drug,fullName,aliasNames]= gatherClinicalEvents(outS,clinDir)
%assigns clinical evennts to ExPANdS output. Two fields are added to each
%input sample:
%%1 - events --> struct with fields:
%%name - name of the event (Surgery or therapy start)
%%start
%%end
%%unit - absolute (that is a date-format) or relative (in days since
%%diagnosis)
%%2 - date --> (one of the above) the event that did result in this sample to be obtained


disease=tdfread([clinDir,filesep,'MergedClinical.txt']);
portions=readtable([clinDir,filesep,'Portion.txt']);
drug=tdfread([clinDir,filesep,'Drugs.txt']);
portions.bcr_portion_barcode = char(portions.bcr_portion_barcode);
portions.bcr_sample_barcode= char(portions.bcr_sample_barcode);
%%unify drug names
newDrugNames=cellstr(lower(drug.drug_name));
newDrugNames=strrep(newDrugNames,'-','');
newDrugNames=strrep(newDrugNames,' ','');
newDrugNames=strrep(newDrugNames,'06gb','06bg');
newDrugNames=strrep(newDrugNames,'5fu','5flu');
newDrugNames=strrep(newDrugNames,'9ac9aminocamplotecian','9aminocamptothecin');
newDrugNames=strrep(newDrugNames,'9immunoaminocamptnetecin','9aminocamptothecin');
newDrugNames=strrep(newDrugNames,'almita','alimta');
newDrugNames=strrep(newDrugNames,'hyroxyurea','hydroxyurea');
newDrugNames=strrep(newDrugNames,'irunotecan','irinotecan');
newDrugNames=strrep(newDrugNames,'ironotecan','irinotecan');
newDrugNames=strrep(newDrugNames,'lumustine','lomustine');
newDrugNames=strrep(newDrugNames,'palixtaxel','paclitaxel');
newDrugNames=strrep(newDrugNames,'sovatenib','sorafenib');
newDrugNames=strrep(newDrugNames,'tanceva','erlotinib');
newDrugNames=strrep(newDrugNames,'tarceva','erlotinib');

xx=char(newDrugNames);
[a,ia]=unique(cellstr(xx(:,1:3)));
fullName=cellstr(xx(ia,:));
aliasNames=cell(1,length(a));
for i=1:length(a)
    idx=find(strcmp(cellstr(xx(:,1:3)),char(a(i))));
    uIDs=unique(cellstr(drug.drug_name(idx,:)));
    aliasNames{i}=cellstr([repmat(',',length(uIDs),1),char(uIDs)]);
    newDrugNames(idx)=fullName(i);
end
drug.drug_name=char(newDrugNames);

radiation=tdfread([clinDir,filesep,'Radiation.txt']);
% else
%     disease=tdfread('D:\SVNImports\Mutations5\Clinical.txt');
%     portions=tdfread('D:\SVNImports\Mutations5\Portion.txt');
%     drug=tdfread('D:\SVNImports\Mutations5\Drug.txt');
%     radiation=tdfread('D:\SVNImports\Mutations5\Radiation.txt');
% end

BEND = strfind(portions.bcr_portion_barcode(1,:), '-')-1;

for i=1:length(outS)
    patient=outS(i).patient;
    
    events=struct('name',[],'start',[],'end',[],'unit',[],'response',[],'drugtype',[]);
    events(1)=[];
    %%portion
    idx=find(strcmp(cellstr(char(portions.bcr_portion_barcode(:,1:BEND(3)))),patient));
    if length(idx)>1
        [~,ia]=unique(cellstr(portions.bcr_portion_barcode(idx,1:BEND(4))));
        idx=idx(ia);
    end
    for ii =1:length(idx)
        sample_type=str2double(portions.bcr_portion_barcode(idx(ii),(BEND(3)+2):(BEND(3)+3)));
        if sample_type>2
            continue;
        end
        n=length(events)+1;
        events(n).name=['Surgery_' num2str(sample_type)];
        events(n).start=datenum(portions.date_of_creation(idx(ii),:));
        events(n).end=events(n).start;
        events(n).unit='absolute';
    end
    
    %%drugs
    idx=find(strcmp(cellstr(char(drug.bcr_patient_barcode)),patient));
    for ii =1:length(idx)
        n=length(events)+1;
        events(n).start=drug.days_to_drug_therapy_start(idx(ii),:);
        events(n).end=drug.days_to_drug_therapy_end(idx(ii),:);
        events(n).unit='relative';
        events(n).name=strtrim(drug.drug_name(idx(ii),:));
        events(n).drugtype=strtrim(drug.therapy_type(idx(ii),:));
        events(n).response=drug.measure_of_response(idx(ii),:);
    end
    
    %%radiation
    idx=find(strcmp(cellstr(char(radiation.bcr_patient_barcode)),patient));
    for ii =1:length(idx)
        n=length(events)+1;
        events(n).start=str2double(radiation.days_to_radiation_therapy_start(idx(ii),:));
        events(n).end=str2double(radiation.days_to_radiation_therapy_end(idx(ii),:));
        events(n).unit='relative';
        events(n).name='Radiation';
        events(n).drugtype='Radiation';
        events(n).response=radiation.measure_of_response(idx(ii),:);
    end
    
    %%dtd, progression
    n=length(events)+1;
    idx=find(strcmp(cellstr(char(disease.bcr_patient_barcode)),patient));
    events(n).start=disease.days_to_death(idx,:);
    if ~ isnumeric(events(n).start)
        events(n).start=str2num(events(n).start);
    end
    events(n).end=events(n).start;
    events(n).unit='relative';
    events(n).name='Death';
    
    %%days to recurrence
    n=length(events)+1;
    idx=find(strcmp(cellstr(char(disease.bcr_patient_barcode)),patient));
    events(n).start=disease.days_to_tumor_recurrence(idx,:);
    if ~ isnumeric(events(n).start)
        events(n).start=str2num(events(n).start);
    end
    events(n).end=events(n).start;
    events(n).unit='relative';
    events(n).name='Recurrence';
    
    
    %     n=length(events)+1;
    %     events(n).name='lastfollowup';
    %     events(n).start=disease.days_to_last_followup(idx,:);
    %     if ~ isnumeric(events(n).start)
    %         events(n).start=str2num(events(n).start);
    %     end
    %     events(n).end=events(n).start;
    %     events(n).unit='relative';
    
    
    %days_to_additional_surgery_locoregional_procedure
    n=length(events)+1;
    events(n).name='new_tumor_event';
    xx=disease.days_to_new_tumor_event_after_initial_treatment(idx,:);
    if ~ isnumeric(xx)
        xx=str2num(xx);
    end
    xx=xx+events(1).start;
    events(n).start=xx;
    events(n).end=events(n).start;
    events(n).unit='absolute';
    
    
    %% Prettify events
    [~,ia]=sort([events.start],'descend');
    for k = ia
        if ~strcmp(events(k).unit,'relative')
            events(k).end= events(k).end- events(1).start;
            events(k).start = events(k).start - events(1).start;
            events(k).unit = 'absolute';
        end
    end
    [~,ia]=sort([events.start]);
    events = events(ia);

    outS(i).events=events;
    idx=find(strcmp(cellstr(char(events.name)),['Surgery_' num2str(outS(i).sample_type)]));
    outS(i).date=events(idx).start;%%%set date of sugery from which this sample has been obtained
end


end
