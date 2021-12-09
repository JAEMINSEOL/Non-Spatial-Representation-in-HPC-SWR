Initial_SWRFilter_common;
thisRegion = 'CA1';
%% export units & ripples table
for clRUN = 1:5
    thisSID.full = ['561-0' num2str(clRUN)];
    load([ROOT.Save '\' thisSID.full '_' thisRegion '.mat'])
    if any(strcmp('FieldPos',fieldnames(UnitsTable)))
        UnitsTable = removevars(UnitsTable, 'FieldPos');
    end
    for i=1:size(UnitsTable,1)
        if UnitsTable.Type(i)>0
            UnitsTable.Field(i)=1;
        else
            UnitsTable.Field(i)=0;
        end
    end
    UnitsTable.Properties.VariableNames{1} = 'Region';
    UnitsTable.Region(:)=1;
    
    
    UnitsTable.PeakArea(UnitsTable.PeakPos<=2) = 1;
    UnitsTable.PeakArea(UnitsTable.PeakPos<29 & UnitsTable.PeakPos>2) = 2;
    UnitsTable.PeakArea(UnitsTable.PeakPos<38 & UnitsTable.PeakPos>28) = 3;
    UnitsTable.PeakArea(UnitsTable.PeakPos>37) = 4;
    
     
    
    UnitsTable = movevars(UnitsTable, 'Region', 'After', 'TT');
    UnitsTable = movevars(UnitsTable, 'PeakArea', 'After', 'Field');
    
    
    for i=1:size(RipplesTable,1)
        id = cell2mat(RipplesTable.RippleID(i));
    RipplesTable.Rat(i) = str2double(id(1:3));
    RipplesTable.Session(i) = str2double(id(5:6));
    RipplesTable.RipID(i) = str2double(id(end-3:end));
    
    if strcmp(RipplesTable.Region(i),'SUB')
        RipplesTable.temp1(i) = 0;
    elseif strcmp(RipplesTable.Region(i),'CA1')
        RipplesTable.temp1(i) = 1;
    elseif strcmp(RipplesTable.Region(i),'CA3')
        RipplesTable.temp1(i) = 3;
    elseif strcmp(RipplesTable.Region(i),'DG')
        RipplesTable.temp1(i) = 4;
    end
    
    
    if strcmp(RipplesTable.Area(i),'stbox')
    RipplesTable.temp(i) = 1;
    elseif strcmp(RipplesTable.Area(i),'stem')
        RipplesTable.temp(i) = 2;
    elseif strcmp(RipplesTable.Area(i),'arm')
        RipplesTable.temp(i) = 4;
    end
end


RipplesTable = removevars(RipplesTable, 'Area');
RipplesTable.Area = RipplesTable.temp;
RipplesTable = removevars(RipplesTable, 'temp');

RipplesTable = removevars(RipplesTable, 'Region');
RipplesTable.Region = RipplesTable.temp1;
RipplesTable = removevars(RipplesTable, 'temp1');

RipplesTable.Epoch = mod(RipplesTable.Trial*10,10);
RipplesTable.Trial = floor(RipplesTable.Trial);

RipplesTable = movevars(RipplesTable, 'Rat', 'Before', 'RippleStartIndex');
RipplesTable = movevars(RipplesTable, 'Session', 'After', 'Rat');
RipplesTable = movevars(RipplesTable, 'Region', 'After', 'Session');
RipplesTable = movevars(RipplesTable, 'RipID', 'After', 'Region');
RipplesTable = movevars(RipplesTable, 'Context', 'After', 'RipID');
RipplesTable = movevars(RipplesTable, 'Trial', 'Before', 'Context');
RipplesTable = movevars(RipplesTable, 'Area', 'After', 'Context');
RipplesTable = movevars(RipplesTable, 'Epoch', 'After', 'Trial');
RipplesTable = movevars(RipplesTable, 'Filter', 'After', 'Area');
RipplesTable = removevars(RipplesTable, 'RippleID');

RipplesTable.Area(RipplesTable.Epoch>=5 | RipplesTable.Epoch==1) = 1;
RipplesTable.Area(RipplesTable.Epoch==0 | RipplesTable.Area==0) = 2;


writetable(UnitsTable,['UnitsTable_r561_s', jmnum2str(clRUN,2), '_CA1.xlsx'])
writetable(RipplesTable,['RipplesTable_r561_s' jmnum2str(clRUN,2) '_CA1.xlsx'])
end

%% make reactivation table

for clRUN = 1:5
    ReactTable_all=[];
    thisSID.full = ['561-0' num2str(clRUN)];
    load([ROOT.Save '\' thisSID.full '_' thisRegion '.mat'])
    for clRip = 1:size(RipplesTable,1)
        ReactTable_Rip=[];
        id = cell2mat(RipplesTable.RippleID(clRip));
        thisUnits = RipplesStruct.CA1.(['ripple' id(end-3:end)]).units;
        for clUnit=1:size(thisUnits,1)
            ReactTable_temp = [];
            if any(strcmp('SpksInRipple',fieldnames(thisUnits)))
            Spks = cell2mat(thisUnits.SpksInRipple(clUnit));
            
            if any(strcmp('celltype',fieldnames(thisUnits)))
                type = thisUnits.celltype(clUnit);
            elseif any(strcmp('Type',fieldnames(thisUnits)))
                type = thisUnits.Type(clUnit);
            else
                type=0;
            end
            
            if size(Spks,1)>=2, burst=1; else, burst=0; end
            
            for i = 1:size(Spks,1)
                temp = [Spks(i), thisUnits.ratID(clUnit), thisUnits.sessionID(clUnit), thisUnits.TT(clUnit), thisUnits.clusterNo(clUnit),1,type, burst];
                ReactTable_temp = [ReactTable_temp; temp];
            end
            end
            ReactTable_Rip = [ReactTable_Rip; ReactTable_temp];
        end
        if ~isempty(ReactTable_Rip)
        ReactTable_Rip = sortrows(ReactTable_Rip,1);
        ReactTable_Rip = [ones(size(ReactTable_Rip,1),1)*str2double(id(end-3:end)) ReactTable_Rip];
        end
        ReactTable_all = [ReactTable_all; ReactTable_Rip];
    end
    
    ReactTable = array2table(ReactTable_all,...
        'VariableNames',{'RipID','SpkTime','RatID','SessionID','TT','UnitID','Region','Type','Bursting'});
    writetable(ReactTable,['ReactTable_r561_s' jmnum2str(clRUN,2) '_CA1.xlsx'])
end
