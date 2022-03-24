
function EventDecoderFunction_forJS(thisRID,thisSID,mother_root)
session_root = [mother_root '\rat' thisRID '\rat' thisRID '-' thisSID];  
cd(session_root);
evtfile_listing = dir('Events.csv');

nTable_temp = readtable([session_root '\Events.csv']);
uTable = readtable([session_root '\Rat' thisRID '_' thisSID '.csv']);
if ~exist('ExtractedEvents','dir')
    mkdir('ExtractedEvents');
end
cd('ExtractedEvents');
%     listing=dir('trial*.mat');[r_trial,c_trial]=size(listing);
% %     if r_trial, return; end

BehStart = find(strcmp("Starting Recording",nTable_temp.Var18));
BehStop= find(strcmp("Stopping Recording",nTable_temp.Var18));

nTable = table;
for n=1:size(BehStart,1)
nTable = [nTable; nTable_temp(BehStart(n)+1:BehStop(n)-1,:)];
end

nTable_decod=table;
for n=1:size(nTable,1)
    nTable_decod.time(n) = nTable.Var4(n);
    nTable_decod.event(n) = hex2dec(nTable.Var18{n}(end-3:end-2));
end

clear nTable nTable_temp

nTable = nTable_decod;

nTable.pos(nTable.event==16) = uTable.Var4(strcmp(uTable.Var1,'X'));
nTable.pos(nTable.event~=16) = interp1(nTable.time(nTable.event==16),nTable.pos(nTable.event==16),nTable.time(nTable.event~=16));
uTable.time(strcmp(uTable.Var1,'X')) = nTable.time(nTable.event==16);

nTable.pos = nTable.pos-fix(nTable.pos/10000)*10000;


figure;
x  = nTable.pos(nTable.event==16);
scatter(x,ones(length(x),1),10,[.8 .8 .8],'filled')
hold on
x  = nTable.pos(nTable.event==4);
scatter(x,ones(length(x),1),10,'k','filled')
x  = nTable.pos(nTable.event==12);
scatter(x,ones(length(x),1),10,'b','filled')

x  = nTable.pos(nTable.event==8);
scatter(x,ones(length(x),1),10,'r','filled')

find(nTable.event ==4)
fclose(xlsID);






