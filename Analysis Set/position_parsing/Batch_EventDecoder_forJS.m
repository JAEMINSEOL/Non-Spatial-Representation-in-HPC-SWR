
Initial_SWRFilter_common;

% mother_root = 'H:\CA1&SUB_SCSM\ephys_analysis';
mother_root = ['F:\EPhysRawData\RawData'];

Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);


%
cd(mother_root);
% [~,inputCSV] = xlsread('D:\HPC-LFP project\Information Sheet\ClusterList.csv');

stCellRun = 1;

for cellRUN = 1 : size(Session_List,1)

    thisRID = Session_List.rat{cellRUN};
    thisSID = jmnum2str(Session_List.session(cellRUN),2);
    
      
    EventDecoderFunction_forJS(thisRID,thisSID,mother_root);

end

posFileName = 'VT1.nvt';