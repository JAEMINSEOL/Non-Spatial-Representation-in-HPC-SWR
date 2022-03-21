
Initial_SWRFilter_common;

% mother_root = 'H:\CA1&SUB_SCSM\ephys_analysis';
mother_root = ['F:\EPhysRawData\RawData'];

Cluster_List = readtable([ROOT.Info '\ClusterList_SWR.xlsx']);


%
cd(mother_root);
% [~,inputCSV] = xlsread('D:\HPC-LFP project\Information Sheet\ClusterList.csv');

stCellRun = 1;

for cellRUN = 1 : size(Cluster_List,1)
    if Cluster_List.rat(cellRUN)==27 && Cluster_List.session(cellRUN)==1
    
    thisCLUSTER = Cluster_List.ID{cellRUN};
           
    createParsedSpike(mother_root, thisCLUSTER);
    disp([thisCLUSTER ' has processed.']);
    end
end