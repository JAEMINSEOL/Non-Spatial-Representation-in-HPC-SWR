Initial_SWRFilter_common;

thisRID = 551;
raw_dir = ['Y:\EPhysRawData\CA1,CA3,V2L Recording_VR_JS\Rat' jmnum2str(thisRID,3) ];
save_dir = 'F:\EPhysRawData\RawData';
for thisSID = 1:4
RateReducer(thisRID, thisSID, raw_dir, save_dir,'AG')
end
