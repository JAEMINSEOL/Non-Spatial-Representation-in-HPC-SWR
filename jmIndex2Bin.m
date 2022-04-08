function bin = jmIndex2Bin(index)
bin = zeros(max(index(:,2)),1);
for i=1:size(index,1)
    bin(index(i,1):index(i,2),1)=1;
end