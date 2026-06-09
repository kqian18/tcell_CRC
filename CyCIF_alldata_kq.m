dataCell = cell(length(slideName), 1); 

for i = 1:length(slideName)
    data1 = dataStruct1.(slideName{i});
    data1.slideName = repmat(slideName(i), size(data1,1), 1);
    dataCell{i} = data1;
end

alldata = vertcat(dataCell{:});

clear data1 dataCell i;