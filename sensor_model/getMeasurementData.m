function measurement_data = getMeasurementData()
load('tactip-benchmarks\trainDepthVideo05012246\c01_01');
load('tactip-benchmarks\trainDepthVideo05012246/expt_mapping.mat')
for i = 1:length(data)
    tmp = data{i};
    pins_xy = [];
    for j = 1:size(tmp, 2)
        x = mean(tmp(:,j,1));
        y = mean(tmp(:,j,2));
        pins_xy = [pins_xy, [x;y]];
    end
    pins_xy = pins_xy(:,expt_mapping);
    
    measurement_data{i} = pins_xy;
end
    
end