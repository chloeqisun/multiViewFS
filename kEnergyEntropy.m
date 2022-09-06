function eE = kEnergyEntropy(data)
% 求信号基于emd分解算法的能量熵
% 参考《面向高铁走行部故障诊断算法的研究与实现》
% 输入：
% data：待分析信号
% 输出：
% eE：能量熵值

% 需要使用MATLAB2018a及更新版本
[len,num] = size(data);
for i = 1:num
    imf = emd(data(:,i));
    imfE = sum(imf.^2,2);
    eE(i) = kInformationEntopy(imfE,length(imfE));  %奇异谱熵。分组数等于数据长度。
end
end