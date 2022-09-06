function svdpE = kSingularSpectrumEntropy(data,n)
% Input：
% data：The segment of signal
% n：window length，notes：2<=n<=lenght(data)-1
% Output：
% svdpE：Singular Spectrum Entropy

[len,num] = size(data);
m = len - n - 1; 
svdpE = []; 
for j = 1:num
    A = []; 
    for i = 1:m
        A = [A;data(i:i+length(data)-m,j)'];
    end

    svdVal = svd(A); 

    svdpE(j) = kInformationEntopy(svdVal,length(svdVal));  
end
end
