
function SampEnVal = SampEn(data, m, r)

% Input£º
% data: The segment of signal
% m: refactoring dimension
% r: threshold£¬r=0.1~0.25*Std(data)
% Output£º
% SampEnVal: Sample Entropy

data = data(:)';
N = length(data);
Nkx1 = 0;
Nkx2 = 0;

% Calculate the distance, where x1 is the sequence of length m and x2 is the sequence of length m+1

for k = N - m:-1:1
    x1(k, :) = data(k:k + m - 1);
    x2(k, :) = data(k:k + m);
end

for k = N - m:-1:1
    % calculate sequence x1
    x1temprow = x1(k, :);
    x1temp    = ones(N - m, 1)*x1temprow;
    
    % Calculate the distance. The maximum number of subtracted elements in each row is the distance
    dx1(k, :) = max(abs(x1temp - x1), [], 2)';
    
    % Number of template matches
    Nkx1 = Nkx1 + (sum(dx1(k, :) < r) - 1)/(N - m - 1);  
    
    % calculate sequence x2
    x2temprow = x2(k, :);
    x2temp    = ones(N - m, 1)*x2temprow;
    dx2(k, :) = max(abs(x2temp - x2), [], 2)';
    Nkx2      = Nkx2 + (sum(dx2(k, :) < r) - 1)/(N - m - 1);
end

% Average
Bmx1 = Nkx1/(N - m);
Bmx2 = Nkx2/(N - m);

% Sample Entropy
SampEnVal = -log(Bmx2/Bmx1);

end
