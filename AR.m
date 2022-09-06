function [A,Sgm2] = AR(x,P)

% Input:
% x: The segment of signal
% P: order of AR
% Output:
% A: AR coefficient [a1 a2 ...aP]
% Sgm2: variance
A  = zeros(1,P); 
Rx = Rxx(x,P);

A(1) = - Rx(1+1)/Rx(0+1);
Sgm2  = Rx(0+1)*(1 - A(1)^2);

for p = 1:P-1
    k = 1:p;
    K = -(Rx((p+1)+1) + A(k)*Rx((p+1-k)+1))/Sgm2;
    Sgm2 = Sgm2 * (1 - K*K);
    
    A(k) = A(k) + K * A(p+1-k);
    A(p+1) = K;
end