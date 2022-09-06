function attrVector=attributionVector(S,N,featureNum,classValue)

%%%%%%%%%
% Input:
% S: The segment of signal
% N: The number of sampling points
% featureNum: The number of features
% classValue: Class label, classValue is 0 or 1


% For each selected channel
for channel=1:5        
    
    %% Time domain
    
    % Maximum
    Max=max(S(:,channel));   
	
    % Minimum
    Min=min(S(:,channel));
	
    % Mean
    mu=mean(S(:,channel));
	
    % Variance
    sigma=std(S(:,channel));
	
	% Line length
    LL=sum(abs(S(1:N-1,channel)-S(2:N,channel)));
	
	% Skewness
    sk=skewness(S(:,channel));
	
	% Kurtosis
    ku=kurtosis(S(:,channel));
	
	% Nonlinear Energy
    nonlinearEnergy=S(2:N-1,channel)'*S(2:N-1,channel)-S(1:N-1,channel)'*S(2:N,channel);
    
	% Root Mean Square Amplitude
    rm=rms(S(:,channel));
	
	% The variance of the first difference signal
	firstOrder=S(2:N,channel)-S(1:N-1,channel); 
    sigmaFirst=std(firstOrder); 
	
	% The variance of the second difference signal
	secondOrder=firstOrder(2:length(firstOrder))-firstOrder(1:length(firstOrder)-1);
    sigmaSecond=std(secondOrder);
	
	% Activity
    activity=sigma^2;
	
	% Mobility
    mobility=sigmaFirst^2;
	
	% Complexity
    complexity=sigmaSecond^2;
	
	% Zero Crossing
    for i=1:N
        d(i)=sign(-S(i,channel)*S(i,channel))+1;
        if i~=N
           dFirst(i)=sign(-firstOrder(i)*firstOrder(i))+1;
        end
        if i~=N && i~=N-1
            dSecond(i)=sign(-secondOrder(i)*secondOrder(i))+1;
        end
    end
	% The number of zero crossings of original signal
    zeroCrossNum=sum(d/2);
	% The number of zero crossings of first derivative
    zeroCrossFirstNum=sum(dFirst/2);
	% The number of zero crossings of second derivative
    zeroCrossSecondNum=sum(dSecond/2);
	
	% Number of Maxima and Minima
    maxima=findpeaks(S(:,channel));
    minima=findpeaks(-S(:,channel));
    maxminNum=length(maxima)+length(minima);
	
    % Wave Form Factor
    av=mean(abs(S(:,channel)));		
    formF=rm/av;
	
	
    % Peak Factor 
    pk = Max-Min;			        
    peakF=pk/rm;
	
    % Pulse Factor
    pulseF=pk/av;
	
    % Margin Factor
    marginF=pk/mean(sqrt(abs(S(:,channel)))).^2;
    
    % Errors of AR modeling with order 1–9
    [A1,AME1]=AR(S(:,channel),1);
    [A2,AME2]=AR(S(:,channel),2);
    [A3,AME3]=AR(S(:,channel),3);
    [A4,AME4]=AR(S(:,channel),4);
    [A5,AME5]=AR(S(:,channel),5);
    [A6,AME6]=AR(S(:,channel),6);
    [A7,AME7]=AR(S(:,channel),7);
    [A8,AME8]=AR(S(:,channel),8);
    [A9,AME9]=AR(S(:,channel),9);
    timeFeatures=[Max Min mu formF peakF pulseF marginF sigma sigmaFirst sigmaSecond LL rm ...
        activity mobility complexity zeroCrossNum zeroCrossFirstNum zeroCrossSecondNum ...
        nonlinearEnergy sk ku maxminNum AME1 AME2 AME3 AME4 AME5 AME6 AME7 AME8 AME9];
    
    %% Frequency domain
	
    % Total Power
    power=abs(fft(S(:,channel))).^2/N;
    [p,f] = periodogram(S(:,channel),[],[],256);
    powerT=sum(abs(S(:,channel)).^2);
	
    % Peak Frequency
    peakPower=max(p);
%     peakPower=max(fft(S(:,channel)));

    % Median Frequency
    medianPower=mean(p);  
	
    % Center Frequency
    FC = sum(p.*f)./sum(p);
	
    % Frequency Variance
    VF = sum((f-FC).^2.*p)./sum(p);
	
	% Root Mean Square Frequency
    MSF = sum(f.^2.*p)./sum(p); 
    
    % Wavelet Entropy
    wpt=wpdec(S(:,channel),5,'db4');      
    for j=1:2^5
        E(j)=sum(abs(wprcoef(wpt,[5,j-1])).^2);
    end
    E1=sum(E);
    dim=length(E);
    for j=1:dim
        p(j)=E(j)/E1;
    end
    waveletEntropy=-sum(p.*log(p)); 
	
    % Wavelet Energy
    waveletEnergy=E(5);
	
	% Spectral Edge Frequencies（80%，90%，95%）
    SEF1=getSEF(S(:,channel),0.8);
    SEF2=getSEF(S(:,channel),0.9);
    SEF3=getSEF(S(:,channel),0.95);
	
    frequencyFeatures=[powerT peakPower medianPower MSF FC VF SEF1 SEF2 SEF3 waveletEnergy];
    
    %% Information theory
	
	% Spectrum Entropy
    psdE = kPowerSpectrumEntropy(S(:,channel));
	
	% Singular Value Entropy
    svdpE = kSingularSpectrumEntropy(S(:,channel),N/2);  
	
    % Sample Entropy
    r=0.1*std(S(:,channel));
    sampleEntropy=SampEn(S(:,channel)', 2, r);  
    
    % Energy Entropy
    energyEntropy = kEnergyEntropy(S(:,channel));
    
    infoEntropyFeature=[waveletEntropy sampleEntropy psdE svdpE energyEntropy];
    
    attrVector(1,(channel-1)*featureNum+1:featureNum*channel)=[timeFeatures frequencyFeatures infoEntropyFeature];
end
attrVector(1,channel*featureNum+1)=classValue;
end