function [data, coefs, borig] = gererateDataColNnD(Dictionary,N, L, minx, SNRdB)
%[Aorig, btrain, xorig, borig] = gererateSyntheticDictionaryAndData(numtrain, div, rows, cols,  params.SNRdB);

randn('state',sum(100*clock));
rand('state',sum(100*clock));


[data,coefs] = CreateDataFromDictionarySimple(Dictionary, N, L, minx);
borig= data;

if (SNRdB==0) | (SNRdB == 80) 
    return
else
    noise = randn(size(data));
    actualNoise = calcNoiseFromSNR(SNRdB,data, noise);
    SNR = calcSNR(data, data+actualNoise);
    data =  data + actualNoise*SNR/SNRdB;   
end

function [D,xOrig] = CreateDataFromDictionarySimple(dictionary, numElements, numCoef, minx)
maxRangeOfCoef = 1;
resolution = 0.0001;

xOrig = zeros(size(dictionary,2),numElements);
%vecOfValues = -1*maxRangeOfCoef:resolution:maxRangeOfCoef;
%coefs = randsrc(numCoef,numElements,vecOfValues);
%randn
coefs = randn(numCoef,numElements)*maxRangeOfCoef;
%rand 
%coefs = rand(numCoef,numElements)*2-1; 
coefs = coefs+minx*sign(coefs);
xOrig(1:numCoef,:) = coefs;
for i=1:size(xOrig,2)
    xOrig(:,i) = xOrig(randperm(size(xOrig,1)),i);
end
for i = 1:numElements
    nrm = norm(xOrig(:, i) );
    xOrig(:, i) = xOrig(:, i) / (nrm);
end
%dictionaryElementIndices = randsrc(numCoef*numElements,1,[1:size(dictionary,2)])   ; 
%matrixOfIndices = repmat([1:numElements],numCoef,1);
%xOrig(sub2ind(size(xOrig),dictionaryElementIndices,matrixOfIndices(:))) = coefs;
D = dictionary*xOrig;

function  actualNoise = calcNoiseFromSNR(TargerSNR, signal, randomNoise)
signal = signal(:);
randomNoiseRow = randomNoise(:);
signal_2 = sum(signal.^2);
ActualNoise_2 = signal_2/(10^(TargerSNR/10));
noise_2 = sum(randomNoiseRow.^2);
ratio = ActualNoise_2./noise_2;
actualNoise = randomNoiseRow.*repmat(sqrt(ratio),size(randomNoiseRow,1),1);
actualNoise = reshape(actualNoise,size(randomNoise));

function SNR = calcSNR(origSignal, noisySignal)
errorSignal = origSignal-noisySignal;
signal_2 = sum(origSignal.^2);
noise_2 = sum(errorSignal.^2);

SNRValues = 10*log10(signal_2./noise_2);
SNR = mean(SNRValues);
