%% Load the proxy, with values before the 90s censored to zero

if strcmp(dataFrequency,'M')
    load('../../instrument/OilSurprisesMLog')

    proxyRaw = [oilProxiesWTIM(:,ncontract)]; 
    
    smplStartProxyInd = find(strcmp(sampleDatesProxy,smplStartProxy));
    smplEndProxyInd   = find(strcmp(sampleDatesProxy,smplEndProxy));

    smplStartProxyVARInd = find(strcmp(sampleDates,smplStartProxy));
    smplEndProxyVARInd   = find(strcmp(sampleDates,smplEndProxy));

end

censorDate = '1990M01';
censorDateInd = find(strcmp(sampleDatesProxy,censorDate));
proxyTemp = proxyRaw;
proxyTemp(1:censorDateInd-1,:) = 0;
proxy = proxyTemp(smplStartProxyInd:smplEndProxyInd,:);   

[T,np] = size(proxy);
k = 1; % index of variable(s) to be instrumented