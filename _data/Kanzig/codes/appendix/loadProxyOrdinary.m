%% Load the proxy, only keeping ordinary announcements

if strcmp(dataFrequency,'M')
    load('../../instrument/OilSurprisesMLog')

    proxyRaw = [oilProxiesWTIM(:,ncontract)]; 
    proxyRaw(statementInfoM.type(:,1)~=1) = 0; % only keep ordinary meetings proxyRaw(statementInfoM.type(:,1)~=1) = 0;

    smplStartProxyInd = find(strcmp(sampleDatesProxy,smplStartProxy));
    smplEndProxyInd   = find(strcmp(sampleDatesProxy,smplEndProxy));

    smplStartProxyVARInd = find(strcmp(sampleDates,smplStartProxy));
    smplEndProxyVARInd   = find(strcmp(sampleDates,smplEndProxy));
end

proxy = proxyRaw(smplStartProxyInd:smplEndProxyInd,:);   % we loose the first p values in the VAR

[T,np] = size(proxy);
k = 1; % index of variable(s) to be instrumented