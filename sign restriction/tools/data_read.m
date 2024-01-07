
% Data
filename                    = 'quarterly_24june_2022'; %don-t have {'#N/A'} in the excel, just blank
sheet                       = 'quarterly';
range                       = 'B2:JN254';
range2                      = 'B1:JN1'; % it can be done in one step but it was not working for Marco
[RawData,~]                 = xlsread(filename,sheet,range);
[~,var_names]               = xlsread(filename,sheet,range2); % it can be done in one step but it was not working for Marco
% Assess name to each variable
for i = 1:size(RawData,2)
    eval([var_names{i} ' = RawData(:,i);']);
end

% Variables transformation
logGDP       = log(GDPC1); % log real gdp
GDPG         = [NaN; diff(logGDP)]; % real GDP growth rate
RGDPG        = logGDP - lagmatrix(logGDP,4); % real gdp growth rate
logC         = log(PCECC96); % log total real consumption
logCD        = log(PCDGx); % log real consumption of durables
logCS        = log(PCESVx); % log real consuption of services
logCND       = log(PCNDx); % log of real consumption of non-durables
logCNDS      = log(PCNDx+PCESVx); % log of real consumption of (non-durables plus services)
logI         = log(GPDIC1); % log of real gross private domestic investment
logICD       = log(GPDIC1+PCDGx); % log of real gross private domestic investment plus durable consumption
logPAYEMS    = log(PAYEMS); % All employees: total nonfarm (thousand of persons)
logH         = log(HOANBS); % Nonfarm Business Sector: Hours of All Persons (Index 2012=100)
logPCE       = log(PCECTPI); % log of Personal Consumption Expenditures: Chain-type Price Index (Index 2012=100)
logCPI       = log(CPIAUCSL); % log of Consumer Price Index for All Urban Consumers: All Items (Index 1982-84=100)
logGDPDEF    = log(GDPDEF); % log of Gross Domestic Product: Implicit Price Deflator (Index 2012=100)
PII          = [NaN; diff(logGDPDEF)]; % Inflation using GDP deflator
PIIy         = logGDPDEF - lagmatrix(logGDPDEF,4);
logGDPCTPI   = log(GDPCTPI); % log of Gross Domestic Product: Chain-type Price Index (Index 2012=100)
logGPDICTPI  = log(GPDICTPI); % log of Gross Private Domestic Investment: Chain-type Price Index (Index 2012=100)
logIPDBS     = log(IPDBS); % log of Business Sector: Implicit Price Deflator (Index 2012=100)
logLP        = log(OPHNFB); % log of nonfarm Business Sector: Real Output Per Hour of All Persons (Index 2012=100)
logINDPRO    = log(INDPRO); % log of Industrial Production Index (Index 2012=100)
logM2REAL    = log(M2REAL); % log of Real M2 Money Stock (Billions of 1982-84 Dollars), deflated by CPI
logTFP       = cumsum(dtfp_util)/400; % log-level of the utilization-adjusted TFP
CASH         = (NCBCDCA + TSDABSNNCB + FDABSNNCB + BOGZ1FL103034000Q/1000)./TABSNNCB; % Cash/Assets as in Bacchetta et al.
CASHT        = (NCBCDCA + TSDABSNNCB + FDABSNNCB + BOGZ1FL103034000Q/1000 + TSABSNNCB)./TABSNNCB; % (Cash+TreasuryHoldings)/Assets. Cash as in Bacchetta et al.
logCASH      = log(CASH); % log of CASH
logCASHT     = log(CASHT); % log of CASHT
logSP500     = log(SP500); % log of stock prices: Standard and Poors' 500



