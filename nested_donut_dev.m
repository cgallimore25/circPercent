%% Sample analysis for nested donut development

load CreditCardData.mat

balnce= data.AMBalance; 
income= data.CustIncome; 
u_rate= data.UtilRate;  

total= sum(balnce); 