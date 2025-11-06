function wealthtax=BM2022_WealthTaxRevenue(h,aprime,a,z,tau_a,xi_a)

wealthtax=a-((1-tau_a)/(1-xi_a))*(a^(1-xi_a));


end