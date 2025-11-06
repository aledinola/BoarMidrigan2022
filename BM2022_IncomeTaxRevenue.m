function incometax=BM2022_IncomeTaxRevenue(h,aprime,a,z,r,tau,xi,iota,delta,alpha)

w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));

pretaxincome=w*h*z+r*a;

incometax=pretaxincome-((1-tau)/(1-xi))*(pretaxincome^(1-xi)) - iota;

end