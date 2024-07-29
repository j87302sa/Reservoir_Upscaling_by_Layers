function sig = error_fun(k1,k2,k3,k4,k5,k6)
sig = (((k1.*k4)-(k2.*k3))./(4*(3*(k1+k3).*(k2+k4))+(0.5*(k5+k6)))).*(((2*(k1.*k4-k2.*k3).*(k5+k6))./((k1+k3).*(k2+k4)))+(k5-k6));
end