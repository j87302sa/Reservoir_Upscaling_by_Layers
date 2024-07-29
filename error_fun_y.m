function error_y = error_fun_y(k1,k2,k3,k4,k5,k6)
error_y = ((k1.*k4 - k2.*k3).^2)./(((k1+k3).*(k2+k4))+0.333.*(k5+k6));
end