u1 <- c(    0.999997,
            0.745196,
            0.474838,
            0.330855,
            0.369445,
            0.617118,
            0.64645,
            0.846729,
            0.0413635,
            0.673069,
            0.549574
)
u2 <- c(  0.288161,
          0.100817,
          0.117162,
          0.990211,
          0.654474,
          0.998242,
          0.502823,
          0.918107,
          0.0234147,
          0.596867,
          0.149753)
sum(BiCopDeriv(u1, u2, cop, deriv = "par") /  BiCopPDF(u1, u2, cop))
sum(log(BiCopPDF(u1, u2, cop)))
theta = 50
dat <- c(0.846729, 0.918107)
dat <- c(0.673069, 0.596867)

log(BiCopPDF(u1[8], u2[8], cop))    0.2884429
BiCopPDF(0.846729, 0.918107, cop)
f = (theta*(exp(theta)-1.0)
     *exp(theta*dat[2]+theta*dat[1]+theta)
     )/
    (exp(theta*dat[2]+theta*dat[1])-exp(theta*dat[2]+theta)-exp(theta*dat[1]+theta)+exp(theta))^2.0

t1 = theta*theta;
t2 = exp(theta);
t3 = exp(theta)-1.0;
t5 = theta*v[j];
t6 = theta*u[j];
t8 = exp(theta*(v[j]+u[j]+1));
t10 = exp(theta*(v[j]+u[j]));
t12 = exp(theta*(v[j]+1));
t14 = exp(theta*(u[j]+1));
t15 = t10-t12-t14+t2;
t16 = t15*t15;
out[j] = t1*t3*t8/t16-2.0*theta*t3*t8/t16/t15*(theta*t10-theta*t14);

BiCopDeriv(u1, u2, cop, deriv = "u1") /  BiCopPDF(u1, u2, cop)
BiCopDeriv(u1, u2, cop, deriv = "u2") /  BiCopPDF(u1, u2, cop)
BiCopDeriv(u1, u2, cop, deriv = "par",log = T)
BiCopDeriv(u1, u2, cop, deriv = "par")/  BiCopPDF(u1, u2, cop)
cop <- BiCop(family = 5, par = 1)
BiCopDeriv(0.1, 0.5, cop, deriv = "par") / BiCopPDF(0.1, 0.5, cop)
BiCopDeriv(0.1, 0.5, cop, deriv = "par", log = T)

cop <- BiCop(family = 3, par = 50)
log(BiCopPDF(0.1, 0.5, cop))
BiCopDeriv(0.1, 0.5, cop, deriv = "u2")  /  BiCopPDF(0.1, 0.5, cop)
BiCopDeriv(0.1, 0.5, cop, deriv = "par")  /  BiCopPDF(0.1, 0.5, cop)
