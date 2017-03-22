r <- stanc(file = "copula.stan")
r <- stan_model(file = "copula.stan")
writeLines(r$cppcode, "outfile.txt")

expose_stan_functions(r)

sm <- stan_model(model_code = r$cppcode  )


stan_shared("outfile.cpp", "foo.so")
>> setting environment variables:
    PKG_LIBS =  -L'/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/StanHeaders/lib' -lStanHeaders
PKG_CPPFLAGS =   -I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/Rcpp/include/"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/RcppEigen/include/"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/RcppEigen/include/unsupported"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/rstan/include/boost_not_in_BH"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/BH/include"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/StanHeaders/include/src/"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/StanHeaders/include/"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -DBOOST_NO_CXX11_RVALUE_REFERENCES
g++ -I/usr/share/R/include -DNDEBUG
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/Rcpp/include/"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/RcppEigen/include/"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/RcppEigen/include/unsupported"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/rstan/include/boost_not_in_BH"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/BH/include"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/StanHeaders/include/src/"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/StanHeaders/include/"
-I"/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -DBOOST_NO_CXX11_RVALUE_REFERENCES     -fpic  -O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -c /tmp/Rtmpc05UNa/file5e0b2c05ff7c.cpp -o /tmp/Rtmpc05UNa/file5e0b2c05ff7c.o
g++ -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o foo.so /tmp/Rtmpc05UNa/file5e0b2c05ff7c.o -L/home/hoanguc3m/R/x86_64-pc-linux-gnu-library/3.3/StanHeaders/lib -lStanHeaders -L/usr/lib/R/lib -lR
