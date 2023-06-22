
library(metapoppkg)

par_li1 <- c(
alpha_be1 = 0.09389817,
Beta_be1  = 0.77943975,
mu_be1    = 0.96682356,
Z_be1     = 0.11650832,
D_be1     = 0.93467362,
Td_be1    = 0.80000000,
alpha_af1 = 0.46182906,
Beta_af1  = 0.13459503,
mu_af1    = 0.68640173,
Z_af1     = 0.04444878,
D_af1     = 0.12646338,
Td_af1    = 0.20000000,
theta1    = 0.89626861,
tau1      = 0.91718058,
sigma_SE1 = 0.10000000,
E_01      = 0.44362243,
a_01     = 0.98409838
)
par_nat <- par_li(par_li1,U=5,dir="from_li")
par_li2 <- par_li(par_nat,U=5,dir="to_li")

# see if transforms are inverse
any(abs(par_li1-par_li2)>1e-6)

# another set of test parameters 
mle_test_params <- c(
                     alpha_be=0.093898173,
                     Beta_be=0.779439748,
                     alpha_af=0.461829062,
                     Beta_af=0.134595027,
                     mu_be=0.966823556,
                     Z_be=0.116508318,
                     D_be=0.934673622,
                     mu_af=0.686401734,
                     Z_af=0.044448782,
                     D_af=0.126463381,
                     theta=0.89626861,
                     tau=0.917180577,
                     sigma_SE=0.1348371, # when dt=1/4, size_SE=0.550023483, sigma_SE should be 5*0.1348371
                     Td_be=0.8,
                     Td_af=0.2,
                     E_0=0.44362243,
                     a_0=0.984098384
                   )

# Li et al MLE mapped onto unit scale (allegedly)
li_test_params <- c(
                     alpha_be=0.122449,
                     Beta_be=0.4571429,
                     alpha_af=0.6836735,
                     Beta_af=0.15,
                     mu_be=0.4375,
                     Z_be=0.5633333,
                     D_be=0.49,
                     mu_af=0.2875,
                     Z_af=0.4733333,
                     D_af=0.4366667,
                     theta=0.48,
                     tau=0.1,
                     sigma_SE=0.1, 
                     Td_be=0.8,
                     Td_af=0.2,
                     E_0=0.5,
                     a_0=0.5
                   )
