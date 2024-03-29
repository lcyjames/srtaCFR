# srtaCFR #
srtaCFR (which stands for <ins>**s**</ins>tandized<ins> **r**</ins>eal-<ins>**t**</ins>ime <ins>**a**</ins>djusted <ins>**C**</ins>ase <ins>**F**</ins>atality <ins>**R**</ins>ate) is a package that performs estimation of the standardized real-time fatality rates with adjustment for reporting delay in deaths proposed by Qu and Lee (2024+) (Under revision).

**srtaCFR** relies on the R-packages `genlasso` and `Rtools`, which are hosted on CRAN.

# How to import the Functions #
> install.packages("devtools")<br />
> library(devtools) <br /> 
> source_url("https://github.com/lcyjames/srtaCFR/blob/main/CoreFunctions.R?raw=TRUE")

# Usage #
The package contains 3 functions:
|Functions  | Description|
|------------- | -------------|
rtaCFR.EST  | Computation of the rtaCFR as proposed in Qu et al. (2022); usage provided in https://github.com/lcyjames/rtaCFR
srtaCFR.SIM  | Generate a data set according to the simulation study in Qu and Lee (2024+)
srtaCFR.EST  | Computation of the srtaCFR as proposed in Qu and Lee (2024+)

<ins>**srtaCFR.SIM**</ins>

```
srtaCFR.SIM(ct, ct_prop_mat, pt_mat, seed = NA, F_mean = 15.43, F_shape = 2.03)
```
This function generates a data set according to the model in Qu and Lee (2024+) that takes the arguments:
>- `ct` is the number of confirmed cases in the population across the `N` time points, a vector of length `N`
>- `ct_prop_mat` is the proportion of confirmed cases in each category across the time points, a matrix of dimension `N` times `J`; `J` is the total number of categories
>- `pt_mat` is the proportion of confirmed cases who will eventually die from the disease in each category, a matrix of dimension `N` times `J`
>- `F_mean` is the mean of the gamma distribution for the time from disease onset to death
>- `F_shape` is the shape parameter of the gamma distribution for the time from disease onset to death

Take scenario III in the simulation study in Qu and Lee (2024+) as an example:
```
Data <- srtaCFR.SIM(ct=10000-50*abs(100-c(1:200)), ct_prop_mat=cbind(seq(0.2,0.6,length.out=200),0.2,seq(0.6,0.2,length.out=200)),
          pt_mat=cbind(rep(0.01,200),rep(0.02,200),rep(0.06,200))*replicate(3,exp(c(1:200)*0.004)), seed = 1234)

head(Data$ct_mat)
#[1,] 1036  974 3040
#[2,] 1012 1020 3068
#[3,] 1057 1027 3066
#[4,] 1101 1106 2993
#[5,] 1073 1087 3090
#[6,] 1111 1004 3185

head(Data$dt_mat)
#[1,] 0.07254913 0.1305884  1.342159
#[2,] 0.28619788 0.4672737  4.953680
#[3,] 0.60777737 0.9480578 10.409497
#[4,] 1.00222024 1.6281263 17.387538
#[5,] 1.42538446 2.5242256 25.479600
#[6,] 1.86519333 3.5626964 34.244764
```
`ct_mat` contains the observed number of confirmed cases in each category, generated based on the multinomial distribution given `ct` and `ct_prop_mat`
`dt_mat` contains the observed number of deaths in each category, according to `pt_mat` and the reporting delay from disease onset to death (`F_mean' and `F_shape`)

This data structure is as follows:
>- `ct` is the number of confirmed cases
>- `dt` is the number of deaths with reporting delay from disease onset to death

<ins>**rtaCFR.EST**</ins>

```
rtaCFR.EST(ct, dt, F_mean = 15.43, F_shape = 2.03, maxsteps = 10000)
```
This function computes the rtaCFR as proposed in Qu et al. (2022). The details of the arguments are as follows:
>- `ct` is the number of confirmed cases
>- `dt` is the number of deaths
>- `F_mean` is the mean of the gamma distribution for the time from disease onset to death
>- `F_shape` is the shape parameter of the gamma distribution for the time from disease onset to death
>- `maxsteps` is an integer specifying the maximum number of steps for the fused lasso to take before termination

Example:
```
Data <- rtaCFR.SIM(ct = 3000-5*abs(100-c(1:200)), pt = 0.01*exp(0.012*c(1:200)), seed = 1)
rt_fit <- rtaCFR.EST(ct = Data$ct, dt = Data$dt)

round(head(rt_fit$p_hat),4)
# [1] 0.0088 0.0096 0.0107 0.0131 0.0087 0.0134
round(tail(rt_fit$p_hat),4)
# [1] 0.1030 0.1131 0.1018 0.1135 0.1102 0.1024
```

```
plot(rt_fit$p_hat, type="b", pch = 19, ylab = "Fatality rates", xlab = "Time", col = "red", cex = 0.6)
lines(c(1:200), 0.01*exp(0.012*c(1:200)), lwd = 2)
legend("topleft", legend = c("rtaCFR", "true"), col = c("red", "black"), lty = c(1,1),lwd = c(1:2), pch = c(19,NA), cex = 0.8)

```
<img src="https://github.com/lcyjames/rtaCFR/blob/main/example.png" width="600"/>

# Contact #
Lee Chun Yin, James <<james-chun-yin.lee@polyu.edu.hk>>

# Reference #
Qu, Y., Lee, C. Y., and Lam, K. F. (2022). A novel method to monitor COVID-19 fatality rate in real-time, a key metric to guide public health policy. Scientific Reports, 12(1), 18277. <DOI: [10.1038/s41598-022-23138-4](https://doi.org/10.1038/s41598-022-23138-4)>
