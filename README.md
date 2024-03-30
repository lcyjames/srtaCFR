# srtaCFR #
srtaCFR (which stands for <ins>**s**</ins>tandized<ins> **r**</ins>eal-<ins>**t**</ins>ime <ins>**a**</ins>djusted <ins>**C**</ins>ase <ins>**F**</ins>atality <ins>**R**</ins>ate) is a package that performs estimation of the standardized real-time fatality rates with adjustment for reporting delay in deaths proposed by Qu and Lee (2024+) (Under Revision).

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
Data <- srtaCFR.SIM(ct = 10000-50*abs(100-c(1:200)),
                    ct_prop_mat = cbind(seq(0.2, 0.6, length.out = 200), 0.2, seq(0.6, 0.2, length.out = 200)),
                    pt_mat = cbind(rep(0.01, 200), rep(0.02, 200), rep(0.06, 200))*replicate(3, exp(c(1:200)*0.004)), seed = 1)

head(Data$ct_mat)
# [1,] 1010 1003 3037
# [2,]  993 1085 3022
# [3,] 1098 1024 3028
# [4,] 1105 1031 3064
# [5,] 1077 1105 3068
# [6,] 1160 1047 3093

head(Data$dt_mat)
# [1,] 0.08705896 0.1233335  1.313139
# [2,] 0.28249618 0.4836344  4.808730
# [3,] 0.51127176 1.0832953  9.966669
# [4,] 0.79844654 1.8896861 16.430421
# [5,] 1.14744447 2.8823098 24.019644
# [6,] 1.56236396 4.0187415 32.432631
```
This data structure is as follows:
`ct_mat` contains the observed confirmed cases in each category, generated via the multinomial distribution given `ct` and `ct_prop_mat`
`dt_mat` contains the observed deaths in each category, according to `pt_mat` and the reporting delay from disease onset to death, controlled by `F_mean` and `F_shape`

<ins>**srtaCFR.EST**</ins>

```
srtaCFR.EST(ct_mat, dt_mat, q_mat=NA, F_mean = 15.43, F_shape = 2.03, maxsteps = 10000)
```
This function computes the srtaCFR as proposed by Qu and Lee (2024). The details of the arguments are as follows:
>- `ct_mat` is the number of confirmed cases in each category across the time points, a matrix of dimension `N` times `J` with `J>2`
>- `dt_mat` is the number of deaths in each category across the time points, a matrix of dimension `N` times `J`
>- `q_mat` is the prespecified distribution used to standardize the fatality rates for each category, a matrix of dimension `N` times `J`, set to an even distribution across the time points as default
>- `F_mean` is the mean of the gamma distribution for the time from disease onset to death
>- `F_shape` is the shape parameter of the gamma distribution for the time from disease onset to death
>- `maxsteps` is an integer specifying the maximum number of steps for the fused lasso to take before termination

Example:
```
Data <- srtaCFR.SIM(ct = 10000-50*abs(100-c(1:200)),
                    ct_prop_mat = cbind(seq(0.2, 0.6, length.out = 200), 0.2, seq(0.6, 0.2, length.out = 200)),
                    pt_mat = cbind(rep(0.01, 200), rep(0.02, 200), rep(0.06, 200))*replicate(3, exp(c(1:200)*0.004)), seed = 1)
srt_fit <-srtaCFR(ct_mat = Data$ct_mat, dt_mat = Data$dt_mat)

round(head(srt_fit$p_gp_spec), 4) #Group-specific fatality rates
# [1,] 0.0119 0.0169 0.0596
# [2,] 0.0060 0.0184 0.0549
# [3,] 0.0046 0.0244 0.0591
# [4,] 0.0100 0.0233 0.0581
# [5,] 0.0074 0.0272 0.0675
# [6,] 0.0121 0.0248 0.0605

round(head(srt_fit$p_hat_std), 4) #Standardized fatality rate
#[1] 0.0295 0.0265 0.0294 0.0304 0.0340 0.0325
```

```
plot(srt_fit$p_gp_spec[,1], ylab = "Group-specific fatality rates", xlab = "Time", type="l",
     lwd = 2, col = "darkgreen", xlim = c(0,200), ylim = c(0,0.15))
lines(srt_fit$p_gp_spec[,2], lty = 2, type = "l", lwd = 2, col = "red")
lines(srt_fit$p_gp_spec[,3], lty = 3, type = "l", lwd = 2, col = "blue")
legend("topleft", col = c("darkgreen", "red", "blue"), lty = c(1, 2, 3), lwd = c(1.5, 1.5, 1.5), bty = "n", x.intersp = 0.5, text.width = 35,
       legend = c(expression(paste("rtaCFR"^"(1)")), expression(paste("rtaCFR"^"(2)")), expression(paste("rtaCFR"^"(3)"))),
       ncol = 3, cex = 1)

```
<img src="https://github.com/lcyjames/srtaCFR/blob/main/example1.png" width="600"/>

```
rt <-rtaCFR.EST(ct = rowSums(Data$ct_mat), dt = rowSums(Data$dt_mat))
plot(srt_fit$p_hat_std, lty = 1, type = "l", lwd = 2, col = "darkgreen", xlab = "Time", ylab = "Fatality rates", xlim = c(0,200), ylim = c(0,0.1))
lines(rt$p_hat, lty = 2, lwd = 2, col = "red")
legend("topleft", col = c("darkgreen","red"), lty = c(1,2), lwd = c(1.5,1.5), bty = "n", x.intersp = 0.5, text.width = 40,
       legend = c(expression(paste("srtaCFR")), expression(paste("rtaCFR (non-stratified)"))), ncol = 3, cex = 1)
```
<img src="https://github.com/lcyjames/srtaCFR/blob/main/example2.png" width="600"/>

# Contact #
Lee Chun Yin, James <<james-chun-yin.lee@polyu.edu.hk>>

# Reference #
Qu, Y., Lee, C. Y., and Lam, K. F. (2022). A novel method to monitor COVID-19 fatality rate in real-time, a key metric to guide public health policy. Scientific Reports, 12(1), 18277. <DOI: [10.1038/s41598-022-23138-4](https://doi.org/10.1038/s41598-022-23138-4)>

Qu, Y., and Lee, C. Y. (2024+). Estimation of standardized real-time fatality rate for ongoing epidemics (Under Revision)
