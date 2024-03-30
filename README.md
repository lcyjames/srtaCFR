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
Data <- srtaCFR.SIM(ct=10000-50*abs(100-c(1:200)),
                    ct_prop_mat=cbind(seq(0.2,0.6,length.out=200),0.2,seq(0.6,0.2,length.out=200)),
                    pt_mat=cbind(rep(0.01,200),rep(0.02,200),rep(0.06,200))*replicate(3,exp(c(1:200)*0.004)), seed = 1)

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
Data <- srtaCFR.SIM(ct=10000-50*abs(100-c(1:200)),
                    ct_prop_mat=cbind(seq(0.2,0.6,length.out=200),0.2,seq(0.6,0.2,length.out=200)),
                    pt_mat=cbind(rep(0.01,200),rep(0.02,200),rep(0.06,200))*replicate(3,exp(c(1:200)*0.004)), seed = 1)
srt_fit <-srtaCFR(ct_mat = Data$ct_mat, dt_mat = Data$dt_mat)

> round(head(srt_fit$p_gp_spec),4) #Group-specific fatality rates
#[1,] 0.0110 0.0182 0.0566
#[2,] 0.0076 0.0163 0.0624
#[3,] 0.0119 0.0214 0.0543
#[4,] 0.0150 0.0159 0.0550
#[5,] 0.0085 0.0215 0.0580
#[6,] 0.0064 0.0220 0.0565
> round(head(srt_fit$p_hat_std),4) #Standardized fatality rate
#[1] 0.0286 0.0288 0.0292 0.0286 0.0293 0.0283
```

```
plot(srt_fit$p_gp_spec[,1], ylab="Group-specific fatality rates", xlab="Time", type="l", #xaxt="n",
     lwd=2, col="darkgreen", xlim=c(0,200), ylim=c(0,0.15))
lines(srt_fit$p_gp_spec[,2], lty=2, type="l", lwd=2, col="red")
lines(srt_fit$p_gp_spec[,3], lty=3, type="l", lwd=2, col="blue")
legend("topleft",col=c("darkgreen","red","blue"),lty=c(1,2,3), lwd = c(1.5,1.5,1.5),bty = "n",x.intersp=0.5,text.width=35,
       legend = c(expression(paste("rtaCFR"^"(1)")), expression(paste("rtaCFR"^"(2)")), expression(paste("rtaCFR"^"(3)"))),
       ncol = 3,cex=1)

```
<img src="https://github.com/lcyjames/srtaCFR/blob/main/example1.png" width="600"/>

```
rt <-rtaCFR.EST(ct = rowSums(Data$ct_mat), dt = rowSums(Data$dt_mat))
plot(srt_fit$p_hat_std, lty=1, type="l", lwd=2, col="darkgreen", xlab="Time", ylab="Fatality rates",
     xlim=c(0,200), ylim=c(0,0.1))
lines(rt$p_hat, lty=2, lwd=2, col="red")
legend("topleft",col=c("darkgreen","red"),lty=c(1,2),
       lwd = c(1.5,1.5),bty = "n",x.intersp=0.5,text.width=40,
       legend = c(expression(paste("srtaCFR")),
                  expression(paste("rtaCFR (non-stratified)"))),ncol = 3,cex=1)
```
<img src="https://github.com/lcyjames/srtaCFR/blob/main/example2.png" width="600"/>

# Contact #
Lee Chun Yin, James <<james-chun-yin.lee@polyu.edu.hk>>

# Reference #
Qu, Y., Lee, C. Y., and Lam, K. F. (2022). A novel method to monitor COVID-19 fatality rate in real-time, a key metric to guide public health policy. Scientific Reports, 12(1), 18277. <DOI: [10.1038/s41598-022-23138-4](https://doi.org/10.1038/s41598-022-23138-4)>
