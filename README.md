
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Installation

<!-- badges: start -->
<!-- badges: end -->

You can install the development version of MaxentDisaggregation from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("simschul/MaxentDisaggregation")
```

Note, this package is under constant development. Together with
co-authors, I’m currently preparing a journal article for describing
more of the background of data disaggregation and showing use cases
within the field of Industrial Ecology.

## Background: Uncertainty propagation involving data disaggregation

The goal of MaxentDisaggregation is an R-package that helps you with
uncertainty propagation when data disaggregation is involved. Data
disaggregation usually involves splitting one data point into several
disaggregates using proxy data. It is a common problem in many different
research disciplines.

```{button-link} https://github.com/jakobsarthur/maxent_disaggregation
:color: info
:expand:
💡 For a **Python** implementation of the package click here!
```

<!-- ![](images/clipboard-133169665.png) -->

```mermaid
flowchart-elk TD
    %% Define node classes
    classDef Aggregate fill:#eeeee4,color:black,stroke:none;
    classDef DisAgg1 fill:#abdbe3,color:black,stroke:none;
    classDef DisAgg2 fill:#e28743,color:black,stroke:none;
    classDef DisAgg3 fill:#abdbe3,color:black,stroke:none;

    agg("Y_0"):::Aggregate
    disagg1("Y_1=x_1 Y_0"):::DisAgg1
    disagg2("Y_2=x_2 Y_0"):::DisAgg1
    disagg3("Y_3=x_3 Y_0"):::DisAgg1
   

    %% Define connections
    agg  --> disagg1
    agg  --> disagg2
    agg  --> disagg3
```

Data disaggregation usually involves an aggregate flow $Y_0$, which is
known, such as the total amount of steel manufactured in a given time
and geography. What we do not know but are interested in are the $K$
disaggregate flows $Y_1,...,Y_K$, such as the different end-use sectors
where the manufactured steel ends up. Even though we do not know the
values of $Y_1, ..., Y_K$, our model structures commonly demands that
the individual $Y_i$’s need to sum to the known aggregate flow $Y_0$ to
respect the mass, energy, stoichiometric or economic balance of the
model

$$
  Y_0 = \sum_{i=1}^{K} Y_i 
$$

This equation, also called an *accounting identity* introduces
dependencies/correlations between the individual disaggregate flows
$Y_i$.

To get estimates for the disaggregate flows, one usually looks for proxy
data. Those proxy data are used to calculate shares (ratios/fractions)
of the respective disaggregate units $x_1, ..., x_K$. To
allocate the entire aggregate flow without leaving any residual (thus to
respect the system balance), those fractions need to sum to one:

$$    \sum_{i=1}^{K} x_i = 1 $$

Disaggregate flows are calculated as

$$ y_i = x_i y_0,  \forall i \in \{1,...,K\}.     $$

## Sampling disaggregates

This package generates a random sample of disaggregates based on the
information provided. The aggregate and the shares are sampled
independently. The distribution from which to sample is determined
internally based on the information provided by the user. This choice of
distribution is mostly based on the principle of Maximum Entropy
(MaxEnt).

The aggregate distribution is determined using the following decision
tree:

<!-- ![](images/clipboard-1519448807.png) -->
```mermaid
flowchart-elk TD
    %% Define node classes
    classDef decision fill:#e28743,color:black,stroke:none;
    classDef distribution fill:#abdbe3,color:black,stroke:none;
    classDef notimplementednode fill:#eeeee4,color:black,stroke:none;

    MeanDecision{{"Best guess/
    mean available?"}}:::decision
    SDDecision{{"Standard deviation available?"}}:::decision
    BoundsDecision1{{"Bounds available?"}}:::decision
    Uniform("Uniform distribution on [a,b]"):::distribution
    GoBackToStart["☠️ !Game Over!
    We suggest to rethink your problem... 🤓"]:::notimplementednode
    BoundsDecision2{{"Bounds available?"}}:::decision
    Normal("Normal distribution"):::distribution
    UnbiasedMean{{"Prefer unbiased mean?"}}:::decision
    TruncNorm("Truncated Normal 
    (Maximum Entropy distribution)"):::distribution
    LogNorm("LogNormal distribution"):::distribution
    LowerBound0{{"Lower bound = 0?"}}:::decision
    Exponential("Exponential"):::distribution
    NotImplemented["Not Implemented"]:::notimplementednode


    %% Define connections
    MeanDecision -- "no" --> BoundsDecision1
    MeanDecision -- "yes" --> SDDecision
    SDDecision -- "yes" --> BoundsDecision2
    BoundsDecision2 -- "yes" --> UnbiasedMean
    UnbiasedMean -- "yes" --> LogNorm
    UnbiasedMean -- "no" --> TruncNorm
    BoundsDecision2 -- "no" --> Normal
    SDDecision -- "no" --> LowerBound0
    LowerBound0 -- "yes" --> Exponential
    LowerBound0 -- "no" --> NotImplemented
    BoundsDecision1 -- "yes" --> Uniform
    BoundsDecision1 -- "no" --> GoBackToStart
```

The shares are sampled from different variants of the Dirichlet
distribution:

<!-- ![](images/clipboard-2622115785.png) -->

```mermaid
flowchart-elk TD
    %% Define node classes
    classDef decision fill:#e28743,color:black,stroke:none;
    classDef distribution fill:#abdbe3,color:black,stroke:none;
    classDef explanationnode fill:#eeeee4,color:black,stroke:none;

    MeanDecision{{"Best guess/mean available?"}}:::decision
    SDDecision{{"Standard deviation available?"}}:::decision
    MaxEntDir("Maximum Entropy Dirichlet"):::distribution
    GenDir("Generalised Dirichlet"):::distribution
    NestedDir("Nested Dirichlet"):::distribution
    UniformDir("Uniform Dirichlet"):::distribution
    

    %% Define connections
    MeanDecision -- "no" --> UniformDir
    MeanDecision -- "yes" --> SDDecision
    MeanDecision -- "paritially" --> NestedDir
    SDDecision -- "no" --> MaxEntDir
    SDDecision -- "yes" --> GenDir
    SDDecision -- "partially" --> NestedDir
```


## How to use

### Sampling disaggregates

The main function is `rdisagg` which creates a random sample of
disaggregates based on the information provided:

``` r
library(MaxentDisaggregation)
#> Loading required package: truncnorm
#> Loading required package: nloptr
#> Loading required package: gtools
#> Loading required package: data.table
#> 
#> Attaching package: 'MaxentDisaggregation'
#> The following object is masked from 'package:gtools':
#> 
#>     rdirichlet
sample <- rdisagg(n = 1000, mean_0 = 100, sd_0 = 5, min = 0, shares = c(0.1, 0.3, 0.6))
head(sample)
#>           [,1]     [,2]     [,3]
#> [1,] 25.542248 15.44779 54.45482
#> [2,]  2.194530 25.33249 79.60015
#> [3,]  1.481192 16.74706 82.64684
#> [4,] 19.181736 43.81459 45.58907
#> [5,] 18.558886 49.36833 36.33028
#> [6,]  2.633181 38.11875 60.79836
```

We can plot the marginal histograms of the sample:

``` r
hist(sample[,1])
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r
hist(sample[,2])
```

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />

``` r
hist(sample[,3])
```

<img src="man/figures/README-unnamed-chunk-3-3.png" width="100%" />

The samples are consistent with all information provided. Thus, summing
the disaggregate samples should provide an aggregate sample consistent
with the information provided (mean: 100, sd: 5):

``` r
sample_agg <- rowSums(sample)
hist(sample_agg)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

And indeed:

``` r
cat('Mean: ', mean(sample_agg), '\n')
#> Mean:  99.91481
cat('SD: ', sd(sample_agg))
#> SD:  5.083025
```

### Sampling aggregates and shares seperatedely

With `MaxentDisaggregation` you can also sample the aggregate and the
shares independently using the `ragg` and `rshares` functions:

``` r
sample_agg <- ragg(1000, mean = 100, sd = 5)
hist(sample_agg)
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

``` r
sample_shares <- rshares(1000, shares = c(0.1, 0.3, 0.6))
boxplot(sample_shares)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />
