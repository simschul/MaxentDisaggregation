library(matrixStats)
library(ggplot2)
library(disaggR)

shares <- c(0.1, 0.4, 0.5)
sds <- c(0.2, 0.3, 0.1)
sds[(sds/shares) > 1] <- shares[(sds/shares) > 1]

rdirichlet_generalised_fixed <- function(n, shares, sds, fix = FALSE) {
  sds[(sds/shares) > 0.5] <- shares[(sds/shares) > 0.5] * 0.5
  rdirichlet_generalised(n, shares, sds)
}

#TODO: check before return samples: rel deviation is larger ...., one possible reason,
# is that sds are too large. try with SDs < shares (CV < 1)

temp_fun <- function(K, max_cv) {
  shares <- rdirichlet_uniform(1, K) %>% as.numeric
  sds <- runif(K, min = 0, max = max_cv) * shares
  sample <- rdirichlet_generalised(1E5, shares, sds)
  sample_fixed <- rdirichlet_generalised_fixed(1E5, shares, sds)
  dt <- data.table(
    K = K,
    max_cv = max_cv,
    sample_mean = colMeans2(sample),
    desired_mean = shares,
    sample_sd = colSds(sample),
    desired_sd = sds,
    fixed_mean = colMeans2(sample_fixed),
    fixed_sd = colSds(sample_fixed)
  )
  return(dt)
}


max_cvs <- c(0.5, 0.6, 0.75, 1)
out_list <- vector('list', length = length(max_cvs))
for (i in max_cvs) {
  cat(i, '')
  out <- lapply(1:50, function(x) temp_fun(K = 10, max_cv = i))
  out_dt <- rbindlist(out, idcol = 'iteration')
  out_list[[as.character(i)]] <- out_dt
}

dt <- rbindlist(out_list)
dt[, reldiff := (sample_mean - desired_mean) / desired_mean]
dt[, reldiff_fixed := (fixed_mean - desired_mean) / fixed_mean]

dt[, .(mean = mean(reldiff),
       p0.1 = quantile(reldiff, probs = c(0.1)),
       p0.9 = quantile(reldiff, probs = c(0.9))

       ), by = max_cv]

ggplot(dt, aes(x = reldiff)) +
  geom_histogram(alpha = 0.7, position = 'identity') +
   geom_histogram(aes(x = reldiff_fixed), alpha = 0.2, position = 'identity',
                  fill = 'red', col = 'red') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~max_cv, labeller = 'label_both', scales = 'free') +
  xlab("reldiff = (sample_mean - desired_mean) / desired_mean") +
  scale_color_viridis_d() +
  theme(legend.position = 'none')


ggplot(dt, aes(x = sample_mean, y = desired_mean, col = as.factor(iteration))) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_abline(slope = 1, intercept = 0, col = 'red') +
  facet_wrap(~max_cv, labeller = 'label_both') +
  scale_color_viridis_d() +
  scale_x_log10() +
  scale_y_log10() +
  theme(legend.position = 'none')

ggplot(dt, aes(x = sample_sd, y = desired_sd, col = as.factor(iteration))) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_abline(slope = 1, intercept = 0, col = 'red') +
  facet_wrap(~max_cv, labeller = 'label_both', scales = 'free') +
  scale_color_viridis_d() +
  theme(legend.position = 'none')




## version 2 with a fixed CV

temp_fun2 <- function(K, cv) {
  shares <- rdirichlet_uniform(1, K) %>% as.numeric
  sds <- cv * shares
  sample <- rdirichlet_generalised(1E5, shares, sds)
  dt <- data.table(
    K = K,
    cv = cv,
    sample_mean = colMeans2(sample),
    desired_mean = shares,
    sample_sd = colSds(sample),
    desired_sd = sds
  )
  return(dt)
}



out_list2 <- list()
for (i in c(0.5, 1,2,3,4)) {
  out <- lapply(1:50, function(x) temp_fun2(K = 10, cv = i))
  out_dt <- rbindlist(out, idcol = 'iteration')
  out_list2[[as.character(i)]] <- out_dt
}

dt2 <- rbindlist(out_list2)
dt2[, reldiff := (sample_mean - desired_mean) / desired_mean]


ggplot(dt2, aes(x = reldiff)) +
  geom_histogram(alpha = 0.3, position = 'identity') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  facet_wrap(~cv, labeller = 'label_both', scales = 'free') +
  xlab("reldiff = (sample_mean - desired_mean) / desired_mean") +
  scale_color_viridis_d() +
  theme(legend.position = 'none')



gamma = 2
small_value = 0.000001
shares = c(small_value, 0.3, 0.5, 0.2-small_value)
samples <- rdirichlet_standard(1E5, alpha = gamma * shares)
samples <- gtools::rdirichlet(1E5, alpha = gamma * shares)
samples[samples[,1] == 0]
samples %>%
  as.data.table() %>%
  melt %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 200) +
  facet_wrap(~variable, scales = 'free')

means <- colMeans2(samples)
means
sds <- colSds(samples)
sds
sds/means
colMaxs(samples)
