#### ASSESSMENT ####

obs_count <- select(ebird_split$test, obs = observation_count)

# presence probability is on the complimentary log-log scale
# we can get the inverse link function with
inv_link <- binomial(link = "cloglog")$linkinv

# combine ziplss presence and count predictions
m_nb_pred <- predict(m_nb, ebird_split$test, type = "response") %>% 
    tibble(family = "Negative Binomial", pred = .) %>% 
    bind_cols(obs_count)

m_tw_pred <- predict(m_tw, ebird_split$test, type = "response") %>% 
    tibble(family = "Tweedie", pred = .) %>% 
    bind_cols(obs_count)

# combine predictions from all three models
test_pred <- bind_rows(m_nb_pred, m_tw_pred) %>% 
    mutate(family = as_factor(family))

# spearman’s rank correlation
test_pred %>% 
    group_by(family) %>% 
    summarise(rank_cor = cor.test(obs, pred, 
                                  method = "spearman", 
                                  exact = FALSE)$estimate) %>% 
    ungroup()

#### SELECT THE MODEL ####
# set negative binomial model to be used for predictions
pred_model <- m_tw

## this model merged and it was the best predicting model for the book example as well 

#### use this model to map GHCF relative abundance IN SWG ####

### first we need to bring effort variables into this prediction surface
## to the pediction surface we created  - 
## standard eBird checklist: a 1 km, 1 hour traveling count at the peak time of day for detecting this species##

#### Find the peak time of the day to observe the gray-headed canary flycatcher

# create a dataframe of covariates with a range of start times
seq_tod <- seq(0, 24, length.out = 300)
tod_df <- ebird_split$train %>% 
    # find average pland habitat covariates
    select(starts_with("pland"), elevation_median) %>% 
    summarize_all(mean, na.rm = TRUE) %>% 
    ungroup() %>% 
    # use standard checklist
    mutate(day_of_year = yday(ymd(str_glue("{max_lc_year}-06-15"))),
           duration_minutes = 60,
           effort_distance_km = 1,
           number_observers = 1,
           protocol_type = "Traveling") %>% 
    cbind(time_observations_started = seq_tod)

# predict at different start times
pred_tod <- predict(pred_model, newdata = tod_df, 
                    type = "link", 
                    se.fit = TRUE) %>% 
    as_tibble() %>% 
    # calculate backtransformed confidence limits
    transmute(time_observations_started = seq_tod,
              pred = pred_model$family$linkinv(fit),
              pred_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
              pred_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit))

# find optimal time of day
t_peak <- pred_tod$time_observations_started[which.max(pred_tod$pred_lcl)]

# plot the partial dependence plot
ggplot(pred_tod) +
    aes(x = time_observations_started, y = pred,
        ymin = pred_lcl, ymax = pred_ucl) +
    geom_ribbon(fill = "grey80", alpha = 0.5) +
    geom_line() +
    geom_vline(xintercept = t_peak, color = "blue", linetype = "dashed") +
    labs(x = "Hours since midnight",
         y = "Predicted relative abundance",
         title = "Effect of observation start time on Wood Thrush reporting",
         subtitle = "Peak detectability shown as dashed blue line")

#### assessing covariates ####

# ggplot function
plot_gam <- function(m, title = NULL, ziplss = c("presence", "abundance")) {
    # capture plot
    tmp <- tempfile()
    png(tmp)
    p <- plot(m, pages = 1)
    dev.off()
    unlink(tmp)
    
    # drop addition models in ziplss
    if (m$family$family == "ziplss") {
        is_presence <- map_lgl(p, ~ str_detect(.$ylab, "^s\\.1"))
        if (ziplss == "presence") {
            p <- p[is_presence]  
        } else {
            p <- p[!is_presence]
        }
    }
    
    # extract data
    p_df <- map_df(p, ~ tibble(cov = rep(.$xlab, length(.$x)),
                               x = .$x, fit = .$fit, se = .$se))
    
    # plot
    g <- ggplot(p_df) +
        aes(x = x, y = fit,
            ymin = fit - se, ymax = fit + se) +
        geom_ribbon(fill = "grey80") +
        geom_line(col = "blue") +
        facet_wrap(~ cov, scales = "free_x") +
        labs(x = NULL,
             y = "Smooth function",
             title = title)
    print(g)
    invisible(p_df)
}

plot_gam(m_nb, title = "Tweedie Distribution GAM")

# plot predicted vs. observed
ticks <- c(0, 1, 10, 100, 1000)
mx <- round(max(test_pred$obs))
ggplot(test_pred) +
    aes(x = log10(obs + 1), 
        y = log10(pred + 1)) +
    geom_jitter(alpha = 0.2, height = 0) +
    # y = x line
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    # area where counts off by a factor of 10
    geom_area(data = tibble(x = log10(seq(0, mx - 1) + 1), 
                            y = log10(seq(0, mx - 1) / 10 + 1)),
              mapping = aes(x = x, y = y),
              fill = "red", alpha = 0.2) +
    # loess fit
    geom_smooth(method = "loess", 
                method.args = list(span = 2 / 3, degree = 1)) +
    scale_x_continuous(breaks = log10(ticks + 1), labels = ticks) +
    scale_y_continuous(breaks = log10(ticks + 1), labels = ticks) +
    labs(x = "Observed count",
         y = "Predicted count") +
    facet_wrap(~ family, nrow = 1)
test_pred %>% 
    group_by(family) %>% 
    summarize(n = sum(obs / pred > 10),
              pct = mean(obs / pred > 10))

test_pred %>% 
    group_by(family) %>% 
    summarise(mad = mean(abs(obs - pred), na.rm = TRUE)) %>% 
    ungroup()
