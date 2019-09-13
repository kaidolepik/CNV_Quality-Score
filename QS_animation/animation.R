library(ggplot2)
library(gganimate)
library(gifski)
library(png)
library(dplyr)
library(transformr)

rfalse <- function(n = 1) {
    if (n <= 0)
        return(numeric(0))
    
    u <- runif(n)
    steps <- ifelse(u <= 0.02*25.5, list(data.frame(min = 0, max = 0.02)), list(data.frame(min = 0.02, max = 1))) %>%
        bind_rows()
    
    runif(n, min = steps$min, max = steps$max)
}

rtrue <- function(n = 1) {
    if (n <= 0)
        return(numeric(0))
    
    u <- runif(n)
    steps <- ifelse(u <= 0.02*10.8, list(data.frame(min = 0.98, max = 1)), list(data.frame(min = 0, max = 0.98))) %>%
        bind_rows()
    
    runif(n, min = steps$min, max = steps$max)
}

N <- 33000
freq_ADR <- 0.05
p_controls <- 0.005
OR <- 5
p01 <- 0.008
p10 <- 0.4

n_cases <- round(N * freq_ADR)
n_controls <- N - n_cases
y <- c(rep(1, n_cases), rep(0, n_controls))

p_cases = p_controls * OR / (p_controls * (OR - 1) + 1)
p00 <- 1 - p01
p11 <- 1 - p10
p_cases_pennCNV <- p_cases*p11 + (1 - p_cases)*p01
p_controls_pennCNV <- p_controls*p11 + (1 - p_cases)*p01

x_true <- c(rbinom(n_cases, 1, p_cases), rbinom(n_controls, 1, p_controls))
x_penn <- x_true
x_penn[x_true == 1 & runif(length(x_true)) < p10] <- 0
x_penn[x_true == 0 & runif(length(x_true)) < p01] <- 1
x_qs <- rep(0, length(x_penn))
x_qs[x_true == 1 & x_penn == 1] <- rtrue(sum(x_true == 1 & x_penn == 1))
x_qs[x_true == 0 & x_penn == 1] <- rfalse(sum(x_true == 0 & x_penn == 1))

tmp <- data.frame(y = y, x_true = x_true, x_penn = x_penn, x_qs = x_qs, nr = 1:N) %>%
    mutate(type = ifelse(x_true == x_penn, "True", 
                         ifelse(x_true == 1 & x_penn == 0, "FN",
                                ifelse(x_true == 0 & x_penn == 1, "FP", "False"))))
scaling_factor = N / 0.02

add_state_pos <- function(state) {
    counts <- state %>%
        group_by(x, y) %>%
        summarize(n = n()) %>%
        as.data.frame()
    rownames(counts) = paste0(counts$x, "_", counts$y)
    
    state <- state %>%
        mutate(ix = paste0(x, "_", y),
               x_pos = x + rnorm(N, sd = N * (counts[ix, "n"] / N)^(0.3) / scaling_factor), 
               y_pos = y + rnorm(N, sd = N * (counts[ix, "n"] / N)^(0.3) / scaling_factor)) %>%
        select(-ix)
    
    return(state)
}

state0 <- data.frame(x = rep(0.5, N), y = rep(0.5, N), type = "True", group = tmp$nr, state = "Initialization") %>% add_state_pos()
state1 <- data.frame(x = tmp$x_true, y = tmp$y, type = tmp$type, group = tmp$nr, state = "True CNV vs adverse drug reaction") %>% add_state_pos()
state2 <- data.frame(x = tmp$x_penn, y = tmp$y, type = tmp$type, group = tmp$nr, state = "PennCNV vs adverse drug reaction") %>% add_state_pos()
state3 <- data.frame(x = tmp$x_qs, y = tmp$y, type = tmp$type, group = tmp$nr, state = "CNV QS vs adverse drug reaction") %>% 
    mutate(x_pos = ifelse(x == 0, state2$x_pos, x), y_pos = state2$y_pos)
state_dat <- rbind(state0, state1, state2, state3)

gg <- ggplot(state_dat) +
    geom_point(aes(x = x_pos, y = y_pos, colour = type, group = group), alpha = 0.1) +
    geom_smooth(aes(x = x, y = y), method = lm, se = FALSE, colour = "black") +
    #geom_smooth(aes(x = x, y = y), method = "glm", method.args = list(family = "binomial"), se = FALSE) +
    scale_colour_manual(values = c("black", "dodgerblue3", "firebrick3")) +
    scale_x_continuous(breaks = c(0, 1), limits = c(-0.1, 1.1)) +
    scale_y_continuous(breaks = c(0, 1), limits = c(-0.1, 1.1)) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 16)) +
    coord_fixed() +
    transition_states(state, wrap = FALSE) +
    enter_fade() +
    labs(title = "{closest_state}", x = "Copy number variation", y = "Adverse drug reaction")
anim_save("/Users/kaidolepik/Desktop/Work/PROJECTS/CNVs_baka/QS/QS_animation.gif", gg)
