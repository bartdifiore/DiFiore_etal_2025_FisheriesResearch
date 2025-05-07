#--------------------------------------
## Setup
#--------------------------------------

library(tidyverse)
library(emmeans)
library(ggeffects)
library(glmmTMB)
library(DHARMa)
library(gt)
source("Code/theme.R")


#--------------------------------------
## Get data
#--------------------------------------

df <- read.csv("Data/Rawdata_20240118.csv") %>% 
  janitor::clean_names() %>%
  select(trap_pair, buoy_number, bait_type, date_set, date_retrieved, deployment_nights, lobster, num_short, num_legalsize, num_vnotch, num_egged, total_crab, rock_crab, green_crab, bait_remaining) %>%
  mutate(bait_type = case_when(bait_type == "Herring" ~ "HERRING", 
                               bait_type == "50:50:00" ~ "50:50",
                               bait_type == "STUCK" ~ NA,
                               .default = bait_type),
         bait_type = forcats::fct_reorder(bait_type, lobster, .fun = mean, .na_rm = T, .desc = T),
         alt_trad = ifelse(bait_type != "HERRING", "b_alt", "a_herring"), 
         bait_remaining = bait_remaining / 100) %>%
  drop_na(alt_trad) %>%
  filter(deployment_nights <= 10)

#----------------------------------------
## Summary stats for manuscript
#----------------------------------------

dim(df)
sum(df$lobster, na.rm = T)
sum(df$rock_crab, na.rm =T)

df %>%
  group_by(alt_trad) %>%
  summarize(total = sum(lobster, na.rm = T), 
            n = n()) %>%
  mutate(cpue = total / n)


df %>% 
  group_by(bait_type) %>%
  drop_na(lobster) %>%
  summarize(sum_lobster = sum(lobster, na.rm = T), 
            n = n(), 
            mean = mean(lobster, na.rm = T),
            se = sd(lobster, na.rm = T)/n) %>% 
  mutate(cpue = sum_lobster/n,
         catch_relative_to_herring = cpue/(1279/300))


df %>% 
  group_by(bait_type) %>%
  drop_na(lobster) %>% 
  summarize(min = min(deployment_nights), 
            max = max(deployment_nights)) %>% 
  mutate(soak_time = paste(min, max, sep = "-"))

table_1 <- df %>% 
  drop_na(lobster) %>%
  select(-c(bait_remaining, total_crab, green_crab)) %>%
  rename("All lobster" = lobster, 
         "Legal lobster" = num_legalsize, 
         "Sub-legal lobster" = num_short, 
         "V-notched lobster" = num_vnotch, 
         "Egged lobster" = num_egged, 
         "Cancer spp." = rock_crab) %>%
  pivot_longer(cols = c(`All lobster`:`Cancer spp.`)) %>%
  mutate(name = factor(name, levels = c("All lobster","Legal lobster", "Sub-legal lobster", "V-notched lobster", "Egged lobster", "Cancer spp." ))) %>%
  group_by(bait_type, name) %>%
  summarize(sum = sum(value, na.rm = T), 
            n = n(), 
            mean = mean(value, na.rm = T),
            se = sd(value, na.rm = T)/n) %>%
  ungroup() %>%
  mutate(mean.ch = as.character(round(mean, 2)), 
         se.ch = as.character(round(se, 3)), 
         value = paste(mean.ch, se.ch, sep = " ± ")) %>%
  select(name, bait_type, value) %>%
  pivot_wider(values_from = value, names_from = bait_type) |>
  ungroup() |>
  gt(rowname_col = "name") |>
  gt::rows_add(
    .before = 1,
    name = "Sample size",
    HERRING = "298" ,
    `50:50` = "168",
    `100T` = "46",
    `100F` = "36"
  )|>
  gt::rows_add(
    .before = 2,
    name = "Soak time range (nights)",
    HERRING = "3-8" ,
    `50:50` = "3-8",
    `100T` = "3-6",
    `100F` = "3-7"
  ) |>
  gt::rows_add(
    .before = 1,
    name = "Bait description",
    HERRING = "Traditional bait" ,
    `50:50` = "50:50% blend of TVP soaked in hydrolysate and gurry",
    `100T` = "100% TVP soaked in hydrolysate",
    `100F` = "100% gurry"
  ) |>
  cols_label(
    HERRING = md("**Herring**"),
    `50:50` = md("**50:50**"),
    `100T` = md("**100T**"), 
    `100F` = md("**100F**") 
  ) |>
  tab_row_group(
    label = md("**Summary**"),
    id = "1", 
    rows = 1:3) |>
  tab_row_group(
    label = md("**CPUE** (ind. per soak)"),
    id = "2",
    rows = 4:9) |>
  row_group_order(groups = c("1", "2")) |>
  cols_width(everything() ~ px(150)) |>
  cols_align(
    align = "center",
    columns = everything()
  )

table_1 |> gtsave("Figures/tab_1.docx")

df %>% 
  group_by(bait_type) %>%
  drop_na(lobster) %>%
  summarize(sum_lobster = sum(lobster, na.rm = T), 
            n = n(), 
            mean = mean(lobster, na.rm = T),
            se = sd(lobster, na.rm = T)/n) %>% 
  mutate(cpue = sum_lobster/n,
         catch_relative_to_herring = cpue/(1279/307))

df %>% 
  group_by(alt_trad) %>% 
  summarize(sum_lobster = sum(lobster, na.rm = T), 
            n = n()) %>% 
  mutate(cpue = sum_lobster/n, 
         catch_relative_to_herring = cpue/(1279/307))


df %>% 
  drop_na(lobster) %>%
  select(-c(bait_remaining, total_crab, green_crab)) %>%
  rename("All lobster" = lobster, 
         "Legal lobster" = num_legalsize, 
         "Sub-legal lobster" = num_short, 
         "V-notched lobster" = num_vnotch, 
         "Egged lobster" = num_egged, 
         "Cancer spp." = rock_crab) %>%
  pivot_longer(cols = c(`All lobster`:`Cancer spp.`)) %>%
  mutate(name = factor(name, levels = c("All lobster","Legal lobster", "Sub-legal lobster", "V-notched lobster", "Egged lobster", "Cancer spp." ))) %>%
  group_by(name) %>%
  summarize(sum = sum(value, na.rm = T), 
            n = n(), 
            mean = mean(value, na.rm = T),
            se = sd(value, na.rm = T)/n)


#----------------------------------------
## Analysis of 2023 trials
#----------------------------------------
df_long <- df %>% 
  select(-alt_trad, -total_crab, -green_crab) %>% 
  pivot_longer(cols = c(lobster:rock_crab))

# Full model
mod1 <- glmmTMB(value ~ bait_type*name*deployment_nights + (1|trap_pair), data = df_long, family = nbinom2(link = "log"), ziformula = ~1)
summary(mod1)
simulated_residuals <- DHARMa::simulateResiduals(mod1)
plot(simulated_residuals)

# Model selection
mod1b <- glmmTMB(value ~ bait_type*name + deployment_nights + (1|trap_pair), data = df_long, family = nbinom2(link = "log"), ziformula = ~1)
summary(mod1b)

sjPlot::tab_model(mod1b)

anova(mod1, mod1b) # Little evidence that the interaction improves the model.

post_hoc2 <- emmeans(mod1b, ~bait_type | name, condition = "response")

pairs(post_hoc2, reverse = T)
temp_cont <- contrast(post_hoc2, method= "revpairwise", type = "response", adjust = "bonferroni")
# temp_cont <- contrast(post_hoc2, method= "revpairwise", adjust = "bonferroni")

data.frame(temp_cont) %>%
  mutate(name = case_when(name == "lobster" ~ "All lobster",
                          name == "num_legalsize" ~ "Legal lobster", 
                          name == "num_short" ~ "Sub-legal lobster", 
                          name == "num_vnotch" ~ "V-notched lobster", 
                          name == "num_egged" ~ "Egged lobster", 
                          name == "rock_crab" ~ "Cancer spp.")) %>%
  group_by(name) %>%
  select(-df, -null) %>%
  mutate(across(c(ratio, SE, z.ratio), \(x) round(x,2))) %>%
  pivot_wider(names_from = name,
              values_from = c(ratio, SE, z.ratio, p.value),
              names_glue = "{name}_{.value}") %>%
  select(sort(tidyselect::peek_vars())) %>%
  gt(rowname_col = "contrast", groupname_col = "name") %>%
  gt::tab_spanner_delim(delim = "_") %>%
  # gt::tab_spanner(label = "test", columns = c(1:3))
  fmt_number(columns = matches(c("z.ratio", "SE", "ratio")), n_sigfig = 2) %>%
  fmt_scientific(columns = matches("p.value"),
    exp_style = "e1",
    force_sign_n = TRUE
  ) %>%
gtsave("Figures/tab_2.docx")

confint(post_hoc2, type = "response")

predictions <- ggeffects::ggpredict(mod1b, terms = ~bait_type*name, condition = c(deployment_nights = 5))
plot(predictions, add.data = F)
pred_df <- as.data.frame(predictions) %>%
  mutate(conf.high = ifelse(is.infinite(conf.high), 0, conf.high), 
         conf.low = ifelse(conf.low < 1e-15, 0, conf.low))

cols_ramp <- colorRampPalette(c("#F7C7BD", "#CD3718"))
alt_cols <- rev(cols_ramp(3))

# 
# p3 <- df %>%
#   ggplot(aes(x = bait_type, y = lobster))+
#   geom_jitter(aes(fill = bait_type), show.legend = F, height = 0, width = 0.2, alpha = 0.5, pch = 21, color = "grey", size = 1.9)+
#   scale_fill_manual(values = c("#18AFCD", alt_cols))+
#   geom_pointrange(data = pred_df, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high))+
#   annotate(geom = "text", x = 1:4, y = c(15, 11, 7, 7), label = c("A", "B", "B", "C"))+
#   labs(y = expression(paste("Lobster catch (ind. ", haul^-1, ")")), x = "")+
#   theme_bd()

pred_df %>%
  ggplot()+
  geom_point(aes(x = x, y = predicted, color = x), size = 5, show.legend = F)+
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = x), width = 0.5, show.legend = F)+
  scale_color_manual(values = c("#18AFCD", alt_cols))+
  labs(y = expression(paste("Lobster catch (ind. ", haul^-1, ")")), x = "")+
  facet_wrap(~group, scales = "free")+
  theme_bd()

ggsave("Figures/figure_1.svg", width = 9.5, height = 7)

#---------------------------
## Figure for Reviewer 2
#---------------------------


rev2 <- df %>% 
  select(trap_pair, date_retrieved, bait_type, lobster) %>%
  mutate(id = 1:n()) %>%
  group_by(trap_pair, date_retrieved, bait_type) %>%
  summarize(total_catch = sum(lobster, na.rm = T)) %>%
  pivot_wider(names_from = bait_type, values_from = total_catch) %>%
  pivot_longer(cols = c(`50:50`, `100F`, `100T`)) %>% 
  drop_na(value)

rev2 %>% 
  mutate(diff = value - HERRING) %>% 
  ggplot()+
  geom_histogram(aes(x = diff))+
  facet_wrap(~name)

rev2 %>%
  group_by(HERRING, value) %>%
  mutate(freq = n()) %>%
  ggplot(aes(x = HERRING, y = value))+
  geom_point(aes(color = name, size = freq))+
  geom_abline(slope = 1, intercept = 0, linetype = 4)+
  geom_smooth(aes(color= name), method = "lm" ) +
  labs(x = "CPUE (ind. per soak) w/ Herring", y = "CPUE (ind. per soak) w/ alternative")+
  theme_classic()

rev2 %>%
  group_by(HERRING, value) %>%
  mutate(freq = n()) %>%
  ggplot(aes(x = HERRING, y = value))+
  geom_jitter(aes(fill = name), width = 0.2, height = 0, pch = 21, size = 2.5, alpha = 0.75, color = "black")+
  geom_abline(slope = 1, intercept = 0, linetype = 4)+
  geom_smooth(aes(color= name), method = "lm" ) +
  scale_color_manual(values = rev(alt_cols))+
  scale_fill_manual(values = rev(alt_cols))+
  labs(x = "All lobster CPUE\n(ind. per soak) w/ Herring", y = "All lobster CPUE\n(ind. per soak) w/ alternative", color = "Alternative\nbait type", fill = "Alternative\nbait type")+
  theme_classic()
ggsave("Figures/figure_reviewer2.png", width = 8, height = 4, dpi = 600)

lm0<- lm(value ~ HERRING *name, rev2)
summary(lm0)


lm1<- MASS::glm.nb(value ~ HERRING *name, rev2)
summary(lm1)

plot(ggeffects::ggpredict(lm1, terms = ~HERRING*name), add.data = T)


#--------------------------------------
# Deployment length analysis
#--------------------------------------

mod3 <- glmmTMB(lobster ~ 0 +bait_type*deployment_nights + (1|trap_pair), data = df, family = nbinom2(link = "log"), ziformula = ~1)
summary(mod3)
simulated_residuals <- DHARMa::simulateResiduals(mod3)
plot(simulated_residuals)

sjPlot::tab_model(mod3, file = "Figures/mod3_table.html")

plot(ggpredict(mod3, terms = ~ deployment_nights*bait_type), add.data = T)

post_hoc <- emtrends(mod3, ~ bait_type, var = "deployment_nights", adjust = "bonferroni")

contrast(post_hoc)
confint(post_hoc)

pairs(post_hoc, type = "response", adjust = "bonferroni") %>%
  data.frame() %>%
  select(-df) %>%
  mutate(across(c(estimate, SE, z.ratio), \(x) round(x,2))) %>%
  gt() %>%
  fmt_number(columns = matches(c("z.ratio", "SE", "ratio")), n_sigfig = 2) %>%
  fmt_scientific(columns = matches("p.value"),
                 exp_style = "e1",
                 force_sign_n = TRUE
  ) %>%
  gtsave("Figures/tab_3.docx")


# mod3b <- glmmTMB(lobster ~ bait_type + deployment_nights + (1|trap_pair), data = df, family = nbinom2(link = "log"))
# summary(mod3b)
# plot(ggpredict(mod3b, terms = ~ deployment_nights), add.data = T)
# 
# mod3c <- glmmTMB(lobster ~ bait_type + (1|trap_pair), data = df, family = nbinom2(link = "log"))
# 
# anova(mod3, mod3b)
# anova(mod3, mod3c)
# anova(mod3b, mod3c)


predictions <- ggeffects::ggpredict(mod3, terms = ~deployment_nights*bait_type)
plot(predictions, add.data = T)

pred_df <- as.data.frame(predictions)

p4.b <- df %>% 
  ggplot(aes(x = deployment_nights, y = lobster))+
  geom_jitter(aes(fill = alt_trad), height = 0.1, width = 0.1, alpha = 0.5, pch = 21, color = "grey", size = 1.9, show.legend = F)+
  scale_fill_manual(values = c("#18AFCD", alt_cols))+
  geom_ribbon(data = pred_df, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group), alpha = 0.2, show.legend = F)+
  geom_line(data = pred_df, aes(x = x, y = predicted, color = group), linewidth = 1.5, show.legend = F)+
  scale_color_manual(values = c("#18AFCD", alt_cols))+
  labs(y = expression(paste("Catch per unit effort (ind. ", soak^-1, ")")), x = "Length of deployment (days)", color = "Bait type")+
  theme_bd()+
  theme(legend.position = c(0.8,0.85), legend.background = element_blank())


dummy <- data.frame(distinct(df, bait_type), deployment_nights = 0, mean_remaining = 100, n = 10, se = NA)


p4.a <- df %>% 
  group_by(bait_type, deployment_nights) %>% 
  summarize(mean_remaining = mean(bait_remaining, na.rm = T)*100, 
            n = n(),
            se = mean_remaining/n) %>%
  bind_rows(dummy) %>%
  filter(deployment_nights < 10) %>%
  ggplot(aes(x = deployment_nights, y = mean_remaining, group = bait_type))+
  geom_pointrange(aes(ymin = mean_remaining - se, ymax = mean_remaining + se, color = bait_type))+
  geom_line(aes(color = bait_type), linewidth = 1.0)+
  scale_color_manual(values = c("#18AFCD", alt_cols))+
  labs(x = "Length of deployment (nights)", y = "Bait remaining (%)", color = "Bait type")+
  theme_bd()+
  theme(legend.position = c(0.8, 0.8))

cowplot::plot_grid(p4.a, p4.b, align = "h", nrow = 1, labels = "AUTO")

ggsave("Figures/figure_2.png", width = 8, height = 4, dpi = 600)


#----------------------------------------
## Temperature analysis
#----------------------------------------

df_2023_temp <- read.csv("Data/Fieldtrials_2023_wtemp.csv") %>% as_tibble() %>%
  drop_na(avg_temp_c) %>% 
  filter(deployment_nights <= 10) %>%
  filter(bait_type %in% c("HERRING", "50:50"))

df_2023_temp %>% 
  group_by(bait_type) %>% 
  count()


mod_temp1 <- glmmTMB(lobster ~ bait_type*avg_temp_c + deployment_nights + (1|trap_pair), data = df_2023_temp, family = nbinom2(link = "log"), ziformula = ~1)
summary(mod_temp1)

# mod_temp2 <- glmmTMB(lobster ~ alt_trad*avg_temp_c + deployment_nights + (1|trap_pair), data = df_2023_temp, family = nbinom2(link = "log"), ziformula = ~1)
# summary(mod_temp2)

# AIC(mod_temp1, mod_temp2)


post_hoc <- emtrends(mod_temp1, ~ bait_type, var = "avg_temp_c", adjust = "bonferroni")

contrast(post_hoc)
confint(post_hoc)

pairs(post_hoc, type = "response", adjust = "bonferroni") 


predictions <- ggeffects::ggpredict(mod_temp1, terms = ~avg_temp_c*bait_type, condition = c(deployment_time = 5))
plot(predictions, add.data = T) # Some evidence, albeit weak that the alternative fished better at warmer temperatures. 
pred_df <- as.data.frame(predictions)


temp_p2 <- df_2023_temp %>% 
  ggplot(aes(x = avg_temp_c, y = lobster))+
  geom_jitter(aes(fill = bait_type), height = 0.1, width = 0.1, alpha = 0.5, pch = 21, color = "grey", size = 1.9, show.legend = F)+
  scale_fill_manual(values = c("#CD3718", "#18AFCD"))+
  geom_ribbon(data = pred_df, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group), alpha = 0.2, show.legend = F)+
  geom_line(data = pred_df, aes(x = x, y = predicted, color = group), linewidth = 1.5)+
  scale_color_manual(values = c("#CD3718", "#18AFCD"))+
  labs(y = expression(paste("Catch per unit effort (ind. ", soak^-1, ")")), x = expression(paste("Temperature (", degree, "C)")), color = "Bait type")+
  theme_bd()+
  theme(legend.position = c(0.8, 0.8), legend.background = element_blank())

source("Code/pull_temp_data.R")

cowplot::plot_grid(temp_plot, temp_p2, nrow = 2, align = "v", labels = "AUTO")

ggsave("Figures/figure3.png", width = 8, height = 7.5, dpi = 600)

predict(mod_temp1, newdata = data.frame(bait_type = "50:50", deployment_nights = 5, avg_temp_c = c(12, 15)), type = "response", re.form = NA)

# #----------------------------------------
# ## Differences by trap style
# #----------------------------------------

df3 <- df %>% 
  mutate(trap_style = ifelse(grepl("N", buoy_number) == T, "rec", "science"))

review.mod <- glmmTMB::glmmTMB(lobster ~ bait_type*trap_style + (1|trap_pair), df3, family = nbinom2(link = "log"), ziformula = ~1)
pairs(emmeans::emmeans(review.mod, specs = "trap_style", by = "bait_type"), adjust = "bonferroni")

df3 <- df %>% 
  mutate(trap_style = ifelse(grepl("N", buoy_number) == T, "rec", "science")) %>%
  select(-total_crab, -green_crab) %>% 
  pivot_longer(cols = c(lobster:rock_crab))

mod1b <- glmmTMB(value ~ 0 + bait_type*name + deployment_nights + (1|trap_pair), data = df3, family = nbinom2(link = "log"), ziformula = ~1)

review.mod <- glmmTMB::glmmTMB(value ~ 0 + bait_type*name + trap_style + deployment_nights + (1|trap_pair), df3, family = nbinom2(link = "log"), ziformula = ~1)
summary(review.mod)

sjPlot::tab_model(mod1b, review.mod, show.aic = T, show.aicc = T, show.loglik = T, show.intercept = T, dv.labels = c("Main model", "Model w/ Trap style"), file = "Figures/mod_reviewer_forcomparison.html")

anova(mod1b, review.mod)

review.mod <- MASS::glm.nb(lobster ~ trap_style, df3)
summary(review.mod)

pairs(emmeans::emmeans(review.mod, specs = "trap_style"))


df3 <- df %>% 
  mutate(trap_style = ifelse(grepl("N", buoy_number) == T, "rec", "science")) %>% 
  filter(trap_style == "science")

# Full model
mod1 <- glmmTMB(lobster ~ bait_type + deployment_nights + (1|trap_pair), data = df3, family = nbinom2(link = "log"), ziformula = ~1)
summary(mod1)
simulated_residuals <- DHARMa::simulateResiduals(mod1)
plot(simulated_residuals)

sci.only <- data.frame(contrast(emmeans::emmeans(mod1, specs = "bait_type"), method = "revpairwise", type = "response", adjust = "bonferroni"))

df3 <- df %>% 
  mutate(trap_style = ifelse(grepl("N", buoy_number) == T, "rec", "science")) %>% 
  filter(trap_style == "rec")

# Full model
mod1 <- glmmTMB(lobster ~ bait_type + deployment_nights + (1|trap_pair), data = df3, family = nbinom2(link = "log"))
summary(mod1)
simulated_residuals <- DHARMa::simulateResiduals(mod1)
plot(simulated_residuals)

rec.only <- data.frame(contrast(emmeans::emmeans(mod1, specs = "bait_type"), method = "revpairwise", type = "response", adjust = "bonferroni"))

sci.only %>% 
  select(-c(df, null)) %>%
  mutate(data = "hake") %>% 
  bind_rows(rec.only %>% select(-c(df, null)) %>% mutate(data = "hoop")) %>%
  gt(groupname_col = "data")
  
  
#----------------------------
## Herring analysis
#----------------------------

me <- read.csv("Data/maine_landings.csv")

me %>% 
  ggplot(aes(x = as.integer(year), y = pounds))+
  geom_line()
me %>% 
  ggplot(aes(x = as.integer(year), y = price_per_pound))+
  geom_line()


scale_factor <- max(me$pounds) / max(me$price_per_pound)

me$millions_of_lbs <- me$pounds / 10^6
me$metric_tons <- me$pounds * 0.000453592
me$price_per_kg <- me$price_per_pound *2.2046
scale_factor <- max(me$metric_tons) / max(me$price_per_kg)

ggplot(me, aes(x = year)) +
  geom_line(aes(y = metric_tons), color = "blue", size = 1, linetype = 2) +
  geom_line(aes(y = price_per_kg * scale_factor), color = "red", size = 1) +
  scale_y_continuous(
    name = "Herring Catch (metric tons)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Price per kg (USD)")
  ) +
  theme_bd()

ggsave("Figures/herring.svg", width = 6, height = 4)













