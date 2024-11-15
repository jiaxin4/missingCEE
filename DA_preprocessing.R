
rm(list = ls())
library(tidyverse)


suggest <- as.tibble(read.csv("suggestions.csv"))
suggest <- suggest %>% 
  filter(!is.na(decision.index.nogap))

### assign a new decision index variable bc there were gaps in the old one from travelling
suggest_tdr <- suggest %>% 
  group_by(user.index) %>% 
  mutate(max.decision.index.nogap = n(),
         decision.index.nogap.new = 0:(unique(max.decision.index.nogap) - 1)) %>% 
  ungroup()

PRINT_SANITY_CHECK <- TRUE

if (PRINT_SANITY_CHECK) {
  # make sure primary keys are unique
  suggest_tdr %>% 
    count(user.index, decision.index.nogap.new) %>% 
    filter(n > 1) # no duplicate!
}

if (PRINT_SANITY_CHECK) {
  # make sure all user has decision.index.nogap from 0 to max consecutively, no gap in between
  print(suggest_tdr %>% group_by(user.index) %>%
          summarise(n_dp = length(decision.index.nogap.new),
                    dp_min = min(decision.index.nogap.new),
                    dp_max = max(decision.index.nogap.new),
                    no_gap_in_dp = ((dp_max - dp_min + 1) == n_dp)), n = Inf)
}
# # A tibble: 37 × 5
# user.index  n_dp dp_min dp_max no_gap_in_dp
# <int> <int>  <int>  <int> <lgl>       
#   1          1   178      0    177 TRUE        
# 2          2   209      0    208 TRUE        
# 3          3   215      0    214 TRUE        
# 4          4   219      0    218 TRUE        
# 5          5   217      0    216 TRUE        
# 6          6   182      0    181 TRUE        
# 7          7   216      0    215 TRUE        
# 8          8   221      0    220 TRUE        
# 9          9   207      0    206 TRUE        
# 10         10   215      0    214 TRUE        
# 11         11   220      0    219 TRUE        
# 12         12   244      0    243 TRUE        
# 13         13   190      0    189 TRUE        
# 14         14   210      0    209 TRUE        
# 15         15   256      0    255 TRUE        
# 16         16   199      0    198 TRUE        
# 17         17   210      0    209 TRUE        
# 18         18   211      0    210 TRUE        
# 19         19   210      0    209 TRUE        
# 20         20   235      0    234 TRUE        
# 21         21   235      0    234 TRUE        
# 22         22   214      0    213 TRUE        
# 23         23   225      0    224 TRUE        
# 24         24   175      0    174 TRUE        
# 25         25   223      0    222 TRUE        
# 26         26   212      0    211 TRUE        
# 27         27   211      0    210 TRUE        
# 28         28   210      0    209 TRUE        
# 29         29   211      0    210 TRUE        
# 30         30   210      0    209 TRUE        
# 31         31   184      0    183 TRUE        
# 32         32   225      0    224 TRUE        
# 33         33   228      0    227 TRUE        
# 34         34   212      0    211 TRUE        
# 35         35   213      0    212 TRUE        
# 36         36   202      0    201 TRUE        
# 37         37   230      0    229 TRUE        
# > 
dat <- suggest_tdr %>% 
  transmute(user = user.index,
            decision.index = decision.index.nogap,
            avail = case_when(avail == "True" ~ 1,
                               avail == "False"~ 0),
            send = case_when(send == "True" ~ 1,
                             send == "False"|send == "" ~ 0),
            send.active = case_when(send.active == "True" ~ 1,
                                    send.active == "False"|send.active == "" ~ 0),
            send.sedentary = case_when(send.sedentary == "True" ~ 1,
                                       send.sedentary == "False"|send.sedentary == "" ~ 0),
            steps30 = jbsteps30,
            steps60 = jbsteps60,
            steps30pre.zero = jbsteps30pre.zero,
            steps30pre.zero.bin_56 = ifelse(steps30pre.zero >= 56, 1, 0),
            steps30pre.log = log(steps30pre.zero + 0.5),
            dec.location.category = dec.location.category,
            homework.location = ifelse(dec.location.category %in% c("home", "work"), 1, 0),
            weekday =  weekdays(ymd_hms(sugg.decision.utime)),
            is.weekday = ifelse(weekday %in% c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"), 1, 0),
            steps30_observe_indicator = ifelse(is.na(steps30), 0, 1))

dat$logsteps30 <- log(dat$steps30 + 0.5)
dat$logsteps60 <- log(dat$steps60 + 0.5)
#documentation for missing: Jawbone steps, but NA was imputed with 0. NAs occurred in a user's Jawbone data e.g. if the user just didn't wear the Jawbone for a whole day, or if she was traveling internationally and left the tracker at home.

saveRDS(dat, "suggestion_cleaned.RDS")

