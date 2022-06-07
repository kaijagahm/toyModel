# Model testing and results
source("agentBased.R")

### Test-run the simulation and plot the results
set.seed(333)
test <- runSim(tmax = 10, n = 15, allowIsolated = FALSE, verbose = FALSE)
plotSim(test$gs) 
# keeping the same coordinates throughout to make it clearer how the edges 
# are changing

### Another test with a longer time frame: 50 days
set.seed(333)
test500 <- runSim(tmax = 500, n = 15, allowIsolated = FALSE, verbose = FALSE)

### Examine temporal dynamics (500 days)
stats <- getNodeStats(test500$gs, type = "df")
head(stats)
stats %>%
  ggplot(aes(x = Day, y = jitter(deg)))+
  geom_point(alpha = 0.05)+
  geom_smooth(size = 2)+
  theme_minimal() +
  ylab("Degree")
# Interesting, the oscillations seem to dampen as time goes on. 
# That makes sense. I think I've created solid baseline dynamics!

### Allow isolated nodes--does behavior change? (500 days)
set.seed(333)
test500 <- runSim(tmax = 500, n = 15, allowIsolated = TRUE, verbose = FALSE)
stats <- getNodeStats(test500$gs, type = "df")
head(stats)
stats %>%
  ggplot(aes(x = Day, y = jitter(deg)))+
  geom_point(alpha = 0.05)+
  geom_smooth(size = 2)+
  theme_minimal()+
  ylab("Degree")
# Approximately the same dynamics. Maybe slightly more fluctuations, but maybe not.
