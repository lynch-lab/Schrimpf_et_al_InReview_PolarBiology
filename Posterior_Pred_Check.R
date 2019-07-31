# Code to examine the posterior predictive check results

# The model code includes the steps for calculating the T_actual and T_predicted test statistics.

# Load those results for a single species:
load(file="./Model_Output/GEPE.out.ppc2.Rdata")


# Bayesian P-value:
sum((t.pred - t.actual) > 0) / length(t.actual)

# Function to make the figure:
PostPredFig <- function(t.actual, t.pred, main) {
  plot(t.pred[which(t.actual > quantile(t.actual, 0.005) &
                      t.pred > quantile(t.pred, 0.005) &
                      t.actual < quantile(t.actual, 0.995) &
                      t.pred < quantile(t.pred, 0.995))],
       t.actual[which(t.actual > quantile(t.actual, 0.005) &
                        t.pred > quantile(t.pred, 0.005) &
                        t.actual < quantile(t.actual, 0.995) &
                        t.pred < quantile(t.pred, 0.995))],
       xlab = expression('T'['Predicted']),
       ylab = expression('T'['Actual']),
       pch = 16, cex = 0.5, col = rgb(0,0,0,0.3),
       main = main)
  mtext("Posterior Predictive Check")
  abline(a=0, b=1)
}

# Draw the figure:
PostPredFig(t.actual = t.actual, t.pred = t.pred, main = "GEPE")