#!/usr/bin/env Rscript
# Generate PDF and PNG plots of stdpopsim demographic models
# Uses the rewritten PlotModel() function from PipeMaster

library(PipeMaster)
load("data_to_test/test_models.RData")

set.seed(123)

# ============================================================
# Combined PDF with all models
# ============================================================
pdf("data_to_test/stdpopsim_model_plots.pdf", width = 12, height = 15)
par(mfrow = c(3, 2))

PlotModel(Vaquita2Epoch, average.of.priors = TRUE, pop.labels = c("Vaquita"))
mtext("Vaquita2Epoch", side = 3, line = 0.5, col = "white", cex = 1.0, font = 2)

PlotModel(Africa_1T12, average.of.priors = TRUE, pop.labels = c("AFR"))
mtext("Africa_1T12", side = 3, line = 0.5, col = "white", cex = 1.0, font = 2)

PlotModel(PonAbe_TwoSpecies, average.of.priors = TRUE, pop.labels = c("Sumatran", "Bornean"))
mtext("PonAbe TwoSpecies", side = 3, line = 0.5, col = "white", cex = 1.0, font = 2)

PlotModel(OutOfAfrica_2T12, average.of.priors = TRUE, pop.labels = c("AFR", "EUR"))
mtext("OutOfAfrica_2T12", side = 3, line = 0.5, col = "white", cex = 1.0, font = 2)

PlotModel(OutOfAfrica_3G09, average.of.priors = TRUE, pop.labels = c("YRI", "CEU", "CHB"))
mtext("OutOfAfrica_3G09", side = 3, line = 0.5, col = "white", cex = 1.0, font = 2)

dev.off()
cat("Saved data_to_test/stdpopsim_model_plots.pdf\n")

# ============================================================
# Individual PNGs for tutorial embedding
# ============================================================
png_plot <- function(model, filename, pop.labels, title_text, width = 600, height = 500) {
  png(file.path("data_to_test", filename), width = width, height = height, res = 100)
  PlotModel(model, average.of.priors = TRUE, pop.labels = pop.labels)
  mtext(title_text, side = 3, line = 0.5, col = "white", cex = 1.0, font = 2)
  dev.off()
  cat(sprintf("  Saved data_to_test/%s\n", filename))
}

png_plot(Vaquita2Epoch, "model_Vaquita2Epoch.png",
         c("Vaquita"), "Vaquita2Epoch")

png_plot(Africa_1T12, "model_Africa_1T12.png",
         c("AFR"), "Africa_1T12")

png_plot(PonAbe_TwoSpecies, "model_PonAbe_TwoSpecies.png",
         c("Sumatran", "Bornean"), "PonAbe TwoSpecies")

png_plot(OutOfAfrica_3G09, "model_OutOfAfrica_3G09.png",
         c("YRI", "CEU", "CHB"), "OutOfAfrica_3G09")

cat("\nDone!\n")
