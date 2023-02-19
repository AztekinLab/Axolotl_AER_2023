
## confusion matrix
library(cvms)
library(dplyr)
library(ggplot2)

confusion_matrix_cal <- function(targets, predictions) {
    targets <- as.character(targets)
    predictions <- as.character(predictions)
    res <- confusion_matrix(targets = targets, predictions = predictions)

    conf_mat <- res[["Confusion Matrix"]][[1]]
    conf_mat <- conf_mat %>%
        group_by(Target) %>%
        mutate(Percent = 100 * N / sum(N))

    # inconsistent cell types in tgts and pred will create NA
    conf_mat <- na.omit(conf_mat)
    return(conf_mat)
}
confusion_matrix_pl <- function(conf_mat, prefix,
                                title = "Confusion matrix", xlab = "Prediction", ylab = "Target",
                                save = TRUE, return = TRUE) {
    p <- ggplot(conf_mat, aes(Prediction, Target)) +
        geom_tile(aes(fill = Percent)) +
        geom_text(aes(label = round(Percent, 1)), size = 5) +
        ggtitle(label = title) +
        scale_fill_distiller(palette = "RdPu", direction = 1) +
        # theme_ipsum(base_family = "sans") +
        xlab(xlab) +
        ylab(ylab) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
            axis.title.x.bottom = element_text(hjust = 0.5, face = "bold", size = 15),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title.y.left = element_text(vjust = -0.5, face = "bold", size = 15)
        )
    if (save) {
        fn <- paste0(prefix, "_confusion_mat.pdf")
        ggsave(plot = p, filename = fn, width = 10, height = 9)
    }
    if (return) {
        return(p)
    }
}
