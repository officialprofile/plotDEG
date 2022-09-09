plotVolcano = function(data_results,
                       FCcolname = "log2FoldChange",
                       FCcutoff = 1,
                       FCbreaks = c(-100, -50, -20,-10,-4, -2, -1, 0, 1, 2, 4, 10, 20, 50, 100),
                       Pcolname = "padj",
                       Pcutoff = 0.05,
                       Pbreaks = c(1, 10, 25, 50, 100, 200, 500),
                       IDcolname = "ENSMBL",
                       title = " ",
                       pch = 21,
                       alpha = 1,
                       showDEGnumbers = F,
                       labels = c("Downregulated", "Not significant", "Upregulated"),
                       colors = c("#EE4000", "#EDEDED", "#4876FF"),
                       labelTop = 10,
                       labelColname = "SYMBOL",
                       ExprColname = "baseMean",
                       Exprcutoff = 10,
                       xlim = F)
{
  data_results = filter(data_results, !is.na(!!sym(FCcolname)), !is.na(!!sym(Pcolname)), !is.na(!!sym(IDcolname)))
  data_results[,".group"] = cut(data_results[[FCcolname]],
                                breaks = c(-Inf, -abs(FCcutoff), abs(FCcutoff), Inf),
                                labels = labels, right = FALSE)
  data_results[, ".radius"] = sqrt(data_results[,FCcolname]**2 + (-log(data_results[,Pcolname], 10))**2)*sign(data_results[,FCcolname])

  .labelsUP = filter(data_results, abs(!!sym(FCcolname)) > abs(FCcutoff), !!sym(Pcolname) < 0.05, !!sym(ExprColname) > Exprcutoff) %>%
    arrange(., -abs(!!sym(FCcolname))) %>% arrange(., .radius) %>% head(., labelTop)
  .labelsDOWN = filter(data_results, abs(!!sym(FCcolname)) > abs(FCcutoff), !!sym(Pcolname) < 0.05, !!sym(ExprColname) > Exprcutoff) %>%
    arrange(., -abs(!!sym(FCcolname))) %>%  arrange(., -.radius) %>% head(., labelTop)
  .labels = rbind(.labelsUP, .labelsDOWN) %>% filter(., !duplicated(!!sym(labelColname)))
  .maxfc = max(-min(data_results[,FCcolname]), max(data_results[,FCcolname]))
  .nUP = length(filter(data_results, !!sym(FCcolname) >= abs(FCcutoff))[,IDcolname])
  .nDOWN = length(filter(data_results, !!sym(FCcolname) <= -abs(FCcutoff))[,IDcolname])

  p <- ggplot(data_results, aes(x = !!sym(FCcolname), y = -log(!!sym(Pcolname), 10), col = .group)) +
    geom_point(pch = pch, cex = 1.5, alpha = alpha) +
    scale_color_manual(values = colors, labels = labels) +
    scale_x_continuous(breaks = FCbreaks, limits = c(-.maxfc, .maxfc)) +
    scale_y_continuous(breaks = Pbreaks, limits = c(0, max(-log(data_results[,Pcolname], 10)))) +
    theme_test() +
    theme(legend.position = "none", legend.title = element_blank(), text = element_text(size = 13)) +
    labs(x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR"), subtitle = title,
         caption = ifelse(showDEGnumbers, paste0('Upregulated genes: ', .nUP, ", downregulated: ", .nDOWN), " ")) +
    geom_vline(xintercept = c(-abs(FCcutoff), abs(FCcutoff)), col = "#444444", lty = 3) +
    geom_hline(yintercept = -log10(Pcutoff), col="#444444", lty = 3) +
    geom_label_repel(data = .labels, aes(label = !!sym(labelColname)),
                     size = 3, box.padding = unit(0.35, "lines"), point.padding = unit(0.1, "lines"))

  return(p)
}
