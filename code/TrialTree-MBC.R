# =============================================================================
# TrialTree-MBC: Median PFS & OS Prediction Tools
# Advanced Breast Cancer Phase II Trials — Decision Tree Models
# PFS Model: n = 868; 24 terminal nodes; MAE = 1.92; RMSE = 2.51
# OS Model:  n = 860; 21 terminal nodes; MAE = 4.54; RMSE = 5.86
# =============================================================================

library(shiny)
library(bslib)

# =============================================================================
# HELPER: return list for approximate (internal node) predictions
# =============================================================================
make_approx <- function(yval, n_pct, path, node_desc) {
  path <- c(path, "\u26a0 Variable not available [Stopped at internal node]")
  list(prediction = yval, n = NA, pct = n_pct, path = path,
       node_id = "internal", node_desc = node_desc, approximate = TRUE)
}

# =============================================================================
# PFS PREDICTION (24 terminal nodes + internal node fallbacks)
# =============================================================================
predict_mpfs <- function(orr, wecog, line_of_therapy, subtype, treatment_size, platinum) {
  path <- c()
  N <- 868
  
  # ROOT: ORR < 40.0
  if (is.na(orr)) return(make_approx(6.70, "100%", path, "Root node (all trials)"))
  
  if (orr < 40.0) {
    path <- c(path, "ORR < 40.0 [Yes]")
    if (orr < 23.7) {
      path <- c(path, "ORR < 23.7 [Yes]")
      if (orr < 13.0) {
        path <- c(path, "ORR < 13.0 [Yes]")
        if (is.na(subtype)) return(make_approx(3.04, "12.3%", path, "Low ORR (<13.0%), subtype unknown"))
        if (subtype %in% c("HER2+", "Mixed", "TNBC")) {
          path <- c(path, paste0("Subtype = ", subtype, " [Yes: HER2/Mixed/TNBC]"))
          if (orr < 6.8) {
            path <- c(path, "ORR < 6.8 [Yes]")
            n <- 28
            return(list(prediction = 2.05, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 32, node_desc = "Low ORR (<6.8%), HER2+/Mixed/TNBC", approximate = FALSE))
          } else {
            path <- c(path, "ORR < 6.8 [No]")
            n <- 24
            return(list(prediction = 3.00, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 33, node_desc = "Low ORR (6.8\u201313.0%), HER2+/Mixed/TNBC", approximate = FALSE))
          }
        } else {
          path <- c(path, paste0("Subtype = ", subtype, " [No: HR+]"))
          if (is.na(wecog)) return(make_approx(3.56, "6.3%", path, "Low ORR (<13.0%), HR+, wECOG unknown"))
          if (wecog >= 0.46) {
            path <- c(path, "wECOG >= 0.46 [Yes]")
            n <- 28
            return(list(prediction = 3.13, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 34, node_desc = paste0("Low ORR (<13.0%), HR+, ", classify_wecog(wecog)$label), approximate = FALSE))
          } else {
            path <- c(path, "wECOG >= 0.46 [No]")
            n <- 27
            return(list(prediction = 4.01, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 35, node_desc = paste0("Low ORR (<13.0%), HR+, ", classify_wecog(wecog)$label), approximate = FALSE))
          }
        }
      } else {
        path <- c(path, "ORR < 13.0 [No]")
        if (is.na(wecog)) return(make_approx(4.51, "16.8%", path, "Low ORR (13.0\u201323.7%), wECOG unknown"))
        if (wecog >= 0.52) {
          path <- c(path, "wECOG >= 0.52 [Yes]")
          if (orr < 16.7) {
            path <- c(path, "ORR < 16.7 [Yes]")
            n <- 35
            return(list(prediction = 3.35, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 36, node_desc = paste0("Low ORR (13.0\u201316.7%), ", classify_wecog(wecog)$label), approximate = FALSE))
          } else {
            path <- c(path, "ORR < 16.7 [No]")
            n <- 48
            return(list(prediction = 4.25, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 37, node_desc = paste0("Low ORR (16.7\u201323.7%), ", classify_wecog(wecog)$label), approximate = FALSE))
          }
        } else {
          path <- c(path, "wECOG >= 0.52 [No]")
          if (is.na(treatment_size)) return(make_approx(5.35, "7.3%", path, "Low ORR, better PS, treatment unknown"))
          if (treatment_size == "1-Agent") {
            path <- c(path, "Treatment Size = 1-Agent [Yes]")
            n <- 39
            return(list(prediction = 4.81, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 38, node_desc = paste0("Low ORR (13.0\u201323.7%), ", classify_wecog(wecog)$label, ", monotherapy"), approximate = FALSE))
          } else {
            path <- c(path, paste0("Treatment Size = ", treatment_size, " [No: combination]"))
            n <- 24
            return(list(prediction = 6.24, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 39, node_desc = paste0("Low ORR (13.0\u201323.7%), ", classify_wecog(wecog)$label, ", combination therapy"), approximate = FALSE))
          }
        }
      }
    } else {
      path <- c(path, "ORR < 23.7 [No]")
      if (is.na(line_of_therapy)) return(make_approx(6.24, "25.5%", path, "Low\u2013Moderate ORR (23.7\u201340.0%), line unknown"))
      if (line_of_therapy == "Pretreated") {
        path <- c(path, "Line of Therapy = Pretreated [Yes]")
        if (is.na(treatment_size)) return(make_approx(5.92, "21.0%", path, "Low\u2013Moderate ORR, pretreated, treatment unknown"))
        if (treatment_size == "1-Agent") {
          path <- c(path, "Treatment Size = 1-Agent [Yes]")
          if (is.na(wecog)) return(make_approx(5.17, "9.2%", path, "Low\u2013Moderate ORR, pretreated, mono, wECOG unknown"))
          if (wecog >= 0.37) {
            path <- c(path, "wECOG >= 0.37 [Yes]")
            n <- 66
            return(list(prediction = 4.90, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 40, node_desc = paste0("Low\u2013Moderate ORR (23.7\u201340.0%), pretreated, monotherapy, ", classify_wecog(wecog)$label), approximate = FALSE))
          } else {
            path <- c(path, "wECOG >= 0.37 [No]")
            n <- 14
            return(list(prediction = 6.45, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 41, node_desc = paste0("Low\u2013Moderate ORR (23.7\u201340.0%), pretreated, monotherapy, ", classify_wecog(wecog)$label), approximate = FALSE))
          }
        } else {
          path <- c(path, paste0("Treatment Size = ", treatment_size, " [No: combination]"))
          if (is.na(platinum)) return(make_approx(6.51, "11.8%", path, "Low\u2013Moderate ORR, pretreated, combo, platinum unknown"))
          if (platinum == 1) {
            path <- c(path, "Platinum-containing = Yes")
            n <- 15
            return(list(prediction = 5.11, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 42, node_desc = "Low\u2013Moderate ORR, pretreated, platinum combo", approximate = FALSE))
          } else {
            path <- c(path, "Platinum-containing = No")
            n <- 87
            return(list(prediction = 6.75, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 43, node_desc = "Low\u2013Moderate ORR, pretreated, non-platinum combo", approximate = FALSE))
          }
        }
      } else {
        path <- c(path, "Line of Therapy = Pretreated [No: First-line]")
        n <- 39
        return(list(prediction = 7.76, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 11, node_desc = "Low\u2013Moderate ORR (23.7\u201340.0%), first-line", approximate = FALSE))
      }
    }
  } else {
    path <- c(path, "ORR < 40.0 [No]")
    if (orr < 59.1) {
      path <- c(path, "ORR < 59.1 [Yes]")
      if (is.na(wecog)) return(make_approx(7.81, "24.1%", path, "Moderate\u2013High ORR (40\u201359%), wECOG unknown"))
      if (wecog >= 0.75) {
        path <- c(path, "wECOG >= 0.75 [Yes]")
        if (is.na(subtype)) return(make_approx(6.60, "6.0%", path, "Moderate\u2013High ORR, Worse PS, subtype unknown"))
        if (subtype %in% c("HER2+", "Mixed", "TNBC")) {
          path <- c(path, paste0("Subtype = ", subtype, " [Yes: HER2/Mixed/TNBC]"))
          n <- 24
          return(list(prediction = 5.69, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 24, node_desc = paste0("Moderate\u2013High ORR (40\u201359%), ", classify_wecog(wecog)$label, ", HER2+/Mixed/TNBC"), approximate = FALSE))
        } else {
          path <- c(path, paste0("Subtype = ", subtype, " [No: HR+]"))
          n <- 28
          return(list(prediction = 7.37, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 25, node_desc = paste0("Moderate\u2013High ORR (40\u201359%), ", classify_wecog(wecog)$label, ", HR+"), approximate = FALSE))
        }
      } else {
        path <- c(path, "wECOG >= 0.75 [No]")
        if (is.na(treatment_size)) return(make_approx(8.21, "18.1%", path, "Moderate\u2013High ORR, better PS, treatment unknown"))
        if (treatment_size %in% c("1-Agent", "2-Agent")) {
          path <- c(path, paste0("Treatment Size = ", treatment_size, " [Yes: 1-/2-Agent]"))
          if (is.na(subtype)) return(make_approx(8.02, "15.7%", path, "Moderate\u2013High ORR, mono/doublet, subtype unknown"))
          if (subtype == "Mixed") {
            path <- c(path, "Subtype = Mixed [Yes]")
            n <- 45
            return(list(prediction = 7.18, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 52, node_desc = paste0("Moderate\u2013High ORR (40\u201359%), ", classify_wecog(wecog)$label, ", mono/doublet, Mixed"), approximate = FALSE))
          } else {
            path <- c(path, paste0("Subtype = ", subtype, " [No: not Mixed]"))
            n <- 91
            return(list(prediction = 8.43, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 53, node_desc = paste0("Moderate\u2013High ORR (40\u201359%), ", classify_wecog(wecog)$label, ", mono/doublet, non-Mixed"), approximate = FALSE))
          }
        } else {
          path <- c(path, paste0("Treatment Size = ", treatment_size, " [No: 3-Agent]"))
          n <- 21
          return(list(prediction = 9.47, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 27, node_desc = paste0("Moderate\u2013High ORR (40\u201359%), ", classify_wecog(wecog)$label, ", triplet"), approximate = FALSE))
        }
      }
    } else {
      path <- c(path, "ORR < 59.1 [No]")
      if (is.na(wecog)) return(make_approx(10.23, "19.0%", path, "High ORR (\u226559%), wECOG unknown"))
      if (wecog >= 0.43) {
        path <- c(path, "wECOG >= 0.43 [Yes]")
        if (wecog < 0.81) {
          path <- c(path, "wECOG < 0.81 [Yes]")
          if (wecog >= 0.64) {
            path <- c(path, "wECOG >= 0.64 [Yes]")
            n <- 23
            return(list(prediction = 7.98, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 56, node_desc = paste0("High ORR (\u226559%), ", classify_wecog(wecog)$label), approximate = FALSE))
          } else {
            path <- c(path, "wECOG >= 0.64 [No]")
            n <- 55
            return(list(prediction = 9.63, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 57, node_desc = paste0("High ORR (\u226559%), ", classify_wecog(wecog)$label), approximate = FALSE))
          }
        } else {
          path <- c(path, "wECOG < 0.81 [No]")
          n <- 23
          return(list(prediction = 10.77, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 29, node_desc = paste0("High ORR (\u226559%), ", classify_wecog(wecog)$label), approximate = FALSE))
        }
      } else {
        path <- c(path, "wECOG >= 0.43 [No]")
        if (orr < 74.5) {
          path <- c(path, "ORR < 74.5 [Yes]")
          if (wecog >= 0.35) {
            path <- c(path, "wECOG >= 0.35 [Yes]")
            n <- 15
            return(list(prediction = 9.67, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 60, node_desc = paste0("High ORR (59\u201374.5%), ", classify_wecog(wecog)$label), approximate = FALSE))
          } else {
            path <- c(path, "wECOG >= 0.35 [No]")
            n <- 30
            return(list(prediction = 10.99, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 61, node_desc = paste0("High ORR (59\u201374.5%), ", classify_wecog(wecog)$label), approximate = FALSE))
          }
        } else {
          path <- c(path, "ORR < 74.5 [No]")
          n <- 19
          return(list(prediction = 13.31, n = n, pct = paste0(round(100*n/N,1),"%"), path = path, node_id = 31, node_desc = paste0("High ORR (\u226574.5%), ", classify_wecog(wecog)$label), approximate = FALSE))
        }
      }
    }
  }
}

# =============================================================================
# OS PREDICTION (21 terminal nodes + internal node fallbacks)
# =============================================================================
predict_mos <- function(median_pfs, orr, wecog, subtype, line_of_therapy,
                        aromatase_inhibitor, her2_kinase_inhibitor,
                        antiestrogen, therapy) {
  path <- c()
  N <- 860
  
  # ROOT: Median_PFS < 7.2
  if (is.na(median_pfs)) return(make_approx(19, "100%", path, "Root node (all trials)"))
  
  if (median_pfs < 7.2) {
    path <- c(path, "Median PFS < 7.2 [Yes]")
    if (median_pfs < 4.2) {
      path <- c(path, "Median PFS < 4.2 [Yes]")
      if (is.na(aromatase_inhibitor)) return(make_approx(12, "27%", path, paste0(classify_mpfs(median_pfs)$label, " PFS (<4.2 mo), AI status unknown")))
      if (aromatase_inhibitor == 0) {
        path <- c(path, "Aromatase Inhibitor = 0 [Yes]")
        if (is.na(subtype)) return(make_approx(12, "26%", path, paste0(classify_mpfs(median_pfs)$label, " PFS (<4.2 mo), no AI, subtype unknown")))
        if (subtype %in% c("HR+", "Mixed", "TNBC")) {
          path <- c(path, paste0("Subtype = ", subtype, " [Yes: HR/Mixed/TNBC]"))
          if (median_pfs < 1.9) {
            path <- c(path, "Median PFS < 1.9 [Yes]")
            n <- round(0.03 * N)
            return(list(prediction = 8.7, n = n, pct = "3%", path = path, node_id = 1, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (<1.9 mo), no AI, HR/Mixed/TNBC"), approximate = FALSE))
          } else {
            path <- c(path, "Median PFS < 1.9 [No]")
            if (is.na(wecog)) return(make_approx(12, "20%", path, paste0(classify_mpfs(median_pfs)$label, " PFS (1.9\u20134.2 mo), no AI, HR/Mixed/TNBC, wECOG unknown")))
            if (wecog >= 0.62) {
              path <- c(path, "wECOG >= 0.62 [Yes]")
              n <- round(0.10 * N)
              return(list(prediction = 11, n = n, pct = "10%", path = path, node_id = 2, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (1.9\u20134.2 mo), no AI, HR/Mixed/TNBC, ", classify_wecog(wecog)$label), approximate = FALSE))
            } else {
              path <- c(path, "wECOG >= 0.62 [No]")
              n <- round(0.10 * N)
              return(list(prediction = 13, n = n, pct = "10%", path = path, node_id = 3, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (1.9\u20134.2 mo), no AI, HR/Mixed/TNBC, ", classify_wecog(wecog)$label), approximate = FALSE))
            }
          }
        } else {
          path <- c(path, paste0("Subtype = ", subtype, " [No: not HR/Mixed/TNBC]"))
          n <- round(0.03 * N)
          return(list(prediction = 15, n = n, pct = "3%", path = path, node_id = 4, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (<4.2 mo), no AI, HER2+"), approximate = FALSE))
        }
      } else {
        path <- c(path, "Aromatase Inhibitor = 0 [No: AI present]")
        n <- round(0.02 * N)
        return(list(prediction = 21, n = n, pct = "2%", path = path, node_id = 5, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (<4.2 mo), AI present"), approximate = FALSE))
      }
    } else {
      path <- c(path, "Median PFS < 4.2 [No]")
      if (is.na(wecog)) return(make_approx(17, "34%", path, paste0(classify_mpfs(median_pfs)$label, " PFS (4.2\u20137.2 mo), wECOG unknown")))
      if (wecog >= 0.71) {
        path <- c(path, "wECOG >= 0.71 [Yes]")
        n <- round(0.10 * N)
        return(list(prediction = 13, n = n, pct = "10%", path = path, node_id = 6, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (4.2\u20137.2 mo), ", classify_wecog(wecog)$label), approximate = FALSE))
      } else {
        path <- c(path, "wECOG >= 0.71 [No]")
        if (is.na(her2_kinase_inhibitor)) return(make_approx(19, "24%", path, paste0(classify_mpfs(median_pfs)$label, " PFS, better PS, HER2 KI status unknown")))
        if (her2_kinase_inhibitor == 0) {
          path <- c(path, "HER2 Kinase Inhibitor = 0 [Yes]")
          if (is.na(antiestrogen)) return(make_approx(18, "23%", path, paste0(classify_mpfs(median_pfs)$label, " PFS, no HER2 KI, antiestrogen status unknown")))
          if (antiestrogen == 0) {
            path <- c(path, "Antiestrogen = 0 [Yes]")
            if (is.na(line_of_therapy)) return(make_approx(18, "21%", path, paste0(classify_mpfs(median_pfs)$label, " PFS, no HER2 KI/antiestrogen, line unknown")))
            if (line_of_therapy == "Pretreated") {
              path <- c(path, "Line of Therapy = Pretreated [Yes]")
              if (is.na(therapy)) return(make_approx(17, "16%", path, paste0(classify_mpfs(median_pfs)$label, " PFS, pretreated, therapy class unknown")))
              if (therapy == "Chemo") {
                path <- c(path, "Therapy = Chemo [Yes]")
                if (is.na(orr)) return(make_approx(16, "9%", path, paste0(classify_mpfs(median_pfs)$label, " PFS, pretreated, chemo, ORR unknown")))
                if (orr >= 24) {
                  path <- c(path, "ORR >= 24 [Yes]")
                  n <- round(0.07 * N)
                  return(list(prediction = 15, n = n, pct = "7%", path = path, node_id = 7, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS, pretreated, chemo, Low\u2013Moderate ORR (\u226524%)"), approximate = FALSE))
                } else {
                  path <- c(path, "ORR >= 24 [No]")
                  n <- round(0.02 * N)
                  return(list(prediction = 21, n = n, pct = "2%", path = path, node_id = 8, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS, pretreated, chemo, Low ORR (<24%)"), approximate = FALSE))
                }
              } else {
                path <- c(path, "Therapy = Chemo [No]")
                n <- round(0.08 * N)
                return(list(prediction = 19, n = n, pct = "8%", path = path, node_id = 9, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS, pretreated, non-chemo therapy"), approximate = FALSE))
              }
            } else {
              path <- c(path, "Line of Therapy = Pretreated [No: First-line]")
              n <- round(0.05 * N)
              return(list(prediction = 20, n = n, pct = "5%", path = path, node_id = 10, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS, first-line, no anti-HER2/antiestrogen"), approximate = FALSE))
            }
          } else {
            path <- c(path, "Antiestrogen = 0 [No: antiestrogen present]")
            n <- round(0.02 * N)
            return(list(prediction = 24, n = n, pct = "2%", path = path, node_id = 11, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS, antiestrogen present"), approximate = FALSE))
          }
        } else {
          path <- c(path, "HER2 Kinase Inhibitor = 0 [No: HER2 KI present]")
          n <- round(0.02 * N)
          return(list(prediction = 26, n = n, pct = "2%", path = path, node_id = 12, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS, HER2 kinase inhibitor present"), approximate = FALSE))
        }
      }
    }
  } else {
    path <- c(path, "Median PFS < 7.2 [No]")
    if (median_pfs < 11) {
      path <- c(path, "Median PFS < 11 [Yes]")
      if (is.na(wecog)) return(make_approx(22, "26%", path, paste0(classify_mpfs(median_pfs)$label, " PFS (7.2\u201311 mo), wECOG unknown")))
      if (wecog >= 0.75) {
        path <- c(path, "wECOG >= 0.75 [Yes]")
        n <- round(0.05 * N)
        return(list(prediction = 18, n = n, pct = "5%", path = path, node_id = 13, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (7.2\u201311 mo), ", classify_wecog(wecog)$label), approximate = FALSE))
      } else {
        path <- c(path, "wECOG >= 0.75 [No]")
        if (median_pfs < 9.2) {
          path <- c(path, "Median PFS < 9.2 [Yes]")
          if (median_pfs >= 8.3) {
            path <- c(path, "Median PFS >= 8.3 [Yes]")
            n <- round(0.06 * N)
            return(list(prediction = 21, n = n, pct = "6%", path = path, node_id = 14, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (8.3\u20139.2 mo), ", classify_wecog(wecog)$label), approximate = FALSE))
          } else {
            path <- c(path, "Median PFS >= 8.3 [No]")
            n <- round(0.08 * N)
            return(list(prediction = 23, n = n, pct = "8%", path = path, node_id = 15, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (7.2\u20138.3 mo), ", classify_wecog(wecog)$label), approximate = FALSE))
          }
        } else {
          path <- c(path, "Median PFS < 9.2 [No]")
          if (is.na(orr)) return(make_approx(26, "8%", path, paste0(classify_mpfs(median_pfs)$label, " PFS (9.2\u201311 mo), ORR unknown")))
          if (orr < 48) {
            path <- c(path, "ORR < 48 [Yes]")
            n <- round(0.02 * N)
            return(list(prediction = 22, n = n, pct = "2%", path = path, node_id = 16, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (9.2\u201311 mo), Low\u2013Moderate ORR (<48%)"), approximate = FALSE))
          } else {
            path <- c(path, "ORR < 48 [No]")
            n <- round(0.06 * N)
            return(list(prediction = 27, n = n, pct = "6%", path = path, node_id = 17, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (9.2\u201311 mo), Moderate\u2013High ORR (\u226548%)"), approximate = FALSE))
          }
        }
      }
    } else {
      path <- c(path, "Median PFS < 11 [No]")
      if (is.na(therapy)) return(make_approx(29, "12%", path, paste0(classify_mpfs(median_pfs)$label, " PFS (\u226511 mo), therapy class unknown")))
      if (therapy == "Chemo") {
        path <- c(path, "Therapy = Chemo [Yes]")
        if (is.na(wecog)) return(make_approx(27, "6%", path, paste0(classify_mpfs(median_pfs)$label, " PFS, chemo, wECOG unknown")))
        if (wecog >= 0.43) {
          path <- c(path, "wECOG >= 0.43 [Yes]")
          n <- round(0.04 * N)
          return(list(prediction = 24, n = n, pct = "4%", path = path, node_id = 18, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (\u226511 mo), chemo, ", classify_wecog(wecog)$label), approximate = FALSE))
        } else {
          path <- c(path, "wECOG >= 0.43 [No]")
          n <- round(0.03 * N)
          return(list(prediction = 30, n = n, pct = "3%", path = path, node_id = 19, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (\u226511 mo), chemo, ", classify_wecog(wecog)$label), approximate = FALSE))
        }
      } else {
        path <- c(path, "Therapy = Chemo [No]")
        if (is.na(wecog)) return(make_approx(32, "6%", path, paste0(classify_mpfs(median_pfs)$label, " PFS, non-chemo, wECOG unknown")))
        if (wecog < 0.38) {
          path <- c(path, "wECOG < 0.38 [Yes]")
          n <- round(0.03 * N)
          return(list(prediction = 29, n = n, pct = "3%", path = path, node_id = 20, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (\u226511 mo), non-chemo, ", classify_wecog(wecog)$label), approximate = FALSE))
        } else {
          path <- c(path, "wECOG < 0.38 [No]")
          n <- round(0.03 * N)
          return(list(prediction = 34, n = n, pct = "3%", path = path, node_id = 21, node_desc = paste0(classify_mpfs(median_pfs)$label, " PFS (\u226511 mo), non-chemo, ", classify_wecog(wecog)$label), approximate = FALSE))
        }
      }
    }
  }
}


# =============================================================================
# SHARED UI HELPERS
# =============================================================================

render_decision_path <- function(path_steps, pred, endpoint_label, pred_val) {
  step_items <- lapply(seq_along(path_steps), function(i) {
    is_last <- (i == length(path_steps))
    step_text <- path_steps[i]
    if (grepl("\u26a0", step_text)) { badge_color <- "#8B0000"
    } else if (grepl("\\[Yes", step_text)) { badge_color <- "#27ae60"
    } else if (grepl("\\[No", step_text)) { badge_color <- "#e74c3c"
    } else { badge_color <- "#2980b9" }
    condition <- gsub(" \\[.*\\]", "", step_text)
    direction <- if (grepl("\u26a0", step_text)) "Stopped" else gsub(".*\\[(.*)\\]", "\\1", step_text)
    tags$div(
      style = paste0(
        "display: flex; align-items: center; gap: 12px; padding: 8px 12px; ",
        "border-left: 3px solid ", if (is_last) "#2c3e50" else "#dee2e6", "; ",
        "margin-left: ", (i - 1) * 16, "px; ",
        if (is_last) "background: #f8f9fa; border-radius: 0 4px 4px 0; font-weight: 600;" else ""
      ),
      tags$span(paste0("Step ", i), style = "font-size: 0.7rem; color: #aaa; min-width: 42px;"),
      tags$span(condition, style = "font-size: 0.85rem; color: #2c3e50;"),
      tags$span(direction, style = paste0(
        "font-size: 0.72rem; padding: 2px 8px; border-radius: 10px; ",
        "background: ", badge_color, "; color: white; font-weight: 600;"))
    )
  })
  
  node_label <- if (isTRUE(pred$approximate)) "Internal Node" else paste0("Node ", pred$node_id)
  node_bg <- if (isTRUE(pred$approximate)) "#8B4513" else "#2c3e50"
  
  step_items <- c(step_items, list(
    tags$div(
      style = paste0(
        "display: flex; align-items: center; gap: 12px; padding: 10px 14px; ",
        "margin-left: ", length(path_steps) * 16, "px; ",
        "background: ", node_bg, "; border-radius: 6px; margin-top: 6px;"),
      tags$span(icon(if (isTRUE(pred$approximate)) "exclamation-triangle" else "leaf"),
                style = paste0("color: ", if (isTRUE(pred$approximate)) "#FFA500" else "#27ae60", "; font-size: 1rem;")),
      tags$span(
        paste0(node_label, ": ", endpoint_label, " \u2248 ",
               sprintf("%.2f", pred_val), " months",
               if (isTRUE(pred$approximate)) " (approximate)" else ""),
        style = "font-size: 0.95rem; color: white; font-weight: 700;")
    )
  ))
  do.call(tags$div, step_items)
}

get_pfs_category <- function(val) {
  if (val < 4.1) list(color="#e74c3c", bg="#fdf2f2", border="#e74c3c", label="Likely Poor",
                      interp="Below the 25th percentile (Q1 = 4.1 mo) of the full dataset.")
  else if (val <= 6.3) list(color="#e67e22", bg="#fef9f0", border="#e67e22", label="Moderate",
                            interp="Between Q1 (4.1 mo) and the median (6.3 mo) of the full dataset.")
  else if (val <= 9.1) list(color="#2980b9", bg="#f0f7fc", border="#2980b9", label="Likely Good",
                            interp="Between the median (6.3 mo) and Q3 (9.1 mo) of the full dataset.")
  else list(color="#27ae60", bg="#f0faf4", border="#27ae60", label="Likely Favorable",
            interp="Above the 75th percentile (Q3 = 9.1 mo) of the full dataset.")
}

get_os_category <- function(val) {
  if (val < 12) list(color="#e74c3c", bg="#fdf2f2", border="#e74c3c", label="Likely Poor",
                     interp="Below the 25th percentile (Q1 = 12 mo) of the full dataset.")
  else if (val <= 17) list(color="#e67e22", bg="#fef9f0", border="#e67e22", label="Moderate",
                           interp="Between Q1 (12 mo) and the median (17 mo) of the full dataset.")
  else if (val <= 24) list(color="#2980b9", bg="#f0f7fc", border="#2980b9", label="Likely Good",
                           interp="Between the median (17 mo) and Q3 (24 mo) of the full dataset.")
  else list(color="#27ae60", bg="#f0faf4", border="#27ae60", label="Likely Favorable",
            interp="Above the 75th percentile (Q3 = 24 mo) of the full dataset.")
}

# --- Input classifiers based on dataset quartiles ---
classify_orr <- function(val) {
  if (is.na(val)) return(list(label = "N/A", color = "#999"))
  if (val < 21) list(label = "Low", color = "#e74c3c")
  else if (val <= 37) list(label = "Low\u2013Moderate", color = "#e67e22")
  else if (val <= 55) list(label = "Moderate\u2013High", color = "#2980b9")
  else list(label = "High", color = "#27ae60")
}

classify_wecog <- function(val) {
  if (is.na(val)) return(list(label = "N/A", color = "#999"))
  if (val < 0.38) list(label = "Good PS", color = "#27ae60")
  else if (val <= 0.52) list(label = "Moderate PS", color = "#2980b9")
  else if (val <= 0.71) list(label = "Poor PS", color = "#e67e22")
  else list(label = "Worse PS", color = "#e74c3c")
}

classify_mpfs <- function(val) {
  if (is.na(val)) return(list(label = "N/A", color = "#999"))
  if (val < 4.1) list(label = "Poor", color = "#e74c3c")
  else if (val <= 6.3) list(label = "Moderate", color = "#e67e22")
  else if (val <= 9.1) list(label = "Good", color = "#2980b9")
  else list(label = "Favorable", color = "#27ae60")
}

render_result_card <- function(pred, cat, endpoint_label, mae_str, rmse_str, dist_str, input_summary = NULL) {
  approx_banner <- if (isTRUE(pred$approximate)) {
    tags$div(
      style = "background: #FFF3CD; border: 1px solid #FFC107; border-radius: 6px; padding: 10px; margin-bottom: 12px; font-size: 0.85rem;",
      tags$strong("\u26a0 Approximate Prediction: "),
      "One or more input variables were not available. The prediction was derived from an ",
      tags$strong("internal (non-terminal) node"), " of the decision tree, representing the ",
      "weighted average of all terminal nodes below that split. This estimate is less precise ",
      "than a terminal node prediction."
    )
  }
  
  n_info <- if (isTRUE(pred$approximate)) {
    paste0("Internal node | ", pred$pct, " of training data")
  } else {
    paste0("Node ", pred$node_id, " | n = ", pred$n, " (", pred$pct, " of training data)")
  }
  
  card(
    style = paste0("border-left: 5px solid ", cat$border, "; background: ", cat$bg, ";"),
    card_body(
      approx_banner,
      tags$div(
        style = "display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap;",
        tags$div(
          tags$span(paste0("Predicted ", endpoint_label), style = "font-size: 0.9rem; color: #555; display: block;"),
          tags$span(paste0(sprintf("%.2f", pred$prediction), " months"),
                    style = paste0("font-size: 2.8rem; font-weight: 700; color: ", cat$color, "; line-height: 1.1;")),
          tags$span(paste0("  (", cat$label, ")"),
                    style = paste0("font-size: 1rem; color: ", cat$color, "; font-weight: 500;"))
        ),
        tags$div(
          style = "text-align: right; min-width: 220px;",
          tags$div(style = "font-size: 0.8rem; color: #777;", n_info),
          tags$div(style = "font-size: 0.8rem; color: #777; margin-top: 4px;",
                   paste0("Model error: MAE = ", mae_str, ", RMSE = ", rmse_str)),
          tags$div(style = paste0("margin-top: 8px; font-size: 0.82rem; color: ", cat$color, "; font-weight: 500;"),
                   pred$node_desc)
        )
      ),
      tags$hr(style = "margin: 10px 0;"),
      tags$div(style = "font-size: 1.1rem; color: #555;",
               tags$strong("Interpretation: "), cat$interp, tags$br(),
               tags$span(style = "color: #888; font-size: 0.75rem;", dist_str)),
      if (!is.null(input_summary)) tagList(
        tags$hr(style = "margin: 10px 0;"),
        tags$div(style = "font-size: 0.85rem; color: #444;",
                 tags$strong("Input Summary: "), input_summary,
                 tags$br(),
                 tags$span(style = "color: #999; font-size: 0.72rem; font-style: italic;",
                           "ORR and wECOG classifications derived from dataset quartiles (ORR: Q1=21%, Median=37%, Q3=55%; wECOG: Q1=0.38, Median=0.52, Q3=0.71)."))
      )
    )
  )
}


# =============================================================================
# VISIT COUNTER (file-based, persists between sessions)
# =============================================================================
counter_file <- "visit_count.rds"
if (!file.exists(counter_file)) saveRDS(0L, counter_file)

increment_counter <- function() {
  count <- tryCatch(readRDS(counter_file), error = function(e) 0L)
  count <- count + 1L
  saveRDS(count, counter_file)
  count
}

get_counter <- function() {
  tryCatch(readRDS(counter_file), error = function(e) 0L)
}

# =============================================================================
# UI
# =============================================================================

app_theme <- bs_theme(version = 5, bootswatch = "flatly", primary = "#2c3e50",
                      success = "#27ae60", info = "#2980b9")

na_checkbox_style <- "margin-top: -8px; margin-bottom: 8px; font-size: 0.82rem; color: #888;"

ui <- navbarPage(
  id = "main_nav",
  title = tags$span(HTML('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1350 420" style="height:60px;width:auto;vertical-align:middle;margin-top:-2px;">
    <text x="84" y="224" font-family="Georgia, Times New Roman, serif" font-size="130" font-weight="700" letter-spacing="-0.5" fill="#ffffff">Trial<tspan fill="#a8d5c8">Tree</tspan></text>
    <text x="717" y="224" font-family="Georgia, Times New Roman, serif" font-size="130" font-weight="400" fill="#ffffff">-</text>
    <text x="770" y="224" font-family="Georgia, Times New Roman, serif" font-size="122" font-weight="700" fill="#e8b4c8">MBC</text>
    <line x1="72" y1="270" x2="1085" y2="270" stroke="rgba(255,255,255,0.4)" stroke-width="5"/>
    <text x="72" y="345" font-family="Arial, Helvetica, sans-serif" font-size="44" font-weight="700" fill="rgba(255,255,255,0.6)">Predict Survival Time | Metastatic Breast Cancer</text>
  </svg>')),
  theme = app_theme,
  windowTitle = "TrialTree-MBC | Trial-Level Survival Prediction",
  header = tags$head(
    tags$meta(property = "og:title", content = "TrialTree-MBC"),
    tags$meta(property = "og:description", content = "Trial-Level Survival Prediction for Metastatic Breast Cancer. Interactive decision tree models for median PFS and OS in phase II trials."),
    tags$meta(property = "og:type", content = "website"),
    tags$meta(property = "og:url", content = "https://f-azam.shinyapps.io/TrialTree-MBC/"),
    tags$meta(name = "description", content = "TrialTree-MBC: Interactive decision tree-based prediction tool for median PFS and OS in advanced breast cancer phase II trials."),
    tags$style(HTML("
    .model-card { cursor: pointer; transition: all 0.2s ease; border: 2px solid #dee2e6; border-radius: 12px; padding: 30px; text-align: center; }
    .model-card:hover { border-color: #2c3e50; box-shadow: 0 4px 15px rgba(0,0,0,0.1); transform: translateY(-2px); }
    .navbar { border-bottom: 2px solid #2c3e50; padding-top: 0; padding-bottom: 0; }
    .navbar .navbar-brand { padding-top: 2px; padding-bottom: 2px; }
    .navbar .nav-link { font-size: 1rem; }
  "))),
  
  # ---- Landing Page ----
  tabPanel("Home",
           div(style = "max-width: 900px; margin: 10px auto; padding: 0 20px;",
               # --- Logo and header ---
               tags$div(style = "text-align: center; padding: 10px 20px 25px; border-bottom: 2px solid #2c3e50; position: relative;",
                        # --- Visit counter (top right) ---
                        uiOutput("visit_counter"),
                        tags$h2(style = "font-size: 2.2rem; font-weight: 700; color: #2c3e50; margin: 0; letter-spacing: -0.5px;",
                                "TrialTree-MBC ",
                                tags$span("v1.3", style = "font-size: 1rem; font-weight: 500; color: #888;")),
                        tags$p("Trial-Level Survival Prediction for Metastatic Breast Cancer",
                               style = "font-size: 1.05rem; color: #555; margin-top: 10px; font-weight: 500;"),
                        tags$p(style = "font-size: 0.88rem; color: #999; margin-top: -5px;",
                               "Decision tree-based prediction models for median PFS and OS in advanced breast cancer phase II trials.",
                               tags$br(), "Select a model to begin.")
               ),
               div(style = "display: grid; grid-template-columns: 1fr 1fr; gap: 30px; margin-top: 25px;",
                   div(class = "model-card", onclick = "Shiny.setInputValue('goto_pfs', Math.random())",
                       tags$div(style = "font-size: 2.5rem; margin-bottom: 10px;", icon("chart-line")),
                       tags$h4("Median PFS Model", style = "font-weight: 600; color: #2c3e50;"),
                       tags$p(style = "color: #777; font-size: 0.9rem;", "Predict median progression-free survival"),
                       tags$hr(style = "margin: 15px 40px;"),
                       tags$div(style = "font-size: 0.8rem; color: #888;",
                                "6 predictors | 24 terminal nodes", tags$br(),
                                "Training: n = 868 | MAE = 1.92 mo", tags$br(),
                                "Predictors: ORR, wECOG, subtype, line of therapy, treatment size, platinum")),
                   div(class = "model-card", onclick = "Shiny.setInputValue('goto_os', Math.random())",
                       tags$div(style = "font-size: 2.5rem; margin-bottom: 10px;", icon("chart-line")),
                       tags$h4("Median OS Model", style = "font-weight: 600; color: #2c3e50;"),
                       tags$p(style = "color: #777; font-size: 0.9rem;", "Predict median overall survival"),
                       tags$hr(style = "margin: 15px 40px;"),
                       tags$div(style = "font-size: 0.8rem; color: #888;",
                                "9 predictors | 21 terminal nodes", tags$br(),
                                "Training: n = 860 | MAE = 4.54 mo", tags$br(),
                                "Predictors: mPFS, ORR, wECOG, subtype, therapy class, drug indicators"))
               ),
               tags$div(style = "text-align: center; margin-top: 15px; padding: 15px; background: #f8f9fa; border-radius: 8px; font-size: 0.75rem; color: #8B0000;",
                        tags$strong("Disclaimer: "),
                        "TrialTree-MBC offers population-level benchmarks using aggregate data from phase II breast cancer trials. ",
                        "It predicts expected median PFS and OS months based on patient and trial characteristics. ",
                        "These predictions may serve as contextual references for clinical trial design and endpoint interpretation. ",
                        "This tool is intended for research purposes only and is not a substitute for professional medical advice. ",
                        "Models are trained on published phase II advanced or metastatic breast cancer trials (1992\u20132022) and may not apply to phase III trials, real-world settings, or other cancers."),
               tags$p(style = "text-align: center; margin-top: 2px; font-size: 0.85rem; color: #666;",
                      "For any inquiries, please contact Faruque Azam at ",
                      tags$a(href = "mailto:faruque.azam@bracu.ac.bd", "faruque.azam@bracu.ac.bd",
                             style = "color: #2980b9; text-decoration: none; font-weight: 500;")),
               # --- Logo at bottom ---
               tags$div(style = "text-align: center; margin-top: -15px;",
                        tags$div(style = "display: inline-block; max-width: 420px; margin-left: 65px;",
                                 HTML('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1350 420" style="width:75%;height:auto;">
  <text x="84" y="224" font-family="Georgia, Times New Roman, serif" font-size="130" font-weight="700" letter-spacing="-0.5" fill="#102030">Trial<tspan fill="#7ebba8">Tree</tspan></text>
  <text x="717" y="224" font-family="Georgia, Times New Roman, serif" font-size="130" font-weight="400" fill="#102030">-</text>
  <text x="770" y="224" font-family="Georgia, Times New Roman, serif" font-size="122" font-weight="700" fill="#ba7893">MBC</text>
  <line x1="72" y1="270" x2="1085" y2="270" stroke="#d7dbe1" stroke-width="5"/>
  <text x="72" y="345" font-family="Arial, Helvetica, sans-serif" font-size="44" font-weight="700" fill="#7e8594">Predict Survival Time | Metastatic Breast Cancer</text>
</svg>')))
           )
  ),
  
  # ---- PFS Model Tab ----
  tabPanel("PFS Model",
           layout_sidebar(
             sidebar = sidebar(width = 370,
                               tags$h5(icon("chart-line"), " PFS Model Inputs", style = "font-weight: 600; color: #2c3e50; margin-bottom: 15px;"),
                               tags$hr(style = "margin: 5px 0 15px;"),
                               tags$label("Overall Response Rate (ORR, %)", style = "font-weight: 600; font-size: 0.95rem;"),
                               sliderInput("pfs_orr", NULL, min = 0, max = 100, value = 35, step = 0.5, post = "%"),
                               tags$label("Weighted ECOG Score (wECOG)", style = "font-weight: 600; font-size: 0.95rem;"),
                               tags$p("Trial-weighted average ECOG performance status", style = "font-size: 0.75rem; color: #888; margin-bottom: 6px;"),
                               checkboxInput("pfs_wecog_na", "Not Available", FALSE),
                               conditionalPanel("!input.pfs_wecog_na", sliderInput("pfs_wecog", NULL, min = 0, max = 1.5, value = 0.55, step = 0.01)),
                               tags$hr(style = "margin: 10px 0;"),
                               selectInput("pfs_subtype", "Breast Cancer Subtype", c("HR+"="HR+","HER2+"="HER2+","TNBC"="TNBC","Mixed"="Mixed","Not Available"="NA"), "HR+"),
                               selectInput("pfs_lot", "Line of Therapy", c("Pretreated"="Pretreated","First-line (1L)"="First-line","Not Available"="NA"), "Pretreated"),
                               selectInput("pfs_tx", "Treatment Regimen", c("1-Agent (Mono)"="1-Agent","2-Agent (Doublet)"="2-Agent","3-Agent (Triplet)"="3-Agent","Not Available"="NA"), "1-Agent"),
                               selectInput("pfs_plat", "Platinum-Containing", c("No"=0,"Yes"=1,"Not Available"="NA"), 0),
                               tags$hr(style = "margin: 15px 0 10px;"),
                               actionButton("pfs_predict_btn", "Predict Median PFS", class = "btn-primary btn-lg w-100", icon = icon("calculator")),
                               tags$div(style = "margin-top: -10px; text-align: center;",
                                        tags$div(style = "display: inline-block; max-width: 240px; margin-left: 20px;",
                                                 HTML('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1350 420" style="width:60%;height:auto;">
  <text x="84" y="224" font-family="Georgia, Times New Roman, serif" font-size="130" font-weight="700" letter-spacing="-0.5" fill="#102030">Trial<tspan fill="#7ebba8">Tree</tspan></text>
  <text x="717" y="224" font-family="Georgia, Times New Roman, serif" font-size="130" font-weight="400" fill="#102030">-</text>
  <text x="770" y="224" font-family="Georgia, Times New Roman, serif" font-size="122" font-weight="700" fill="#ba7893">MBC</text>
  <line x1="72" y1="270" x2="1085" y2="270" stroke="#d7dbe1" stroke-width="5"/>
  <text x="72" y="345" font-family="Arial, Helvetica, sans-serif" font-size="44" font-weight="700" fill="#7e8594">Predict Survival Time | Metastatic Breast Cancer</text>
</svg>')))
             ),
             uiOutput("pfs_result"),
             card(card_header(tags$h5(icon("code-branch"), " Decision Path", style = "margin:0; font-weight:600;")),
                  card_body(uiOutput("pfs_path"))),
             tags$p(style = "font-size: 0.85rem; color: #555; margin: 10px 0;",
                    icon("sitemap"), " To visualize the full PFS tree structure, ",
                    tags$a(href = "https://1drv.ms/i/c/1ede7f09fdbafa1e/IQDx_r9xX3t7RYBkTcaCK3dUAWm63U5Zs33ZJU92GPOGxNw?e=K1bmOP",
                           "click here", target = "_blank", style = "color: #2980b9; font-weight: 500;")),
             card(card_header(tags$h5(icon("info-circle"), " PFS Model Information", style = "margin:0; font-weight:600;")),
                  card_body(
                    tags$div(style = "display: grid; grid-template-columns: 1fr 1fr; gap: 20px;",
                             tags$div(tags$h6("Specifications", style = "font-weight:600;"),
                                      tags$table(class = "table table-sm table-borderless", style = "font-size: 0.85rem;",
                                                 tags$tr(tags$td("Algorithm:"), tags$td(tags$strong("CART Regression Tree (rpart)"))),
                                                 tags$tr(tags$td("Training:"), tags$td(tags$strong("n = 868 phase II trials"))),
                                                 tags$tr(tags$td("Max tree depth:"), tags$td(tags$strong("5"))),
                                                 tags$tr(tags$td("Cost-complexity (cp):"), tags$td(tags$strong("1e-4"))),
                                                 tags$tr(tags$td("Min node size:"), tags$td(tags$strong("40"))),
                                                 tags$tr(tags$td("Internal nodes:"), tags$td(tags$strong("23"))),
                                                 tags$tr(tags$td("Terminal nodes:"), tags$td(tags$strong("24"))))),
                             tags$div(tags$h6("Performance", style = "font-weight:600;"),
                                      tags$table(class = "table table-sm table-borderless", style = "font-size: 0.85rem;",
                                                 tags$tr(tags$td("MAE:"), tags$td(tags$strong("1.92 mo"))),
                                                 tags$tr(tags$td("RMSE:"), tags$td(tags$strong("2.51 mo"))),
                                                 tags$tr(tags$td("Population:"), tags$td(tags$strong("Advanced or metastatic breast cancer"))),
                                                 tags$tr(tags$td("Outcome:"), tags$td(tags$strong("Median PFS (mo)"))),
                                                 tags$tr(tags$td("Primary split:"), tags$td(tags$strong("ORR (<40.0%)"))))))
                  ))
           )
  ),
  
  # ---- OS Model Tab ----
  tabPanel("OS Model",
           layout_sidebar(
             sidebar = sidebar(width = 370,
                               tags$h5(icon("chart-line"), " OS Model Inputs", style = "font-weight: 600; color: #2c3e50; margin-bottom: 15px;"),
                               tags$hr(style = "margin: 5px 0 15px;"),
                               tags$label("Median PFS (months)", style = "font-weight: 600; font-size: 0.95rem;"),
                               tags$p("Reported or predicted median PFS from the trial", style = "font-size: 0.75rem; color: #888; margin-bottom: 6px;"),
                               numericInput("os_mpfs", NULL, value = 6.0, min = 0, max = 40, step = 0.1),
                               tags$label("Overall Response Rate (ORR, %)", style = "font-weight: 600; font-size: 0.95rem;"),
                               checkboxInput("os_orr_na", "Not Available", FALSE),
                               conditionalPanel("!input.os_orr_na", sliderInput("os_orr", NULL, min = 0, max = 100, value = 35, step = 0.5, post = "%")),
                               tags$label("Weighted ECOG Score (wECOG)", style = "font-weight: 600; font-size: 0.95rem;"),
                               tags$p("Trial-weighted average ECOG performance status", style = "font-size: 0.75rem; color: #888; margin-bottom: 6px;"),
                               checkboxInput("os_wecog_na", "Not Available", FALSE),
                               conditionalPanel("!input.os_wecog_na", sliderInput("os_wecog", NULL, min = 0, max = 1.5, value = 0.55, step = 0.01)),
                               tags$hr(style = "margin: 10px 0;"),
                               selectInput("os_subtype", "Breast Cancer Subtype", c("HR+"="HR+","HER2+"="HER2+","TNBC"="TNBC","Mixed"="Mixed","Not Available"="NA"), "HR+"),
                               selectInput("os_lot", "Line of Therapy", c("Pretreated"="Pretreated","First-line (1L)"="First-line","Not Available"="NA"), "Pretreated"),
                               selectInput("os_therapy", "Therapy Class", c("Chemotherapy"="Chemo","Non-Chemotherapy"="Non-Chemo","Not Available"="NA"), "Chemo"),
                               tags$hr(style = "margin: 10px 0;"),
                               tags$label("Drug Class Indicators", style = "font-weight: 600; font-size: 0.95rem;"),
                               selectInput("os_ai", "Aromatase Inhibitor", c("No"=0,"Yes"=1,"Not Available"="NA"), 0),
                               selectInput("os_h2ki", "HER2 Kinase Inhibitor", c("No"=0,"Yes"=1,"Not Available"="NA"), 0),
                               selectInput("os_anti", "Antiestrogen", c("No"=0,"Yes"=1,"Not Available"="NA"), 0),
                               tags$hr(style = "margin: 15px 0 10px;"),
                               actionButton("os_predict_btn", "Predict Median OS", class = "btn-primary btn-lg w-100", icon = icon("calculator")),
                               tags$div(style = "margin-top: -10px; text-align: center;",
                                        tags$div(style = "display: inline-block; max-width: 240px; margin-left: 20px;",
                                                 HTML('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1350 420" style="width:60%;height:auto;">
  <text x="84" y="224" font-family="Georgia, Times New Roman, serif" font-size="130" font-weight="700" letter-spacing="-0.5" fill="#102030">Trial<tspan fill="#7ebba8">Tree</tspan></text>
  <text x="717" y="224" font-family="Georgia, Times New Roman, serif" font-size="130" font-weight="400" fill="#102030">-</text>
  <text x="770" y="224" font-family="Georgia, Times New Roman, serif" font-size="122" font-weight="700" fill="#ba7893">MBC</text>
  <line x1="72" y1="270" x2="1085" y2="270" stroke="#d7dbe1" stroke-width="5"/>
  <text x="72" y="345" font-family="Arial, Helvetica, sans-serif" font-size="44" font-weight="700" fill="#7e8594">Predict Survival Time | Metastatic Breast Cancer</text>
</svg>')))
             ),
             uiOutput("os_result"),
             card(card_header(tags$h5(icon("code-branch"), " Decision Path", style = "margin:0; font-weight:600;")),
                  card_body(uiOutput("os_path"))),
             tags$p(style = "font-size: 0.85rem; color: #555; margin: 10px 0;",
                    icon("sitemap"), " To visualize the full OS tree structure, ",
                    tags$a(href = "https://1drv.ms/i/c/1ede7f09fdbafa1e/IQC2GW_bA0eXRpmQ16wuiLALAVMG1vF61Z46EUEcf1laPOI?e=8ffNC2",
                           "click here", target = "_blank", style = "color: #2980b9; font-weight: 500;")),
             card(card_header(tags$h5(icon("info-circle"), " OS Model Information", style = "margin:0; font-weight:600;")),
                  card_body(
                    tags$div(style = "display: grid; grid-template-columns: 1fr 1fr; gap: 20px;",
                             tags$div(tags$h6("Specifications", style = "font-weight:600;"),
                                      tags$table(class = "table table-sm table-borderless", style = "font-size: 0.85rem;",
                                                 tags$tr(tags$td("Algorithm:"), tags$td(tags$strong("CART Regression Tree (rpart)"))),
                                                 tags$tr(tags$td("Training:"), tags$td(tags$strong("n = 860 phase II trials"))),
                                                 tags$tr(tags$td("Max tree depth:"), tags$td(tags$strong("5"))),
                                                 tags$tr(tags$td("Cost-complexity (cp):"), tags$td(tags$strong("3.16e-3"))),
                                                 tags$tr(tags$td("Min node size:"), tags$td(tags$strong("40"))),
                                                 tags$tr(tags$td("Internal nodes:"), tags$td(tags$strong("20"))),
                                                 tags$tr(tags$td("Terminal nodes:"), tags$td(tags$strong("21"))))),
                             tags$div(tags$h6("Performance", style = "font-weight:600;"),
                                      tags$table(class = "table table-sm table-borderless", style = "font-size: 0.85rem;",
                                                 tags$tr(tags$td("MAE:"), tags$td(tags$strong("4.54 mo"))),
                                                 tags$tr(tags$td("RMSE:"), tags$td(tags$strong("5.86 mo"))),
                                                 tags$tr(tags$td("Population:"), tags$td(tags$strong("Advanced or metastatic breast cancer"))),
                                                 tags$tr(tags$td("Outcome:"), tags$td(tags$strong("Median OS (mo)"))),
                                                 tags$tr(tags$td("Primary split:"), tags$td(tags$strong("mPFS (<7.2 mo)"))))))
                  ))
           )
  )
)


# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {
  
  # --- Visit counter: increment once per session ---
  visit_count <- increment_counter()
  
  output$visit_counter <- renderUI({
    tags$div(style = "position: absolute; top: 12px; right: 0; text-align: right;",
             tags$div(style = "display: inline-flex; align-items: center; gap: 6px; background: #f0f4f8; border: 1px solid #dee2e6; border-radius: 20px; padding: 5px 14px;",
                      icon("eye", style = "color: #2980b9; font-size: 0.8rem;"),
                      tags$span(formatC(visit_count, format = "d", big.mark = ","),
                                style = "font-size: 0.85rem; font-weight: 600; color: #2c3e50;"),
                      tags$span("visits", style = "font-size: 0.75rem; color: #888;")
             )
    )
  })
  
  observeEvent(input$goto_pfs, { updateNavbarPage(session, "main_nav", selected = "PFS Model") })
  observeEvent(input$goto_os, { updateNavbarPage(session, "main_nav", selected = "OS Model") })
  
  # --- PFS Model ---
  pfs_pred <- eventReactive(input$pfs_predict_btn, {
    orr_val     <- input$pfs_orr
    wecog_val   <- if (input$pfs_wecog_na) NA_real_ else input$pfs_wecog
    sub_val     <- if (input$pfs_subtype == "NA") NA_character_ else input$pfs_subtype
    lot_val     <- if (input$pfs_lot == "NA") NA_character_ else input$pfs_lot
    tx_val      <- if (input$pfs_tx == "NA") NA_character_ else input$pfs_tx
    plat_val    <- if (input$pfs_plat == "NA") NA_integer_ else as.integer(input$pfs_plat)
    predict_mpfs(orr = orr_val, wecog = wecog_val, line_of_therapy = lot_val,
                 subtype = sub_val, treatment_size = tx_val, platinum = plat_val)
  })
  
  output$pfs_result <- renderUI({
    req(pfs_pred())
    pred <- pfs_pred()
    cat <- get_pfs_category(pred$prediction)
    orr_val <- input$pfs_orr
    wecog_val <- if (input$pfs_wecog_na) NA_real_ else input$pfs_wecog
    orr_cls <- classify_orr(orr_val)
    wecog_cls <- classify_wecog(wecog_val)
    input_summary <- tagList(
      tags$span(paste0("ORR = ", if (is.na(orr_val)) "N/A" else paste0(orr_val, "%"), " "),
                tags$span(paste0("(", orr_cls$label, ")"),
                          style = paste0("color:", orr_cls$color, "; font-weight:600;"))),
      tags$span(" | "),
      tags$span(paste0("wECOG = ", if (is.na(wecog_val)) "N/A" else sprintf("%.2f", wecog_val), " "),
                tags$span(paste0("(", wecog_cls$label, ")"),
                          style = paste0("color:", wecog_cls$color, "; font-weight:600;")))
    )
    render_result_card(pred, cat, "Median PFS", "1.92 mo", "2.51 mo",
                       "Categories derived from dataset distribution: median PFS = 6.3 months (IQR: 4.1\u20139.1). Poor: < Q1; Moderate: Q1\u2013median; Good: median\u2013Q3; Favorable: > Q3.",
                       input_summary = input_summary)
  })
  
  output$pfs_path <- renderUI({
    req(pfs_pred())
    pred <- pfs_pred()
    render_decision_path(pred$path, pred, "mPFS", pred$prediction)
  })
  
  # --- OS Model ---
  os_pred <- eventReactive(input$os_predict_btn, {
    mpfs_val  <- input$os_mpfs
    orr_val   <- if (input$os_orr_na) NA_real_ else input$os_orr
    wecog_val <- if (input$os_wecog_na) NA_real_ else input$os_wecog
    sub_val   <- if (input$os_subtype == "NA") NA_character_ else input$os_subtype
    lot_val   <- if (input$os_lot == "NA") NA_character_ else input$os_lot
    ther_val  <- if (input$os_therapy == "NA") NA_character_ else input$os_therapy
    ai_val    <- if (input$os_ai == "NA") NA_integer_ else as.integer(input$os_ai)
    h2ki_val  <- if (input$os_h2ki == "NA") NA_integer_ else as.integer(input$os_h2ki)
    anti_val  <- if (input$os_anti == "NA") NA_integer_ else as.integer(input$os_anti)
    predict_mos(median_pfs = mpfs_val, orr = orr_val, wecog = wecog_val,
                subtype = sub_val, line_of_therapy = lot_val,
                aromatase_inhibitor = ai_val, her2_kinase_inhibitor = h2ki_val,
                antiestrogen = anti_val, therapy = ther_val)
  })
  
  output$os_result <- renderUI({
    req(os_pred())
    pred <- os_pred()
    cat <- get_os_category(pred$prediction)
    mpfs_val <- input$os_mpfs
    orr_val <- if (input$os_orr_na) NA_real_ else input$os_orr
    wecog_val <- if (input$os_wecog_na) NA_real_ else input$os_wecog
    mpfs_cls <- classify_mpfs(mpfs_val)
    orr_cls <- classify_orr(orr_val)
    wecog_cls <- classify_wecog(wecog_val)
    input_summary <- tagList(
      tags$span(paste0("mPFS = ", if (is.na(mpfs_val)) "N/A" else paste0(mpfs_val, " mo"), " "),
                tags$span(paste0("(", mpfs_cls$label, ")"),
                          style = paste0("color:", mpfs_cls$color, "; font-weight:600;"))),
      tags$span(" | "),
      tags$span(paste0("ORR = ", if (is.na(orr_val)) "N/A" else paste0(orr_val, "%"), " "),
                tags$span(paste0("(", orr_cls$label, ")"),
                          style = paste0("color:", orr_cls$color, "; font-weight:600;"))),
      tags$span(" | "),
      tags$span(paste0("wECOG = ", if (is.na(wecog_val)) "N/A" else sprintf("%.2f", wecog_val), " "),
                tags$span(paste0("(", wecog_cls$label, ")"),
                          style = paste0("color:", wecog_cls$color, "; font-weight:600;")))
    )
    render_result_card(pred, cat, "Median OS", "4.54 mo", "5.86 mo",
                       "Categories derived from dataset distribution: median OS = 17.0 months (IQR: 12\u201324). Poor: < Q1; Moderate: Q1\u2013median; Good: median\u2013Q3; Favorable: > Q3.",
                       input_summary = input_summary)
  })
  
  output$os_path <- renderUI({
    req(os_pred())
    pred <- os_pred()
    render_decision_path(pred$path, pred, "mOS", pred$prediction)
  })
}

shinyApp(ui = ui, server = server)
