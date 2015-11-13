bdsky <-
  function(datadir = "~/Data", resdir = "~/Results", file = "Rparameters.csv") {
    setwd(datadir)
    
    packages <-
      c("multicore", "s20x", "boa", "Hmisc", "miscTools", "foreach", "doMC")
    lapply(packages, require, character.only = T)
    registerDoMC(2)
    
    # input : A matrix M and a column ascii name output : the numeric index of the column
    colnameindex = function(M, colname0) {
      colsnames = names(M[1,])
      theindex = which(colsnames == colname0)
      return(theindex)
    }
    
    # Parameters are read from a csv file in the working directory.
    # This allows for multiple files to be processed at once.
    
    parameters <- read.table(file, sep = ",", header = T)
    
    # Parameters per file are extracted from "Parameters.csv"

    loglist = as.character(parameters[, 1])
    burninpercent = parameters[, 4]
    recent = parameters[, 2]
    gridSize = 100
    origin = parameters[, 3]
    
    cat("Calculating.. ")
    
    foreach(i = 1:dim(parameters)[1]) %do% {
      # /* read and assign file from log list */
      assign(paste("log", i, sep = ""), read.table(loglist[i], header = T))
      attach(get(paste("log", i, sep = "")))
      
      R0_names = names(get(paste("log", i, sep = "")))[which(regexpr("R0\\w+", names(get(
        paste("log", i, sep = "")
      ))) > 0)]
      delta_names = names(get(paste("log", i, sep = "")))[which(regexpr("becomeUninfectiousRateS.\\w+", names(get(
        paste("log", i, sep = "")
      ))) > 0)]
      sampling_names = names(get(paste("log", i, sep = "")))[which(regexpr("samplingProportionS.\\w+", names(get(
        paste("log", i, sep = "")
      ))) > 0)]
      
      nsamples = length(get(R0_names[1]))
      treeheights = get(paste("log", i, sep = ""))[, match("treeheight", tolower(names(get(
        paste("log", i, sep = "")
      ))))]
      orig_roots = rep(origin[i], nsamples)
      width = median(treeheights + orig_roots)
      burnin = round(burninpercent[i] * nsamples / 100)
      
      intervalNumber = length(R0_names)
      if (intervalNumber > gridSize) {
        gridSize = intervalNumber
      }
      
      
      medians = matrix(data = NA, nrow = 1, ncol = gridSize)
      medians_G = matrix(data = NA, nrow = 1, ncol = gridSize)
      medians_H = matrix(data = NA, nrow = 1, ncol = gridSize)
      
      hpd_F = matrix(data = NA, nrow = 2, ncol = gridSize)
      hpd_G = matrix(data = NA, nrow = 2, ncol = gridSize)
      hpd_H = matrix(data = NA, nrow = 2, ncol = gridSize)
      
      F = matrix(data = NA, nrow = nsamples - burnin, ncol = gridSize)  #R0
      G = matrix(data = NA, nrow = nsamples - burnin, ncol = gridSize)  #becomeUninfectiousRateS.
      H = matrix(data = NA, nrow = nsamples - burnin, ncol = gridSize)  #samplingProportionS.
      
      step = width / (gridSize - 1)
      F_times = seq(recent[i] - width, recent[i], step)
      
      for (k in 1:(nsamples - burnin)) {
        time = treeheights[k] + orig_roots[k]
        
        for (l in 1:length(F_times)) {
          currentWidth = time / intervalNumber
          index = ceiling(intervalNumber - (recent[i] - F_times[l]) / currentWidth)
          
          F[k, l] = get(R0_names[max(index, 1)])[k + burnin]
          if (length(delta_names) == length(R0_names))
            G[k, l] = get(paste("becomeUninfectiousRateS.", max(index, 1), sep = ""))[k + burnin]
          else
            G[k, l] = becomeUninfectiousRateS.[k + burnin]
          if (length(sampling_names) == length(R0_names))
            H[k, l] = get(paste("samplingProportionS.", max(index, 1), sep = ""))[k + burnin]
          else
            H[k, l] = samplingProportionS.[k + burnin]
          
        }
      }
      
      foreach(j = 1:gridSize) %do% {
        if (length(which(F[, j] != "NA")) > (nsamples / 10)) {
          medians[1, j] = median(F[, j], na.rm = T)
          medians_G[1, j] = median(G[, j], na.rm = T)
          medians_H[1, j] = median(H[, j], na.rm = T)
          hpd_F[, j] = boa.hpd(F[which(F[, j] != "NA"), j], 0.05)[1:2]
          hpd_G[, j] = boa.hpd(G[which(G[, j] != "NA"), j], 0.05)[1:2]
          hpd_H[, j] = boa.hpd(H[which(H[, j] != "NA"), j], 0.05)[1:2]
        }
      }
      
      layout20x(3, 1)
      seed = strsplit(loglist[i], ".log")
      pdf(
        file = paste(resdir, seed, ".pdf", sep = ""), width = 13, height = 5
      )
      
      
      # /* plot R0 */
      plot(
        1, ylab = expression(R[0]), xlim = c(recent[i] - width, recent[i]), ylim = c(0, max(hpd_F[2,], na.rm = T) * 1.1), xlab = "year", col = "white", main = ""
      )
      minor.tick(nx = 5, ny = 2, tick.ratio = 0.2)
      polygon(c(F_times, rev(F_times)), c(hpd_F[2,], rev(hpd_F[1,])), col = "grey90", border = NA)
      lines(c(F_times), c(medians[1,]), type = "l")
      abline(1, 0, col = "grey")
      
      # /* plot become non-infectious rate */
      plot(
        1, ylab = expression(delta), xlim = c(recent[i] - width, recent[i]), ylim = c(0, max(hpd_G[2,], na.rm = T) * 1.1), xlab = "year", col = "white", main = ""
      )
      minor.tick(nx = 5, ny = 2, tick.ratio = 0.2)
      polygon(c(F_times, rev(F_times)), c(hpd_G[2,], rev(hpd_G[1,])), col = "grey90", border = NA)
      lines(F_times, medians_G[1,], type = "l")
      
      # /* plot samplingProportionS. */
      plot(
        0, ylab = expression(s), xlim = c(recent[i] - width, recent[i]), ylim = c(0, max(hpd_H[2,], na.rm = T) * 1.1), xlab = "year", col = "white", main = ""
      )
      polygon(c(F_times, rev(F_times)), c(hpd_H[2,], rev(hpd_H[1,])), col = "grey90", border = NA)
      lines(F_times, medians_H[1,], type = "l")
      
      
      cat("Plotting finished.")
      
      dev.off()
      
      plot.data <-
        data.frame(
          F_times, UR = hpd_F[2,], MR = medians[1,], LR = hpd_F[1,], Ld = hpd_G[1,], Md = medians_G[1,], Ud = hpd_G[2,], Ls = hpd_H[1,], Ms = medians_H[1,], Us = hpd_H[2,]
        )
      
      write.table(
        plot.data, file = paste0(resdir, seed, ".csv"), row.names = FALSE, sep = "\t"
      )
    }
    
  }