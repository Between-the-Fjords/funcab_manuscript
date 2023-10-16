# Set the model output folder
modelOutputFolder <- "C:/Users/joseph.chipperfield/OneDrive - NINA/Work/FUNCAB/ImputationModel"
# Set the file import locations
imputeDataLocation <- "C:/Users/joseph.chipperfield/OneDrive - NINA/Work/FUNCAB/ImputationModel/CompiledTraitData_B.csv"
responseDataLocation <- "C:/Users/joseph.chipperfield/OneDrive - NINA/Work/FUNCAB/ImputationModel/CompiledResponseData_D.csv"

# Import the data from the csv files
traitData <- read.csv(imputeDataLocation)
plotData <- read.csv(responseDataLocation)
# Remove any empty lines in the csv files
traitData <- traitData[apply(X = as.matrix(traitData), FUN = function(curRow) {!all(is.na(curRow))}, MARGIN = 1), ]
plotData <- plotData[apply(X = as.matrix(plotData), FUN = function(curRow) {!all(is.na(curRow))}, MARGIN = 1), ]

# Extra functional group information not present in the trait and plot data files
extraFunGroupInfo <- as.data.frame(matrix(c(
	"Emp_nig", "woody",
	"Emp_her", "woody",
	"Jun_sp", "woody",
	"Hyp_mac", "forb",
	"Phl_sp", "graminoid",
#	"Pop_tre", "tree",
	"Ste_als", "forb",
	"Eri_bor", "forb",
	"Gen_sp", "forb",
	"Sax_aiz_", "forb",
	"Tri_arv", "forb",
	"Pedicularis", "forb",
	"Arenaria", "forb",
	"Pilosella", "forb",
	"Ste_sp", "forb",
	"Luz_tri", "graminoid",
	"Cre_pal", "forb",
#	"Sch_gig", "graminoid",
	"Porub", "forb"
), ncol = 2, dimnames = list(NULL, c("Species", "FunGroup")), byrow = TRUE))

# Important headers in the trait data
traitData_speciesCol <- "Species"					  		# The name of the species column
traitData_funGroup <- "functionalGroup"					# The name of the functional group column
traitData_siteCol <- "Site"							  			# The name of the site column
traitData_traitCols <- c(                       # The names of the columns denoting traits
	"LDMC", "Lth_ave", "logHeight", "logLA", "SLA", "N", "C", "CN_ratio"
)

# Set the number of traits measured in the analysis
numTraits <- length(traitData_traitCols)

# Important headers in the plot data
plotData_speciesCol <- "species"								# The name of the species column
plotData_funGroup <- "functionalGroup"					# The name of the functional group column
plotData_siteCol <- "Site"											# The name of the site column
plotData_plotCol <- "turfID_Y"									# The name of the plot column
plotData_coverCol <- "cover"										# The name of the cover column

# Replace any species bionomials with '.' seperators with '_' instead
traitData[, traitData_speciesCol] <- as.factor(gsub(".", "_", as.character(traitData[, traitData_speciesCol]), fixed = TRUE))
plotData[, plotData_speciesCol] <- as.factor(gsub(".", "_", as.character(plotData[, plotData_speciesCol]), fixed = TRUE))

# Columns to be used as extra covariates in the model
extraCovariates_plotLevel <- c("temp.c", "p.mm")
# Replace the "woody" functional trait
#traitData[, traitData_funGroup] <- as.factor(ifelse(
#  as.character(traitData[, traitData_funGroup]) == "woody",
#  "forb", as.character(traitData[, traitData_funGroup])))
colnames(traitData)[tolower(colnames(traitData)) %in% extraCovariates_plotLevel] <- tolower(colnames(traitData)[tolower(colnames(traitData)) %in% extraCovariates_plotLevel])
colnames(plotData)[tolower(colnames(plotData)) %in% extraCovariates_plotLevel] <- tolower(colnames(plotData)[tolower(colnames(plotData)) %in% extraCovariates_plotLevel])

# Retrieve the full list of species names (and remove NA values)
speciesNames <- unique(c(as.character(traitData[, traitData_speciesCol]), as.character(plotData[, plotData_speciesCol])))
speciesNames <- speciesNames[!is.na(speciesNames)]
# Retrieve the list of plot IDs
plotNames <- unique(c(as.character(plotData[, plotData_plotCol])))
plotNames <- plotNames[!is.na(plotNames)]
# Retrieve the list of site names
siteNames <- unique(c(as.character(traitData[, traitData_siteCol]), as.character(plotData[, plotData_siteCol])))
siteNames <- siteNames[!is.na(siteNames)]

# Retrieve the functional group for each of the species
funGroupVec <- sapply(X = speciesNames, FUN = function(curSpecies, fullSpecFrame) {
	# Retrieve the functional groups associated with the species
	curFunGroups <- unique(as.character(fullSpecFrame[curSpecies == as.character(fullSpecFrame$Species), "FunGroup"]))
	curFunGroups <- curFunGroups[!is.na(curFunGroups)]
	if(length(curFunGroups) <= 0) {
		warning(paste("No functional group information found for species ", curSpecies, sep = ""))
		curFunGroups <- NA
	} else if(length(curFunGroups) > 1) {
		warning(paste("More than one functional group found for species ", curSpecies, ": ", paste(curFunGroups, collapse = " , "), sep = ""))
		curFunGroups <- curFunGroups[1]
	}
	curFunGroups
}, fullSpecFrame = rbind(data.frame(
	Species = c(as.character(plotData[, plotData_speciesCol]), as.character(traitData[, traitData_speciesCol])),
	FunGroup = c(as.character(plotData[, plotData_funGroup]), as.character(traitData[, traitData_funGroup]))), extraFunGroupInfo)
)
names(funGroupVec) <- speciesNames

# Create a data frame of the plot-level covariates at each site
extraCovariateFrame_plotLevel <- do.call(rbind, lapply(X = siteNames, FUN = function(curSite, fullCovFrame) {
   # Retrieve a matrix of covariate values for the current site
   siteValMatrix <- apply(X = as.matrix(fullCovFrame[as.character(fullCovFrame$Site) == curSite, colnames(fullCovFrame) != "Site"]), FUN = function(curCol) {
      curVals <- unique(curCol)
      if(length(curVals) == 0) {
         curVals <- NA
      }
      if(length(curVals) > 1) {
         warning("multiple unique values found for extra covariates at particular sites: only the first unique value will be used")
      }
      curVals[1]
   }, MARGIN = 2)
   # Set the dimensions and column and row names of the output matrix
   dim(siteValMatrix) <- c(1, length(siteValMatrix))
   colnames(siteValMatrix) <- colnames(fullCovFrame)[colnames(fullCovFrame) != "Site"]
   rownames(siteValMatrix) <- curSite
   as.data.frame(siteValMatrix)
}, fullCovFrame = cbind(
   data.frame(
      Site = c(as.character(traitData[, traitData_siteCol]), as.character(plotData[, plotData_siteCol]))
   ),
   rbind(traitData[, extraCovariates_plotLevel], plotData[, extraCovariates_plotLevel])
)))
# Create a matrix of species composition at each plot (will generate a warning for duplicate cover values of a given species at a plot)
#  Firstly iterate across the plots
plotData_speciesComposition <- t(sapply(X = plotNames, FUN = function(curPlot, speciesNames, plotData, speciesCol, plotCol, coverCol) {
	#  Secondly iterate across the species
	outCovers <- sapply(X = speciesNames, FUN = function(curSpecies, curPlot, plotData, speciesCol, plotCol, coverCol) {
		coverVals <- 0
		# Find the rows that are the current species at the current plot
		rowsToSelect <- as.character(plotData[, speciesCol]) == curSpecies & as.character(plotData[, plotCol]) == curPlot
		if(any(rowsToSelect)) {
			# Select out the current cover values for the species at the plot
			coverVals <- plotData[rowsToSelect, coverCol]
			if(length(coverVals) > 1) {
				warning(paste("more than one cover value found for species", curSpecies, "at plot", curPlot, ": computed value is the mean value found", sep = " "))
				coverVals <- mean(coverVals, na.rm = TRUE)
			}
		}
		coverVals
	}, curPlot = curPlot, plotData = plotData, speciesCol = speciesCol, plotCol = plotCol, coverCol = coverCol)
	# Normalise the cover values
	outCovers / sum(outCovers, na.rm = TRUE)
}, speciesNames = speciesNames, plotData = plotData, speciesCol = plotData_speciesCol, plotCol = plotData_plotCol, coverCol = plotData_coverCol))
rownames(plotData_speciesComposition) <- paste("traitTotal", rownames(plotData_speciesComposition), sep = "_")
# Create a new data frame for the plot data
plotData_newShape <- cbind(
	as.data.frame(matrix(NA, nrow = nrow(plotData_speciesComposition), ncol = numTraits, dimnames = list(rownames(plotData_speciesComposition), traitData_traitCols))),
	data.frame(site = sapply(X = plotNames, FUN = function(curPlot, plotInfo) {
		outSite <- NA
		# Retrieve the current plot values in the information frame
		plotRows <- as.character(plotInfo[, 1]) == curPlot
		if(any(plotRows)) {
			# Lookup the site codes
			outSite <- unique(as.character(plotInfo[plotRows, 2]))
			if(length(outSite) > 1) {
				warning(paste("more than one site value found for plot", curPlot, ": only using the first value", sep = " "))
				outSite <- outSite[1]
			}
		}
		outSite
	}, plotData[, c(plotData_plotCol, plotData_siteCol)]), row.names = rownames(plotData_speciesComposition)),
	as.data.frame(plotData_speciesComposition),
	# Create a data frame of the plot-level extra covariates to include in the model
	do.call(cbind, lapply(X = extraCovariates_plotLevel, FUN = function(curCovName, plotNames, plotData, turfIDCol) {
	   # Retrieve the covariate values for each of the plots
	   curPlotValues <- sapply(X = plotNames, FUN = function(curPlot, curCovName, plotData, turfIDCol) {
	      # Retrieve the current covariate values
	      outVals <- unique(plotData[plotData[, turfIDCol] == curPlot, curCovName], na.rm = TRUE)
	      if(length(outVals) == 0) {
	         outVals <- NA
	      } else if(length(outVals) > 1) {
	         warning(paste("more than one unique value for covariate", curCovName, "for plot", curPlot, ": only using the first value", sep = " "))
	         outVals <- outVals[1]
	      }
	      outVals
	   }, curCovName = curCovName, plotData = plotData, turfIDCol = turfIDCol)
	   # Create a data frame with the respective column and row names
	   outFrame <- data.frame(outCol = curPlotValues)
	   rownames(outFrame) <- paste("traitTotal", plotNames, sep = "_")
	   colnames(outFrame) <- curCovName
	   outFrame
	}, plotNames = plotNames, plotData = plotData, turfIDCol = plotData_plotCol))
)

# Create a matrix of species identifier variables for the trait data (used to define the random effects in the regression model)
traitData_speciesComposition <- sapply(X = speciesNames, FUN = function(curSpecies, traitData, speciesCol) {
	outVals <- ifelse(as.character(traitData[, speciesCol]) == curSpecies, 1.0, 0.0)
}, traitData = traitData, speciesCol = traitData_speciesCol)
rownames(traitData_speciesComposition) <- paste("traitMeasurement", 1: nrow(traitData), "_Site", traitData[, traitData_siteCol], sep = "")
rownames(traitData) <- rownames(traitData_speciesComposition)
# Create a new data frame for the trait data
traitData_newShape <- cbind(
	traitData[, traitData_traitCols],
	data.frame(site = traitData[, traitData_siteCol], row.names = rownames(traitData)),
	as.data.frame(traitData_speciesComposition),
	traitData[, extraCovariates_plotLevel]
)

# Create a matrix of species identifier variables for the species and site prediction
predictData_speciesComposition <- do.call(rbind, lapply(X = siteNames, FUN = function(curSite, speciesNames) {
	# Create a species indicator matrix
	speciesIndicatorMat <- t(sapply(1:length(speciesNames), FUN = function(curIndex, numSpecies) {
		outVals <- rep(0.0, numSpecies)
		outVals[curIndex] <- 1.0
		outVals
	}, numSpecies = length(speciesNames)))
	# Set the column and rows of the species indicator matrix
	colnames(speciesIndicatorMat) <- speciesNames
	rownames(speciesIndicatorMat) <- paste("Site", curSite, "_Species", speciesNames, sep = "")
	speciesIndicatorMat
}, speciesNames))
# Create a new data frame for the prediction
sitePredictVals <- rep(siteNames, rep(length(speciesNames), length(siteNames)))
predictData_newShape <- cbind(
	as.data.frame(matrix(NA, nrow = nrow(predictData_speciesComposition), ncol = numTraits, dimnames = list(rownames(predictData_speciesComposition), traitData_traitCols))),
	data.frame(site = sitePredictVals, row.names = rownames(predictData_speciesComposition)),
	as.data.frame(predictData_speciesComposition),
	as.data.frame(as.matrix(extraCovariateFrame_plotLevel)[sitePredictVals, ], row.names = rownames(predictData_speciesComposition))
)

# Iterate over each of the trait values and run the imputation model
imputationModelOutputs <- lapply(X = traitData_traitCols, FUN = function(curTrait, combinedData, modelOutputFolder, extraCovariates_plotLevel, funGroupVec) {
	cat(paste("Calculating model for trait", curTrait, "...\n", sep = " "))
	# Retrieve the species names
	speciesNames <- names(funGroupVec)
	# Retrieve the functional group frequencies
	funcGroupTable <- table(funGroupVec)
	funcGroupNames <- names(funcGroupTable)
  # See which species actually have trait information measured for them
  hasTraitInfoSpecies <- sapply(X = speciesNames, FUN = function(curSpec, curTrait, combinedData) {
     # Retrieve the trait values for the currenct species
     traitVals <- combinedData[combinedData[, curSpec] > 0, curTrait]
     # See if there are any non-NA values
     sum(ifelse(is.na(traitVals), 0, 1)) > 0
  }, curTrait = curTrait, combinedData = combinedData)
  # Also remove species that are sole representation of their functional group
	hasTraitInfoSpecies <- hasTraitInfoSpecies & !(funGroupVec %in% funcGroupNames[funcGroupTable <= 1])
	# Set the species to use in the model
	useSpecies <- speciesNames[hasTraitInfoSpecies]
	# Remove the first species in every functional group (for parameter identifiability)
	useSpecies <- useSpecies[duplicated(funGroupVec[useSpecies])]
  # Create a data frame of functional groups
  funGroupFrame <- sapply(X = funcGroupNames, FUN = function(curFunGroup, combinedData, funGroupVec) {
  	as.matrix(combinedData[, speciesNames]) %*% ifelse(curFunGroup == funGroupVec, 1, 0)
  }, combinedData = combinedData, funGroupVec = funGroupVec)
  colnames(funGroupFrame) <- funcGroupNames
  # Retrieve the names of the functional groups that have a species with at least one measurement
  useFuncGroup <- colnames(funGroupFrame)[apply(X = as.matrix(funGroupFrame[!is.na(combinedData[, curTrait]), ]), MARGIN = 2, FUN = sum) > 0.0]
  # Create the model formula
  modelFormula <- as.formula(paste(curTrait, paste(c(
  	# No intercept term (this is wrapped up on the functional group factor term)
  	"-1"
  	# Fixed intercept for each functional group
  	, paste(useFuncGroup, collapse = " + ")
  	# Fixed slope for each functional group for each covariate
  	, sapply(X = extraCovariates_plotLevel, FUN = function(curCov, funcGroupNames) {
  		paste(curCov, funcGroupNames, sep = ":", collapse = " + ")
  	}, funcGroupNames = useFuncGroup)
  	# Fixed species effect
  	, paste(useSpecies, collapse = " + ")
  ), collapse = " + "), sep = " ~ "))
	# Make a data frame of indicator variables for the species random effects
	specEffectsFrame <- cbind(
		as.data.frame(replicate(length(extraCovariates_plotLevel), 1:nrow(combinedData), simplify = FALSE), col.names = paste("spec", extraCovariates_plotLevel, sep = "")),
		data.frame(specIntercept = 1:nrow(combinedData)))
	# Fit the model using a mixed effects model (using ridge-regression offsets for the species effect)
	modelOutput <- glm(modelFormula, data = cbind(combinedData, funGroupFrame, specEffectsFrame))
	#fittedValues <- cbind(data.frame(observedValue = combinedData[, curTrait]),
	#  modelOutput$summary.fitted.values)
	# Retrieve the fitted values and the prediction error
	fittedValues <- data.frame(
		predictValue = predict(modelOutput, newdata = cbind(combinedData, funGroupFrame, specEffectsFrame)),
		observedValue = combinedData[, curTrait],
		predictionSE = predict(modelOutput, newdata = cbind(combinedData, funGroupFrame, specEffectsFrame), se.fit = TRUE)$se.fit
	)
	fittedValues <- cbind(fittedValues, data.frame(
		lowerConfInf = fittedValues$predictValue - 1.96 * fittedValues$predictionSE,
		upperConfInf = fittedValues$predictValue + 1.96 * fittedValues$predictionSE
	))
	rownames(fittedValues) <- rownames(combinedData)
	# Designated rows used for the trait totals
	isTraitTotal <- grepl("^traitTotal_", rownames(fittedValues), perl = TRUE)
	# Designated rows used for the site/species prediction
	isSiteSpecPred <- grepl("^Site", rownames(fittedValues), perl = TRUE)
	# Select out the trait table
	traitTotalTable <- fittedValues[isTraitTotal, colnames(fittedValues) != "observedValue"]
	traitTotalTable <- cbind(
		data.frame(Plot = gsub("^traitTotal_", "", rownames(traitTotalTable), perl = TRUE), row.names = rownames(traitTotalTable)),
		traitTotalTable
	)
	# Select out the site/species prediction
	siteSpecPredTable <- fittedValues[isSiteSpecPred, colnames(fittedValues) != "observedValue"]
	# Select out the trait measurements
	traitMeasurementTable <- cbind(data.frame(
		# Retrieve the site information
		site = combinedData[!isSiteSpecPred & !isTraitTotal, "site"],
		# Retrieve the species information
		species = apply(X = as.matrix(combinedData[!isSiteSpecPred & !isTraitTotal, speciesNames]), FUN = function(curRow, speciesNames) {
			paste(speciesNames[curRow > 0], collapse = ",")
		}, MARGIN = 1, speciesNames = speciesNames)
	), fittedValues[!isSiteSpecPred & !isTraitTotal, ])
	# Create new columns to hold the site and species information
	extraSiteSpecInfo <- t(sapply(X = rownames(siteSpecPredTable), FUN = function(curRowName) {
		strsplit(gsub("^Site", "", curRowName, perl = TRUE), "_Species", fixed = TRUE)[[1]]
	}))
	colnames(extraSiteSpecInfo) <- c("Site", "Species")
	siteSpecPredTable <- cbind(extraSiteSpecInfo, siteSpecPredTable)
	# Create a list with the output reshaped for convienient lookup
	outList <- list(
		glmModel = modelOutput,
		plotWeightedMean = traitTotalTable,
		speciesSitePred = siteSpecPredTable,
		traitMeasurements = traitMeasurementTable
	)
	# Delete files if they exist
	if(file.exists(paste("plotWeightedMean_", curTrait, ".csv", sep = ""))) {
	  file.remove(paste("plotWeightedMean_", curTrait, ".csv", sep = ""))
	}
	if(file.exists(paste("speciesSitePredictions_", curTrait, ".csv", sep = ""))) {
	  file.remove(paste("speciesSitePredictions_", curTrait, ".csv", sep = ""))
	}
	if(file.exists(paste("traitMeasurements_", curTrait, ".csv", sep = ""))) {
	  file.remove(paste("traitMeasurements_", curTrait, ".csv", sep = ""))
	}
	# Create files for the prediction tables
	write.csv2(traitTotalTable, file.path(modelOutputFolder, paste("plotWeightedMean_", curTrait, ".csv", sep = "")), row.names = FALSE)
	write.csv2(siteSpecPredTable, file.path(modelOutputFolder, paste("speciesSitePredictions_", curTrait, ".csv", sep = "")), row.names = FALSE)
	write.csv2(traitMeasurementTable, file.path(modelOutputFolder, paste("traitMeasurements_", curTrait, ".csv", sep = "")), row.names = FALSE)
	outList
}, combinedData = rbind(plotData_newShape, traitData_newShape, predictData_newShape), modelOutputFolder = modelOutputFolder,
	extraCovariates_plotLevel = extraCovariates_plotLevel, funGroupVec = funGroupVec)
names(imputationModelOutputs) <- traitData_traitCols
# Save the computed model output
if(file.exists(paste(modelOutputFolder, "imputationModelOutputs.rds", sep = "/"))) {
  file.remove(paste(modelOutputFolder, "imputationModelOutputs.rds", sep = "/"))
}
saveRDS(imputationModelOutputs, file = file.path(modelOutputFolder, "imputationModelOutputs.rds"))
