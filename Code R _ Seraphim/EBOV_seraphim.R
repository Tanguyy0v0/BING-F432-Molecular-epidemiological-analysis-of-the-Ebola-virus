# Continuous phylogeographic reconstruction of EBOV dispersal
# BIOL-F432 Spatial and Molecular Epidemiology
# Authors: Léa Troquet and Tanguy Vandermergel


# ==============================================================================
# USER-DEFINED PARAMETERS
# ==============================================================================

# Set the working directory to the folder containing the input files before running
setwd("path/to/project/folder")

mostRecentSamplingDatum <- 2025.0

mccTreeFile <- "Binome7_MCC_tree.tree"
treesFile   <- "Binome7_trees.trees"

localTreesDirectory <- "Extracted_trees"

# Use 10% burn-in if the .trees file is the raw BEAST output.
# If burn-in was already removed before this step, set burnIn <- 0.
burnIn <- 10

randomSampling <- FALSE
nberOfTreesToSample <- 200

coordinateAttributeName <- "location"

# HPD probability for uncertainty regions
prob <- 0.80

# Temporal precision used for HPD reconstructions
# 1/12 corresponds to one month.
precision <- 1 / 12


# ==============================================================================
# STEP 0: LOAD REQUIRED PACKAGES
# ==============================================================================

required_packages <- c(
  "seraphim",
  "ape",
  "sp",
  "raster",
  "maps",
  "RColorBrewer",
  "diagram"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Package '", pkg, "' is not installed. Please install it before running this script."
      )
    )
  }
}

library(seraphim)
library(ape)
library(sp)
library(raster)
library(maps)
library(RColorBrewer)
library(diagram)


# ==============================================================================
# STEP 1: EXTRACT SPATIO-TEMPORAL INFORMATION FROM THE MCC TREE
# ==============================================================================

mccTreeExtraction_binome7 <- function(mcc_tre, mostRecentSamplingDatum) {
  
  mcc_tab <- matrix(NA, nrow = nrow(mcc_tre$edge), ncol = 7)
  
  colnames(mcc_tab) <- c(
    "node1", "node2", "length",
    "startLon", "startLat", "endLon", "endLat"
  )
  
  mcc_tab[, c("node1", "node2")] <- mcc_tre$edge
  mcc_tab[, "length"] <- mcc_tre$edge.length
  
  # Extract coordinates of descendant nodes from BEAST annotations.
  loc_annos <- mcc_tre$annotations[
    sapply(mcc_tre$annotations, function(x) {
      "location1" %in% names(x) && "location2" %in% names(x)
    })
  ]
  
  for (i in seq_len(min(nrow(mcc_tab), length(loc_annos)))) {
    
    anno <- loc_annos[[i]]
    
    mcc_tab[i, c("endLon", "endLat")] <- c(
      as.numeric(anno$location1),
      as.numeric(anno$location2)
    )
  }
  
  # Infer coordinates of parent nodes from descendant coordinates.
  for (i in seq_len(nrow(mcc_tab))) {
    
    parent <- which(mcc_tab[, "node2"] == mcc_tab[i, "node1"])
    
    if (length(parent) > 0) {
      
      mcc_tab[i, c("startLon", "startLat")] <-
        mcc_tab[parent[1], c("endLon", "endLat")]
      
    } else if (!is.null(mcc_tre$root.annotation$location1) &&
               !is.null(mcc_tre$root.annotation$location2)) {
      
      mcc_tab[i, c("startLon", "startLat")] <- c(
        as.numeric(mcc_tre$root.annotation$location1),
        as.numeric(mcc_tre$root.annotation$location2)
      )
    }
  }
  
  # Convert branch lengths into calendar years.
  dist_to_tip <- node.depth.edgelength(mcc_tre)
  max_depth <- max(dist_to_tip)
  
  mcc_tab <- as.data.frame(mcc_tab)
  
  mcc_tab$startYear <- mostRecentSamplingDatum -
    (max_depth - dist_to_tip[mcc_tab$node1])
  
  mcc_tab$endYear <- mostRecentSamplingDatum -
    (max_depth - dist_to_tip[mcc_tab$node2])
  
  return(mcc_tab)
}


# Read MCC tree and extract branch-level spatio-temporal information.
mcc_tre <- readAnnotatedNexus(mccTreeFile)

mcc_tab <- mccTreeExtraction_binome7(
  mcc_tre = mcc_tre,
  mostRecentSamplingDatum = mostRecentSamplingDatum
)

write.csv(
  mcc_tab,
  "Binome7_MCC.csv",
  row.names = FALSE,
  quote = FALSE
)

# Check missing coordinates.
missing_coordinates <- colSums(
  is.na(mcc_tab[, c("startLon", "startLat", "endLon", "endLat")])
)

print(missing_coordinates)


# ==============================================================================
# STEP 2: EXTRACT POSTERIOR TREES
# ==============================================================================

allTrees <- scan(
  file = treesFile,
  what = "",
  sep = "\n",
  quiet = TRUE
)

treeExtractions(
  localTreesDirectory = localTreesDirectory,
  allTrees = allTrees,
  burnIn = burnIn,
  randomSampling = randomSampling,
  nberOfTreesToSample = nberOfTreesToSample,
  mostRecentSamplingDatum = mostRecentSamplingDatum,
  coordinateAttributeName = coordinateAttributeName
)


# ==============================================================================
# STEP 3: ESTIMATE HPD REGIONS THROUGH TIME
# ==============================================================================

nberOfExtractionFiles <- nberOfTreesToSample

# Start of the inferred dispersal history.
startDatum <- min(mcc_tab$startYear, na.rm = TRUE)

polygons <- suppressWarnings(
  spreadGraphic2(
    localTreesDirectory = localTreesDirectory,
    nberOfExtractionFiles = nberOfExtractionFiles,
    prob = prob,
    startDatum = startDatum,
    precision = precision
  )
)


# ==============================================================================
# STEP 4: CORRECT HPD POLYGON COORDINATES IF NEEDED
# ==============================================================================

swap_spdf_coordinates <- function(spdf) {
  
  # In this dataset, HPD polygon coordinates were returned as latitude/longitude.
  # They are swapped here to obtain the expected longitude/latitude format for mapping.
  
  spdf_data <- spdf@data
  
  new_polygons <- lapply(spdf@polygons, function(poly_group) {
    
    new_poly_list <- lapply(poly_group@Polygons, function(poly) {
      
      coords <- poly@coords
      coords_swapped <- coords[, c(2, 1), drop = FALSE]
      
      Polygon(
        coords_swapped,
        hole = poly@hole
      )
    })
    
    Polygons(
      new_poly_list,
      ID = poly_group@ID
    )
  })
  
  new_sp <- SpatialPolygons(
    new_polygons,
    proj4string = spdf@proj4string
  )
  
  new_spdf <- SpatialPolygonsDataFrame(
    new_sp,
    data = spdf_data,
    match.ID = FALSE
  )
  
  return(new_spdf)
}

# Apply coordinate correction to all HPD polygons.
if (length(polygons) > 0) {
  polygons <- lapply(polygons, swap_spdf_coordinates)
}

# Assign dates to HPD polygons.
polygon_dates <- seq(
  from = startDatum,
  by = precision,
  length.out = length(polygons)
)

names(polygons) <- polygon_dates

cat("Number of HPD polygons available:", length(polygons), "\n")

if (length(polygons) > 0) {
  cat("Bounding box of the first HPD polygon after coordinate correction:\n")
  print(bbox(polygons[[1]]))
}


# ==============================================================================
# STEP 5: DEFINE TEMPORAL COLOUR SCALE
# ==============================================================================

colour_scale <- colorRampPalette(
  brewer.pal(11, "RdYlGn")
)(141)[21:121]

# Minimum and maximum years used for the temporal colour scale.
minYear <- 2023.5
maxYear <- mostRecentSamplingDatum

get_colour_index <- function(year) {
  
  idx <- round((((year - minYear) / (maxYear - minYear)) * 100) + 1)
  idx[idx < 1] <- 1
  idx[idx > 101] <- 101
  
  return(idx)
}

# Colours for MCC nodes.
endYears_indices <- get_colour_index(mcc_tab[, "endYear"])
endYears_colours <- colour_scale[endYears_indices]

# Colours for HPD polygons.
polygons_colours <- rep(NA, length(polygons))

if (length(polygons) > 0) {
  
  for (i in seq_along(polygons)) {
    
    poly_date <- as.numeric(names(polygons)[i])
    polygon_index <- get_colour_index(poly_date)
    
    polygons_colours[i] <- paste0(colour_scale[polygon_index], "70")
  }
}


# ==============================================================================
# STEP 6: MAP HPD REGIONS AND MCC TRAJECTORIES
# ==============================================================================

# Keep only MCC branches with complete coordinates.
mcc_tab_plot <- mcc_tab[
  !is.na(mcc_tab$startLon) &
    !is.na(mcc_tab$startLat) &
    !is.na(mcc_tab$endLon) &
    !is.na(mcc_tab$endLat),
]

# Background raster based on MCC coordinates.
template_raster <- raster(
  xmn = min(c(mcc_tab$startLon, mcc_tab$endLon), na.rm = TRUE) - 1,
  xmx = max(c(mcc_tab$startLon, mcc_tab$endLon), na.rm = TRUE) + 1,
  ymn = min(c(mcc_tab$startLat, mcc_tab$endLat), na.rm = TRUE) - 1,
  ymx = max(c(mcc_tab$startLat, mcc_tab$endLat), na.rm = TRUE) + 1,
  nrows = 100,
  ncols = 100
)

values(template_raster) <- NA

# Define the four temporal panels.
time_slices <- list(
  c(2023.5, 2024.0),
  c(2023.5, 2024.5),
  c(2023.5, 2025.0),
  c(2023.5, 2025.0)
)

slice_titles <- c(
  "2023.5 - 2024.0",
  "2023.5 - 2024.5",
  "2023.5 - 2025.0",
  "Cumulative 2023.5 - 2025.0"
)

# Save figure as PDF for reproducibility.
pdf(
  file = "EBOV_phylogeographic_reconstruction.pdf",
  width = 10,
  height = 8
)

par(
  mfrow = c(2, 2),
  mar = c(1.5, 1.5, 2, 1),
  oma = c(5, 3.5, 1.5, 1),
  mgp = c(0, 0.4, 0),
  lwd = 0.2,
  bty = "o"
)

for (p in seq_along(time_slices)) {
  
  slice_start <- time_slices[[p]][1]
  slice_end <- time_slices[[p]][2]
  
  # Select MCC branches for each temporal panel.
  if (p < 4) {
    
    tab_slice <- mcc_tab_plot[
      mcc_tab_plot$endYear >= slice_start &
        mcc_tab_plot$endYear < slice_end,
    ]
    
  } else {
    
    tab_slice <- mcc_tab_plot[
      mcc_tab_plot$endYear >= slice_start &
        mcc_tab_plot$endYear <= slice_end,
    ]
  }
  
  # Add the oldest branch to the first panel.
  if (p == 1 && nrow(mcc_tab_plot) > 0) {
    
    oldest_idx <- which.min(mcc_tab_plot$endYear)
    tab_slice <- unique(rbind(tab_slice, mcc_tab_plot[oldest_idx, ]))
  }
  
  plot(
    template_raster,
    col = "white",
    box = FALSE,
    axes = FALSE,
    colNA = "grey95",
    legend = FALSE,
    main = slice_titles[p],
    cex.main = 0.9
  )
  
  # Add HPD polygons for the corresponding temporal slice.
  if (length(polygons) > 0) {
    
    for (i in seq_along(polygons)) {
      
      poly_date <- as.numeric(names(polygons)[i])
      
      if (p < 4) {
        keep_poly <- poly_date >= slice_start & poly_date < slice_end
      } else {
        keep_poly <- poly_date >= slice_start & poly_date <= slice_end
      }
      
      if (keep_poly) {
        
        plot(
          polygons[[i]],
          axes = FALSE,
          col = polygons_colours[i],
          add = TRUE,
          border = adjustcolor("gray20", alpha.f = 0.20)
        )
      }
    }
  }
  
  # Add country borders.
  map(
    "world",
    add = TRUE,
    fill = FALSE,
    col = "gray20",
    lwd = 0.25
  )
  
  # Add MCC branches.
  if (nrow(tab_slice) > 0) {
    
    for (i in seq_len(nrow(tab_slice))) {
      
      curvedarrow(
        cbind(tab_slice[i, "startLon"], tab_slice[i, "startLat"]),
        cbind(tab_slice[i, "endLon"], tab_slice[i, "endLat"]),
        arr.length = 0,
        arr.width = 0,
        lwd = 0.1,
        lty = 1,
        lcol = adjustcolor("gray10", alpha.f = 0.55),
        arr.col = NA,
        arr.pos = FALSE,
        curve = 0.1,
        dr = NA,
        endhead = FALSE
      )
    }
    
    # Add MCC node points.
    point_cols <- colour_scale[get_colour_index(tab_slice$endYear)]
    
    points(
      tab_slice[, "endLon"],
      tab_slice[, "endLat"],
      pch = 21,
      bg = point_cols,
      col = "black",
      cex = 1.25,
      lwd = 0.6
    )
  }
  
  # Add root point on the cumulative panel.
  if (p == 4) {
    
    root_idx <- which(!is.na(mcc_tab$startLon) & !is.na(mcc_tab$startLat))[1]
    
    points(
      mcc_tab[root_idx, "startLon"],
      mcc_tab[root_idx, "startLat"],
      pch = 16,
      col = colour_scale[1],
      cex = 0.9
    )
    
    points(
      mcc_tab[root_idx, "startLon"],
      mcc_tab[root_idx, "startLat"],
      pch = 1,
      col = "black",
      cex = 0.9
    )
  }
  
  # Add map frame.
  rect(
    xmin(template_raster),
    ymin(template_raster),
    xmax(template_raster),
    ymax(template_raster),
    xpd = TRUE,
    lwd = 0.3
  )
  
  # Add longitude axis.
  axis(
    1,
    at = pretty(c(xmin(template_raster), xmax(template_raster)), 3),
    pos = ymin(template_raster),
    cex.axis = 0.65,
    lwd = 0,
    lwd.tick = 0.2,
    tck = -0.01,
    col.axis = "gray30"
  )
  
  # Add latitude axis.
  axis(
    2,
    at = pretty(c(ymin(template_raster), ymax(template_raster)), 3),
    pos = xmin(template_raster),
    cex.axis = 0.65,
    lwd = 0,
    lwd.tick = 0.2,
    tck = -0.01,
    col.axis = "gray30",
    las = 1
  )
}


# ==============================================================================
# STEP 7: ADD TEMPORAL COLOUR LEGEND
# ==============================================================================

par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()

legend_raster <- raster(matrix(nrow = 1, ncol = 2))
legend_raster[1] <- minYear
legend_raster[2] <- maxYear

plot(
  legend_raster,
  legend.only = TRUE,
  add = TRUE,
  col = colour_scale,
  legend.width = 0.45,
  legend.shrink = 0.35,
  smallplot = c(0.34, 0.66, 0.005, 0.03),
  horizontal = TRUE,
  legend.args = list(
    text = "Year",
    cex = 0.65,
    line = 0.45,
    col = "black"
  ),
  axis.args = list(
    cex.axis = 0.6,
    lwd = 0,
    lwd.tick = 0.3,
    tck = -0.5,
    col.axis = "black",
    line = 0,
    mgp = c(0, -0.02, 0),
    at = c(2023.5, 2024.0, 2024.5, 2025.0)
  )
)

dev.off()

cat("Phylogeographic reconstruction saved as EBOV_phylogeographic_reconstruction.pdf\n")