# LeadR - a simple Pb-210 modelling tool

# a simple Pb-210 CIC model
leadR_CIC <- function( filename, unsupp="NULL", coredate, export=NULL ) {

  # read in some data
  if ( !file.exists(filename) ) { stop('File doesn\'t exist, check your working directory.') 
    } else { df <- read.csv(file = filename, header = TRUE) }

  # sort everything by depth
  df <- df[ , -order("Depth")]
  
  # If it is a numeric value, subtract that. If it is estimate, do that, but if it is already defined, do nothing.
  if ( is.numeric(unsupp) ) { 
      df$Pb.210_unsupp <- df$Pb.210 - unsupp 
    } else if ( unsupp == "estimate" ) {  
      df$Pb.210_unsupp <- df$Pb.210 - min(df$Pb.210)
    } else if ("Pb.210_unsupp" %in% colnames(df) == FALSE ) {
      stop('You need to define a means of estimating unsupported Pb-210 or define it in the data.')
    } else{}

  # trim everything below the lowest value
  df <- subset( df[1:which.min(df$Pb.210_unsupp ) , ] )
    
  # log transform the data
  df$Pb.210_unsupp_log <- log( df$Pb.210_unsupp )
  
  # run a linear fit model
  df <- subset( df[1:(nrow(df)-1) , ] )
  linearmodel <- lm( Pb.210_unsupp_log ~ Depth, data = df )
  
  # extract the slope co-efficient
  slope <- coef( summary(linearmodel) )[ "Depth", "Estimate" ]
  
  # set the constant half life of Pb-210
  halflife <- 22.2

  # estimate the calendar dates of the samples
  df$date <- coredate - ( (df$Depth/-1) / ( (log(2) / halflife) / slope ) )
  
  # calculate the accumulation rate - there is a better way of doing this
  accrate <- ( df$date[1] - df$date[2] ) / ( df$Depth[1] - df$Depth[2] )
  
  # extract the model error
  slope_1s <- coef( summary(linearmodel) )["Depth", "Std. Error"]
  
  # calculate the estimate specific errors
  df$date_error <- ( df$Depth/-1 ) / ( ( log(2) / halflife ) / slope_1s )
  
  # return the data and model results
  datessubset <- names(df) %in% c( "Depth", "Depth_Span", "date", "date_error" )
  dates <- df[ , datessubset ]
  inputsubset <- names(df) %in% c( "date", "date_error" )
  input <- df[ , !inputsubset ]
  raw <- read.csv(file = filename, header = TRUE)
  
  if( !is.null(export) ) {
    write.csv( df, file = export )
  } else{}
  
  return( list( "raw" = c( raw ),  "input" = c( input ), "output" = c( dates ), "accrate" = c( accrate ) ) )
  # return( list( "input" = c( df ), "accrate" = c( accrate ) ) )
  
}

# a simple Pb-210 CRS model
leadR_CRS <- function( filename, unsupp="estimate", coredate, area="UWITEC", export=NULL ) {
 
  # read in some data
  if ( !file.exists(filename) ) { stop( 'File doesn\'t exist, check your working directory.' ) 
  } else { df <- read.csv(file = filename, header = TRUE) }
  
  # sort everything by depth
  df <- df[ , -order("Depth") ]
  
  # If it is a numeric value, subtract that. If it is estimate, do that, but if it is already defined, do nothing.
  if ( is.numeric( unsupp ) ) { 
    df$Pb.210_unsupp <- df$Pb.210 - unsupp 
  } else if ( unsupp == "estimate" ) {  
    df$Pb.210_unsupp <- df$Pb.210 - min( df$Pb.210 )
  } else if ( "Pb.210_unsupp" %in% colnames(df) == FALSE ) {
    stop( 'You need to define a means of estimating unsupported Pb-210 or define it in the data.' )
  } else{}
  
  # trim everything below the lowest value
  df <- subset( df[ 1:which.min( df$Pb.210_unsupp ) , ] )

  # Multiply the unsupported Pb.210 by the dry mass
  df$Pb.210_unsupp_g <- df$Pb.210_unsupp * df$dry_mass
  
  # Calculate the area of the corer
  if ( area == "UWITEC" ) {
    area <- ( pi * 6/2 ) ^2
  } else if ( is.numeric(area) == FALSE ) {
    stop( 'You must state the diameter of the corer or state it is a standard gauge UWITEC.' )
  } else {
    area <- ( pi * area/2 ) ^2
  }
  
  # Calculate unsupported Pb-210 per unit dry weight per unit area
  df$Pb.210_unsupp_g_a <- df$Pb.210_unsupp_g / area
  
  # Calculate the total inventory
  l <- sum( df$Pb.210_unsupp_g_a )
  
  # Calculate the cumulative inventories
  df$li <- rev( cumsum( df$Pb.210_unsupp_g_a ) )
  
  # Calculate the proportional inventories
  df$l_li <- l / df$li
  
  # Log transform that data
  df$l_li_log <- log( df$l_li )
  
  # Set some constants - date of coring and decay constant
  decay_constant <- 0.03114
  
  # Express the decay constant as divided by 1
  decay_constant_div1 <- 1 / decay_constant
  
  # Calculate the ages
  df$date <- coredate - ( decay_constant_div1 * df$l_li_log )
  
  # Calculate the errors
  # Multiply the unsupported Pb.210 error by the dry mass
  df$Pb.210_1s_unsupp_g <- df$Pb.210_1s * df$dry_mass
  
  # Calculate unsupported Pb-210 error per unit dry weight per unit area
  df$Pb.210_1s_unsupp_g_a <- df$Pb.210_1s_unsupp_g / area
  
  # Calculate the total inventory error
  # Remember, errors add in quadrature (or the square root of the sum of the squares)
  # So square them first
  df$Pb.210_1s_unsupp_g_a_sq <- df$Pb.210_1s_unsupp_g_a ^2
  # And then add them all up, then square root them
  l_1s <- sqrt( sum( df$Pb.210_1s_unsupp_g_a_sq ) )
  
  # Calculate the cumulative inventories error
  # Remember, errors add in quadrature (or the square root of the sum of the squares)
  df$li_1s <- rev( cumsum( df$Pb.210_1s_unsupp_g_a_sq ) )
  df$li_1s <- sqrt( df$li_1s )
  
  # Calculate the proportional inventories error
  # Remember, when errors are multiplied, the fractional uncertainties add in quadrature.
  # Work out the fractional uncertainty of l
  l_1s_prop <- l_1s / l
  # Work out the fractional uncertainty of li
  df$li_1s_prop <- df$li_1s / df$li
  # Work out the square root of the sum of the squares of the fractional uncertainties
  df$l_li_1s_prop <- sqrt ( ( l_1s_prop ^2 ) + ( df$li_1s_prop ^2 ) )
  # Report the error of l/li
  df$l_li_1s <- df$l_li * df$l_li_1s_prop
  
  # Calculate the age error
  df$date_1s <- decay_constant_div1 * df$l_li_1s
  
  # calculate the accumulation rate
  df$accrate <- lapply( as.numeric(rownames(foob)), function(i) (foob$Depth[i] - foob$Depth[i+1]) / (foob$date[i] - foob$date[i+1]) )
  
  # return the data and model results
  datessubset <- names(df) %in% c( "Depth", "Depth_Span", "date", "date_error", "accrate" )
  dates <- df[ , datessubset ]
  inputsubset <- names(df) %in% c( "date", "date_error" )
  input <- df[ , !inputsubset ]
  raw <- read.csv(file = filename, header = TRUE)
  
  if( !is.null( export ) ) {
    write.csv( df, file = export )
  } else{}
  
  return( list( "raw" = c( raw ),  "input" = c( input ), "output" = c( dates ) ) )
}

# a tool for generating graphs from these data

#}