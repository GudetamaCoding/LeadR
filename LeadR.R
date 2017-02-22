# LeadR - a simple Pb-210 modelling tool

# a simple Pb-210 CIC model
leadR_CIC <- function( filename, unsupp=NULL, coredate, export=FALSE ) {

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
    } else if (!exists(df$Pb.210_unsupp)) {
      stop('You need to define a means of estimating unsupported Pb-210 or define it in the data.')
    } else{
      stop('There is something wrong with df$Pb.210_unsupp and/or the definition of the unsupp estimation behaviour.')
    }

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
#leadR_CRS <- function( filename, unsupp="estimate", coredate ) {
 
  
   
#}

# a tool for generating graphs from these data
# leadR_plots <- function(  ) {

#}