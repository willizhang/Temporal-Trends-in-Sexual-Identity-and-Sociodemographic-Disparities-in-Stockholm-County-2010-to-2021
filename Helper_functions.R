# title: "Helper functions"
# author: Guoqiang Zhang
# email: guoqiang.zhang@ki.se


##### Function for Weighted Regression Models #####

### crude analysis ###

fit_model_crude <- function( formula, design ) {
  
  # Poisson regression
  tryCatch(
    svyglm( formula, design = design, family = quasipoisson( link = "log" ) ),
    error = function( e ) {
      message( e$message )
      return ( NULL )
      }
  )
  }


### adjusted analysis ###

fit_models <- function( formula, design ) {
  
  models <- list()
  
  # Poisson regression
  models$poisson <- tryCatch(
    svyglm( formula, design = design, family = quasipoisson( link = "log" ) ),
    error = function( e ) {
      message( e$message )
      return ( NULL )
      }
  )
  
  # log-binomial regression
  models$log_binomial <- tryCatch(
    svyglm( formula, design = design, family = quasibinomial( link = "log" ) ),
    error = function( e ) {
      message( e$message )
      return ( NULL )
    }
  )
  
  # logistic regression
  models$logistic <- tryCatch(
    svyglm( formula, design = design, family = quasibinomial( link = "logit" ) ),
    error = function( e ) {
      message( e$message )
      return ( NULL )
    }
  )
  
  return( models )
}



##### Function to Extract Coefficients and Confidence Intervals (on log scale) #####

### crude analysis ###

extract_results_coef_crude <- function( all_models_list, exposure, outcome, variables ) {
  
  model <- all_models_list[[ exposure ]][[ outcome ]]
  
  # check if model is available (i.e., converged)
  if ( !is.null( model ) ) {
    
    results <- cbind( coef( model )[ variables ], 
                      confint( model, ddf = degf( model$survey.design ) )[ variables, , drop = FALSE ] )
    
    return( results )
  } else {
    return( NULL )
  }
}


### adjusted analysis ###

extract_results_coef <- function( all_models_list, exposure, outcome, model_type, variables ) {
  
  model <- all_models_list[[ exposure ]][[ outcome ]][[ model_type ]]
  
  # check if model is available (i.e., converged)
  if ( !is.null( model ) ) {
    
    results <- cbind( coef( model )[ variables ], 
                      confint( model, ddf = degf( model$survey.design ) )[ variables, , drop = FALSE ] )
    
    return( results )
  } else {
    return( NULL )
  }
}


### sexual identity fluidity ###

# for data SPHC-F 2010-2021
extract_results_coef_fluidity <- function( all_models_list, exposure, model_type, variables ) {
  
  model <- all_models_list[[ exposure ]][[ model_type ]]
  
  # check if model is available
  if ( !is.null( model ) ) {
    
    results <- cbind( coef( model )[ variables ], 
                      coefci( model, vcov = sandwich )[ variables, , drop = FALSE ] )
    
    return( results )
  } else {
    return( NULL )
  }
}



##### Function to Extract Proportion Ratios and Confidence Intervals for Each Model Across Exposures #####

### crude analysis ###

extract_results_pr_crude <- function( model_base, identities, exposure ) {
  results_list <- list()
  
  for ( identity in identities ) {
    
    model_name <- paste0( "fml_", identity, "_", exposure, "_crude" )
    model <- model_base[[ model_name ]]
    
    if ( !is.null( model ) ) {
      exp_model <- exp( model )
      temp_df <- as.data.frame( exp_model )
      
      colnames( temp_df ) <- c( paste0( identity, "_point_estimate_crude" ),
                                paste0( identity, "_lower_ci_crude" ),
                                paste0( identity, "_upper_ci_crude" ) )
      
      results_list[[ identity ]] <- temp_df
    }
  }
  
  names( results_list ) <- NULL
  
  combined_df <- tibble::rownames_to_column(
    do.call( cbind, results_list ),
    var = "Exposure" )
  
  return( combined_df )
}


### adjusted analysis ###

extract_results_pr <- function( model_base, model_type, identities, exposure ) {
  results_list <- list()
  
  for ( identity in identities ) {
    
    model_name <- paste0( "fml_", identity, "_", exposure, "_", model_type )
    model <- model_base[[ model_name ]]
    
    if ( !is.null( model ) ) {
      exp_model <- exp( model )
      temp_df <- as.data.frame( exp_model )
      
      colnames( temp_df ) <- c( paste0( identity, "_point_estimate" ),
                                paste0( identity, "_lower_ci" ),
                                paste0( identity, "_upper_ci" ) )
      
      results_list[[ identity ]] <- temp_df
    }
  }
  
  names( results_list ) <- NULL
  
  combined_df <- tibble::rownames_to_column(
    do.call( cbind, results_list ),
    var = "Exposure" )
  
  return( combined_df )
}



##### Function to Extract Risk Ratios and Confidence Intervals for Each Model Across Exposures #####

extract_rr_fluidity <- function( model_results ) {
  
  results_list <- list()
  
  for ( model_name in names( model_results ) ) {
    
    exp_results <- exp( model_results[[ model_name ]] )
    
    temp_df <- as.data.frame( exp_results )
    
    type <- case_when(
      grepl( "crude", model_name )   ~ "crude",
      grepl( "model_1", model_name ) ~ "model_1",
      grepl( "model_2", model_name ) ~ "model_2" )
    
    colnames( temp_df ) <- c( paste0( "point_estimate_", type ),
                              paste0( "lower_ci_", type ),
                              paste0( "upper_ci_", type ) )
    
    results_list[[ model_name ]] <- temp_df
  }
  
  names( results_list ) <- NULL
  
  combined_df <- tibble::rownames_to_column(
    do.call( cbind, results_list ),
    var = "Exposure" )
  
  return( combined_df )
}



##### Function for Estimation of Proportion of Change in Sexual Identity #####

calculate_proportions <- function( data, variable, design ) {
  prop_var <- svyby( 
    formula = ~I( sexual_identity_fluidity == "changed" ),
    by = ~get( variable ),
    design = subset( design, !is.na( sexual_identity_fluidity ) & !is.na( get( variable ) ) ),
    FUN = svyciprop,
    vartype = "ci",
    method = "beta"
    )
  
  colnames( prop_var ) <- c( "subgroup", "changed_point_estimate", "changed_lower_ci", "changed_upper_ci" )
  
  prop_var <- left_join(
    prop_var,
    data %>% 
      filter( !is.na( sexual_identity_fluidity ) & !is.na( get( variable ) ) ) %>%
      group_by( subgroup = get( variable ) ) %>%
      summarise( sample_size = n() ),
    by = "subgroup"
  ) %>%
    mutate(
      sample_size = prettyNum( sample_size, big.mark = ",", preserve.width = "none" )
    )
  
  return( prop_var )
}


##### Function to Calculate Proportions of Sexual Identities after Imputation #####

# among demographic subgroups
calc_prop_imp_subgroup <- function( implist, design, sexual_identities, demog_var, year ) {
  
  # calculate sample size
  sample_size <- lapply( implist, function( df ) {
    table( df[[ demog_var ]] ) 
  } )
  combined_freqs <- as.data.frame( Reduce( "+", sample_size ) / length( implist ) )
  colnames( combined_freqs ) <- c( "subgroup", "sample_size" )
  combined_freqs$sample_size <- prettyNum( 
    round( combined_freqs$sample_size, 0 ), big.mark = ",", preserve.width = "none" )
  names( combined_freqs )[ names( combined_freqs ) == "sample_size" ] <- paste0( "sample_size_", year )
  
  # calculate proportions
  combined_results_list <- list()
  
  for ( identity in sexual_identities ) {
    combined_summary <- summary( MIcombine( 
      with( design,
            svyby( formula = as.formula( paste0( "~ I( sexual_identity_", year, " == '", identity, "')" ) ),
                   by = as.formula( paste0( "~", demog_var ) ),
                   FUN = svyciprop,
                   method = "beta" ) ) ) )
    
    results_df <- rownames_to_column( combined_summary[ , c( "results", "(lower", "upper)" ) ], var = "subgroup" )
    colnames( results_df ) <- c( "subgroup", paste0( identity, "_point_estimate_", year ), paste0( identity, "_lower_ci_", year ), paste0( identity, "_upper_ci_", year ) )
    
    combined_results_list[[ identity ]] <- results_df
  }
  
  final_combined_df <- Reduce( function( x, y ) {
    merge( x, y, by = "subgroup" ) 
  }, 
  combined_results_list )
  
  final_combined_df <- merge( final_combined_df, combined_freqs, by = "subgroup" )
  
  return( final_combined_df )
}


# in Stockholm County
calc_prop_imp_overall <- function( design, sexual_identities, year ) {
  
  # calculate proportions
  combined_results_list <- list()
  
  for ( identity in sexual_identities ) {
    combined_summary <- summary( MIcombine( 
      with( design,
            svyciprop( formula = as.formula( paste0( "~ I( sexual_identity_", year, " == '", identity, "')" ) ),
                       method = "beta" ) ) ) )
    
    results_df <- rownames_to_column( combined_summary[ , c( "results", "(lower", "upper)" ) ], var = "subgroup" )
    colnames( results_df ) <- c( "subgroup", paste0( identity, "_point_estimate_", year ), paste0( identity, "_lower_ci_", year ), paste0( identity, "_upper_ci_", year ) )
    results_df[ 1, 1 ] <- "Stockholm County"
    
    combined_results_list[[ identity ]] <- results_df
  }
  
  final_combined_df <- Reduce( function( x, y ) {
    merge( x, y, by = "subgroup" ) 
  }, 
  combined_results_list )
  
  final_combined_df$sample_size <- prettyNum( 
    nrow( get( paste0( "d_", year ) ) ), big.mark = ",", preserve.width = "none" )
  names( final_combined_df )[ names( final_combined_df ) == "sample_size" ] <- paste0( "sample_size_", year )
  
  return( final_combined_df )
}


# by varying age cut-offs
calc_prop_imp_iteration_age <- function( min_age, sexual_identities, year ) {
  
  for ( cut_off in ( min_age ):100 ) {
    implist_iteration_age <- lapply( get( paste( "implist", year, "transformed", sep = "_" ) ), function( df ) {
      within( df, {
        age_iteration <- ifelse( age <= cut_off, "young", "old" )
        } )
      } )
  
  survey_design_iteration_age_imp <- svydesign( ids = ~ 1,
                                                strata = ~ sampling_strata_region,
                                                weights = ~ calibrated_weight,
                                                fpc = ~ no.of.population,
                                                data = imputationList( implist_iteration_age )
                                                )
  
  combined_results_list <- list()
  
  for ( identity in sexual_identities ) {
    combined_summary <- summary( MIcombine( 
      with( survey_design_iteration_age_imp,
            svyby( formula = as.formula( paste0( "~ I( sexual_identity_", year, " == '", identity, "')" ) ),
                   by = ~ age_iteration,
                   FUN = svyciprop,
                   method = "beta" ) ) ) )
    
    results_df <- rownames_to_column( combined_summary[ , "results", drop = FALSE ], var = "subgroup" )
    colnames( results_df ) <- c( "subgroup", paste0( identity, "_point_estimate_", year ) )
    
    combined_results_list[[ identity ]] <- results_df
  }
  
  final_combined_df <- Reduce( function( x, y ) {
    merge( x, y, by = "subgroup" ) 
  }, 
  combined_results_list )
  
  final_combined_df$cut_off <- cut_off
  
  results_summary <- rbind( results_summary, final_combined_df )
  }
  
  return( results_summary )
  }



##### Function for Weighted Poisson Regression Model after Imputation #####

# Poisson regression
fit_model_imp <- function( formula, design, year ) {
  
  model_result <- summary( MIcombine(
    with( design,
          svyglm( formula, family = quasipoisson( link = "log" ) ) )
    ) )
  
  model_result_selected <- rownames_to_column( model_result[ , c( "results", "(lower", "upper)" ) ], var = "subgroup" )
  
  model_result_selected[ , -1 ] <- exp( model_result_selected[ , -1 ] )
  
  colnames( model_result_selected ) <- c( "subgroup", paste0( "point_estimate_", year ), paste0( "lower_ci_", year ), paste0( "upper_ci_", year ) )
  
  return( model_result_selected )
}


# extract results
extract_model_imp <- function( model_results, exposure, sexual_identities, model_type, variables ) {
  
  pr_imp_list <- list()
  
  for( identity in sexual_identities ){
    results <- model_results[[ exposure ]][[ paste( "fml", identity, exposure, model_type, sep = "_" ) ]]
    results_selected <- subset( results, subgroup %in% variables )
    
    colnames( results_selected )[-1] <- paste( identity, colnames( results_selected )[-1], model_type, sep = "_" )
    
    pr_imp_list[[ identity ]] <- results_selected
  }
  
  combined_df <- Reduce( function( x, y ) { 
    merge( x, y, by = "subgroup" ) 
    }, 
    pr_imp_list )
  
  return( combined_df )
  }


# by sex
extract_model_imp_age_by_sex <- function( model_results, exposure, sexual_identities, model_type, variables ) {
  
  pr_imp_age_by_sex_list <- list()
  
  for( identity in sexual_identities ){
    results <- model_results[[ paste( "fml", identity, exposure, model_type, sep = "_" ) ]]
    results_selected <- subset( results, subgroup %in% variables )
    
    colnames( results_selected )[-1] <- paste( identity, colnames( results_selected )[-1], model_type, sep = "_" )
    
    pr_imp_age_by_sex_list[[ identity ]] <- results_selected
  }
  
  combined_df <- Reduce( function( x, y ) { 
    merge( x, y, by = "subgroup" ) 
  }, 
  pr_imp_age_by_sex_list )
  
  return( combined_df )
}


##### Function to Calculate Proportion of Change in Sexual Identity after Imputation #####

# among demographic subgroups
calc_prop_fluidity_imp_subgroup <- function( implist, design, demographic_vars, year ) {
  
  combined_results_list <- list()
  
  for( demog_var in demographic_vars ) {
    
    # calculate sample size
    sample_size <- lapply( implist, function( df ) {
    table( df[[ demog_var ]] ) } )
    
    combined_freqs <- as.data.frame( Reduce( "+", sample_size ) / length( implist ) )
    colnames( combined_freqs ) <- c( "subgroup", "sample_size" )
    combined_freqs$sample_size <- prettyNum( 
      round( combined_freqs$sample_size, 0 ), big.mark = ",", preserve.width = "none" )
    names( combined_freqs )[ names( combined_freqs ) == "sample_size" ] <- paste0( "sample_size_", year )
    
    # calculate proportions
    combined_summary <- summary( MIcombine( 
      with( design,
            svyby( formula = ~ I( sexual_identity_fluidity_cat == "changed" ),
                   by = as.formula( paste0( "~", demog_var ) ),
                   FUN = svyciprop,
                   method = "beta" ) ) ) )
    
    results_df <- rownames_to_column( combined_summary[ , c( "results", "(lower", "upper)" ) ], var = "subgroup" )
    colnames( results_df ) <- c( "subgroup", paste0( "changed_point_estimate_", year ), paste0( "changed_lower_ci_", year ), paste0( "changed_upper_ci_", year ) )
    
    combined_df <- merge( results_df, combined_freqs, by = "subgroup" )
    
    combined_results_list[[ demog_var ]] <- combined_df
    
    }
  
  final_combined_df <- do.call( rbind, combined_results_list )
  rownames( final_combined_df ) <- NULL
  
  return( final_combined_df )
  }


# in Stockholm County
calc_prop_fluidity_imp_overall <- function( design, year ) {
  
  combined_summary <- summary( MIcombine( 
    with( design,
          svyciprop( formula = ~ I( sexual_identity_fluidity_cat == "changed" ),
                     method = "beta" ) ) ) )
  
  results_df <- rownames_to_column( combined_summary[ , c( "results", "(lower", "upper)" ) ], var = "subgroup" )
  colnames( results_df ) <- c( "subgroup", paste0( "changed_point_estimate_", year ), paste0( "changed_lower_ci_", year ), paste0( "changed_upper_ci_", year ) )
  results_df[ 1, 1 ] <- "Stockholm County"
  
  results_df$sample_size <- prettyNum( 
    nrow( get( paste0( "d_", year ) ) ), big.mark = ",", preserve.width = "none" )
  names( results_df )[ names( results_df ) == "sample_size" ] <- paste0( "sample_size_", year )
  
  return( results_df )
}


##### Function to Extract Risk Ratio for Change in Sexual Identity after Imputation #####
extract_fluidity_model_imp <- function( model_results, exposures, model_type ) {
  
  rr_imp_list <- list()
  
  for( cat in exposures ){
    results <- model_results[[ cat ]][[ paste( "fml", cat, model_type, sep = "_" ) ]]
    results_selected <- subset( results, subgroup %in% get( paste( "exposures", cat, sep = "_" ) ) )
    
    colnames( results_selected )[-1] <- paste( colnames( results_selected )[-1], model_type, sep = "_" )
    
    rr_imp_list[[ cat ]] <- results_selected
  }
  
  combined_df <- do.call( rbind, rr_imp_list )
  rownames( combined_df ) <- NULL
  
  return( combined_df )
}
