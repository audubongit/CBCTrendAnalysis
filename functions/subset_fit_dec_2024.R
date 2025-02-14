# model_fit is the fit object from `cmdstanr`
# target_params argument ignores square brackets (i.e,. []).
# target_params = 'alpha' will return 'alpha[1]', 'alpha[2]', etc.
# thin argument saves every X iteration in csv
subset_cmdstanr <- function(model_fit, target_params, thin = 1)
{
  #get list of temp cmdstan files
  cmdstan_files <- model_fit$output_files()
  
  #get param names
  all_params <- model_fit$metadata()$model_params
  
  #remove []
  all_params_ISB <- vapply(strsplit(all_params,
                                    split = "[", fixed = TRUE), 
                           `[`, 1, FUN.VALUE=character(1))
  
  #get idx for target params
  f_ind <- which(all_params_ISB %in% target_params)
  
  #add 6 to all indices except 1, if it exists to account for:
  #accept_stat__,stepsize__,treedepth__,n_leapfrog__,divergent__,energy__
  #which are in csv but not all_params
  #add 1 to beginning of vec if not there
  f_ind2 <- f_ind
  if (sum(f_ind2 == 1))
  {
    f_ind2[-1] <- f_ind2[-1] + 6
  } else {
    f_ind2 <- c(1, f_ind2 + 6)
  }
  
  #add 2-7 (metadata)
  f_ind2 <- append(f_ind2, 2:7, after = 1)
  
  #when feeding >0.5 million args (columns) to cut, there are issues when passing
  #too many indices. Get areas where range ('-') can be used
  #get diff between indices
  del <- diff(f_ind2)
  #add max indx to vec
  del_idx <- c(which(del > 1), length(f_ind2))
  
  #start with first var
  fv <- paste0('-f', f_ind2[1])
  #if indices are a series, inset - max idx
  if (sum(del > 1) == 0)
  {
    fv <- paste0(fv, '-', tail(f_ind2, 1))
  } else {
    while (length(del) > 0)
    {
      #find first diff > 1
      if (sum(del > 1) != 0)
      {
        md <- min(which(del > 1))
        
        # add to call
        if (md > 1)
        {
          fv <- paste0(fv, '-', f_ind2[md], ',', f_ind2[md+1])
        } else {
          fv <- paste0(fv, ',', f_ind2[md+1])
        }
        
        # reorg del and f_ind2
        del <- del[-c(1:md)]
        f_ind2 <- f_ind2[-c(1:md)]
      } else {
        #if series to end
        fv <- paste0(fv, '-', tail(f_ind2, 1))
        del <- del[-c(1:length(del))]
      }
    } 
  }
  
  #function to subset draws in csv to df
  awk_fun <- function(file, IDX, thin)
  {
    print(paste0('subsetting ', file))
    #remove header
    p1_awk_call <- paste0("gawk -F: '/^[^#]/ {print}' ", file, 
                          # select params of interest
                          " | cut -d \",\" ", IDX) 
    if (thin == 1)
    {
      awk_call <- paste0(p1_awk_call, 
                         " > ",
                         strsplit(file, '.csv')[[1]], '-draws-subset.csv')
    } else {
      #thin after keeping first line
      #https://unix.stackexchange.com/questions/648113/how-to-skip-every-three-lines-using-awk
      awk_call <- paste0(p1_awk_call, 
                         " | gawk 'NR%", thin, "==1' > ", 
                         strsplit(file, '.csv')[[1]], '-draws-subset.csv')
    }
    system(awk_call)
    
    #get header only and write to file
    call2 <- paste0("head -n 47 ", file, " > ", 
                    strsplit(file, '.csv')[[1]], '-header.csv')
    system(call2)
    
    #get elapsed time
    call3 <- paste0("tail -n 5 ", file, " > ", 
                    strsplit(file, '.csv')[[1]], '-time.csv')
    system(call3)
    
    print(paste0('combining and writing to file'))
    #combine and write to file
    awk_call4 <- paste0('cat ', 
                        strsplit(file, '.csv')[[1]], '-header.csv ', 
                        strsplit(file, '.csv')[[1]], '-draws-subset.csv ',
                        strsplit(file, '.csv')[[1]], '-time.csv > ',
                        strsplit(file, '.csv')[[1]], '-subset-comb.csv')
    system(awk_call4)
  }
  
  #run awk fun and put all chains into list
  invisible(lapply(cmdstan_files, FUN = function(x) awk_fun(file = x, 
                                                            IDX = fv,
                                                            thin = thin)))
  
  #list files in dir
  lf <- list.files(dirname(cmdstan_files)[1], full.names = TRUE)
  #get names of comb files
  new_files <- grep('-subset-comb.csv', lf, value = TRUE)
  #make sure only files pertaining to this run
  tt <- strsplit(basename(cmdstan_files), '.csv')[[1]]
  new_files2 <- grep(substring(tt, nchar(tt)-5, nchar(tt)), new_files, value = TRUE)
  
  print(paste0('creating new cmdstanr object'))
  model_fit2 <- cmdstanr::as_cmdstan_fit(new_files2)
  
  return(model_fit2)
}