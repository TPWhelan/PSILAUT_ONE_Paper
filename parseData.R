parseData <- function(path, names){
  
  temp <- read.csv(path, header = FALSE, fileEncoding = 'UTF-8-BOM') %>% 
    as.matrix()
  
  
  rownames(temp) <- colnames(temp) <- names
  
  temp %>% 
    reshape2::melt() %>% 
    rowwise() %>% 
    dplyr::mutate(Var1 = stringr::str_split(Var1, pattern = '_')[[1]][2], 
                  Var2 = stringr::str_split(Var2, pattern = '_')[[1]][2]) %>% 
    
    dplyr::mutate(Var1 = dplyr::recode(Var1, 
                                       'Cont' = 'FPN',
                                       'Default' = 'DMN',
                                       'DorsAttn' = 'DAN', 
                                       'Limbic' = 'LN', 
                                       'SalVentAttn' = 'VAN', 
                                       'SomMot' = 'SMN', 
                                       'Vis' = 'VN'),
                  Var2 = dplyr::recode(Var2, 
                                       'Cont' = 'FPN',
                                       'Default' = 'DMN',
                                       'DorsAttn' = 'DAN', 
                                       'Limbic' = 'LN', 
                                       'SalVentAttn' = 'VAN', 
                                       'SomMot' = 'SMN', 
                                       'Vis' = 'VN')) %>% 
    dplyr::group_by(Var1, Var2) %>% 
    dplyr::summarise(mean_z = mean(value, na.rm = TRUE)) %>% 
    dplyr::ungroup() %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(pair = sort(c(Var1, Var2)) %>% 
                    paste(collapse = ",")) %>% 
    dplyr::group_by(pair) %>% 
    dplyr::distinct(pair, .keep_all = T) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-pair) %>% 
    return()
  
}