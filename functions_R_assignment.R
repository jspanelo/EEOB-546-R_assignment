#####  Functions for Filtering Exercise - R assignment #######

## Two functions with three arguments for filtering each file. 
## Commons: save the output into the output folder created before hand
## Names automatically each output file as 'key'_c'Chrom_number'.txt; second one includes reverse_ before 'key'

##Arguments (shared)
# data_frame: Data frame to be filtered
# Chrom_number: chromosome to be selected in each filtering step
# key: for naming the output


write_chrom <- function(data_frame, Chrom_number, key) { # Name the function and assign names for the three arguments  
  write_delim(                                           # write_delim will write each output file:
                                                           # first component for write_delim: data file 
    arrange(                                                  # arrange() will sort the input file based on 'Position'
      filter(data_frame,                                      # filter() for split the input file (data_frame): 
                                                                      # data_frame and Chrom_number are arguments from the function
             Chromosome == Chrom_number),                       # Column to be filtered by filter()
      Position),                                              # Column to be sorted by arrange()
    paste0("./output/",                                     # second component for wrire_delim: output name
                                                              # paste0() to ensemble the file name 
           key,                                                 # argument for the function write_chrom()
           "_c",                                                
           Chrom_number,                                        # second argument for the function write_chrom()
           ".txt"),                                             # to be sure is writting a .txt file
    col_names = T,                                           # third component for write_delim(): include col names 
    delim = "\t")                                            # fourth component for write_delim(): column separator is a tab
}                                                          # end


write_chrom_reverse <- function(data_frame, Chrom_number, key) {  # Name the function and assign names for the three arguments 
  write_delim(                                                    # same as write_chrom()
    arrange(                                                      # same as write_chrom()
      filter(                                                     # same as write_chrom()
        as_tibble(                                                   # input for filter() function, the input is a tibble() 
                                                                       # as_tibble() takes the output of the upcoming function to be handled as that
          lapply(data_frame, function(x) {                             # lapply() iterates through the data_frame, applying the function gsub()
            gsub("[?]", "-", x)                                        # gsub() iterates replacing all "?" by "-" 
            })),                                                       # close lapply() function
        Chromosome == Chrom_number),                               # same as write_chrom() 
      desc(Position)),                                             # Column to be sorted in descendent order using desc(), by arrange() function 
    paste0("./output/reverse_",                                    # same as write_chrom() up to the end
           key,
           "_c", 
           Chrom_number, 
           ".txt"), 
    col_names = T, 
    delim = "\t")
}                                                                   # end
