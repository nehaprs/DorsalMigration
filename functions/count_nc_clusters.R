#count foxd3 nd zic3 clusters from nc_Candidates_reshaped

make_input_counted <- function(input_path = "input.xlsx",
                               output_path = "input_counted.xlsx") {
  stopifnot(file.exists(input_path))
  suppressPackageStartupMessages({
    library(readxl)
    library(writexl)
    library(dplyr)
    library(purrr)
    library(stringr)
  })
  
  sheets <- readxl::excel_sheets(input_path)
  
  out <- purrr::map_dfr(sheets, function(sh) {
    df <- readxl::read_excel(input_path, sheet = sh, col_names = TRUE)
    
    # normalize to character, trim spaces, lowercase
    df_chr <- dplyr::mutate(df, dplyr::across(
      dplyr::everything(),
      ~ tolower(trimws(as.character(.x)))
    ))
    
    has_a <- purrr::map_lgl(df_chr, ~ any(.x == "foxd3", na.rm = TRUE))
    has_b <- purrr::map_lgl(df_chr, ~ any(.x == "zic3", na.rm = TRUE))
    
    cols_a_and_b <- names(df_chr)[has_a & has_b]
    cols_a_only  <- names(df_chr)[has_a & !has_b]
    
    tibble::tibble(
      name    = sh,
      `foxd3 and zic3` = paste(cols_a_and_b, collapse = ", "),
      `foxd3 only`  = paste(cols_a_only,  collapse = ", ")
    )
  })
  
  writexl::write_xlsx(out, output_path)
  invisible(out)
}
