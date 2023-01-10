library(lubridate)
library(tidyverse)

dir_data <- "data"
dir_read <- "properties"

all_chuck <- bind_rows(
  lapply(
    list.files(
      file.path(dir_data, dir_read),
      pattern = ".csv",
      full.names = TRUE
    ),
    read_csv
  )
) |>
  suppressMessages()

chuck_data <- all_chuck |>
  mutate(Date = format(mdy_hm(Date), "%Y-%m-%d"),
         Month = month(Date),
         year = year(Date))

min_date <- ymd(min(chuck_data$Date))
max_date <- ymd(max(chuck_data$Date))


# Month3 <- c("May 20", "Jul-Sep 20", "Oct-Dec 20", "Jan-Mar 21", "Apr 21", "Jul-Sep 21", "Oct-Nov 21", "Jan-Feb 22", "Jun 22")

interval <- 3 # number of MONTHS that comprise one 'primary period'

start_dates <- seq(min_date, max_date, by = paste(interval, "month"))
end_dates <- c(start_dates[-1] - 1, max_date)

timestep_df <- tibble(start_dates, end_dates) %>%
  mutate(timestep = 1:n())
timestep_df$month <- month(timestep_df$end_dates)
timestep_df$year <- year(timestep_df$end_dates)

# for each row in the merged data, insert the integer primary period timestep
chuck_data$timestep <- NA
pb <- txtProgressBar(max = nrow(chuck_data), style = 3)
for (i in 1:nrow(chuck_data)) {
  after_start <- which(timestep_df$start_dates <= chuck_data$Date[i]) %>% max
  before_end <- which(timestep_df$end_dates >= chuck_data$Date[i]) %>% min
  if (after_start == before_end) {
    # then the start and end date is contained within a primary period
    chuck_data$timestep[i] <- timestep_df$timestep[before_end]
  } # otherwise, timestep[i] will be left as NA and filtered out later
  setTxtProgressBar(pb, i)
}
close(pb)

pp_table <- tibble(
  pp_start_date = start_dates,
  pp_end_date = end_dates
) |>
  arrange(pp_start_date, pp_end_date) |>
  mutate(PPNum = 1:n())

all_chuck_data <- chuck_data |>
  mutate(PPNum = timestep) |>
  left_join(pp_table)

end_file <- "all_chuck_data.csv"
write_csv(all_chuck_data, file.path(dir_data, end_file))



