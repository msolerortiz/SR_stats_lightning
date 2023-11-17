######### SETUP ####
# Package names
# packages <- c("xts", "chron", "zoo", "imputeTS", "gridExtra", "slider", "cowplot",
#               "astsa", "vars", "forecast", "R.matlab", "ggplot2", "tsibble",
#               "feasts","seasonal", "lubridate", "dplyr", "tidyr", "tibbletime",
#               "fpp3", "viridis", "stringr", "twosamples", "purrr", "Matching", "pander")

packages <- c("R.matlab", "tsibble", "cowplot")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading""
invisible(lapply(packages, library, character.only = TRUE)) # GBM

######### DATA CONFIG ####

# configuration variables
G_sample_frequency = 187
G_debug = FALSE
setwd("/home/msolerortiz/workspace_R/SR_stats_lightning")

######### FUNCTIONS ####
get_datetime_from_file <- function(file) {
  raw_info = unlist( strsplit(file, "_") )
  file_date = paste(raw_info[2],raw_info[3],raw_info[4], sep="/")
  file_datetime = paste(file_date, raw_info[5], sep=" ")
  #file_name is arrival date, 30 minutes shifting to get start date.
  date_info = ymd_hms(file_datetime) - minutes(30)
  return(date_info)
}

char_and_add_zero <-function(value) {
  if (value > 9) {
    value_as_char = as.character(value)
  } 
  
  else {
    value_as_char = paste( "0", as.character(value), sep="" )
  }
  
  return(value_as_char)  
}

#using year-month-day-hour we can find the files where the data sought might be.
#this function creates
get_search_pattern_from_datetime <- function(date) {
  #file_name is arrival date, date must be shifted backwards 30 minutes.
  date = date + minutes(30)
  yr = as.character( year(date) )
  mo = char_and_add_zero( month(date) )
  dy = char_and_add_zero( day(date) )
  ho = char_and_add_zero( hour(date) )
  file_pattern = paste("c", yr, mo, dy, ho, sep="_")
  return(file_pattern)
}

get_filename_from_datetime <- function(date, or) {
  #file_name is arrival date, date must be shifted backwards 30 minutes.
  date = date + minutes(30)
  yr = as.character( year(date) )
  mo = char_and_add_zero( month(date) )
  dy = char_and_add_zero( day(date) )
  ho = char_and_add_zero( hour(date) )
  mi = char_and_add_zero( minute(date) )
  se = char_and_add_zero( second(date) )
  file_time = paste(ho,mi,se, sep = "")
  file_date = paste("c", yr, mo, dy, file_time, or, "cal.mat", sep="_")
  return(file_date)
}

date_difference_as_samples <- function(date_post, date_prev) {
  seconds_diff = as.numeric(date_post - date_prev, units = "secs")
  samples_diff = duration_to_samples(seconds_diff)
  return(samples_diff)
}

retrieve_file_list_by_orientation <- function(file_path, orientation) {
  search_path = dir( file_path )
  is_oriented = str_detect(search_path, orientation)
  is_filtered = str_detect(search_path, "fil")
  return( search_path[is_oriented & is_filtered] )
}

load_SR_data <- function(data_path, file) {
  data <- unlist( readMat( paste( data_path, file, sep="/" ) ), use.names = FALSE )
  return( data )
}


load_matlab_as_tsibble <- function(data_path) { 
  #ToDo: remove the hardcoded 60 seconds to work with segments with different separation.
  #read .mat data into an R tsibble
  matlabFile <- readMat(data_path)
  
  dts <- dates("01/01/2016")
  tms <- times(c("00:00:00"))
  start_time <- chron(dates = dts, times = tms)
  
  data_frame <- tsibble(
    time = as_datetime(as.numeric(matlabFile$index) * 60, origin = start_time),
    nor = as.numeric(matlabFile$is.nor),
    lap = as.numeric(matlabFile$is.lap),
    log = as.numeric(matlabFile$is.log),
    avg = as.numeric(matlabFile$mean),
    std = as.numeric(matlabFile$std),
    skew = as.numeric(matlabFile$skew),
    kurt = as.numeric(matlabFile$kurt),
    index = time 
  )
  
  data_frame = data_frame %>%
    fill_gaps( .full = TRUE)  
  
  rm(matlabFile)
  return(data_frame)
}

mean_na_handling <- function(values, threshold){
  #calculates the mean of the array given as values if the % of na is under threshold
  na_ind = is.na( values )
  na_ratio = sum( na_ind ) / length(values)
  
  if ( na_ratio <= threshold ) {
    mean_value = mean(values[ !na_ind ])
  } else {
    mean_value = NA
  }
  
  return(mean_value)
}

ma_with_na_handling <- function(values, threshold, w_before, w_after) { 
  #moving average using mean_na_handling; interpolates linearly any na elements
  #left after the filtering.
  averaged = slider::slide_dbl(values, mean, .before = w_before, .after = w_after, 
                               .complete = FALSE, .size = length( values ) )
  na_ind = which( is.na(averaged) )
  
  for (val in na_ind) {
    start_val = val - w_before;
    end_val = val + w_after;
    
    if ( end_val > length(values) ) {
      end_val = length(values)
    }
    
    if ( start_val >= 1 ) {
      averaged[val] = mean_na_handling( values[ start_val:end_val ], threshold )
    }
  }
  
  averaged = na.approx(averaged, na.rm = FALSE)
  return(averaged)
}

######### FUNCTIONS DEPENDENT ON GLOBAL VARIABLES ####

duration_to_samples <- function(duration) {
  return(duration * G_sample_frequency)
}

save_fig <- function(name, type) {
  #aliases for ggsave in different sizes.
  switch(type,
         wide = ggsave(name, device= "png", path = getwd(), 
                       width = 21, height = 12, units = "cm"),
         double_wide = ggsave(name, device= "png", path = getwd(), 
                              width = 21, height = 10, units = "cm"), 
         square1 = ggsave(name, device= "png", path = getwd(), 
                          width = 11, height = 9, units = "cm"),
         square2 = ggsave(name, device= "png", path = getwd(), 
                          width = 9, height = 9, units = "cm"),
         slim = ggsave(name, device= "png", path = getwd(), 
                       width = 21, height = 4, units = "cm")
  )
}

######### DATA CONFIG  ####
na_ratio_threshold = 0.75
my_color_scale = c("springgreen", "lightskyblue", "darkolivegreen4", "lightskyblue4", "darkgreen", "darkblue")
#years_color_scale = c("gold2", "goldenrod",  "gold4", "saddlebrown", "bisque4")
years_color_scale = c("#F8766D", "#00B9E3", "#D39200", "#800080", "#00BA38")
my_elements = c("mnor", "cnor", "mmonor", "cmonor", "myenor", "cyenor")
my_labels = c("10 min raw", "1 min raw", "10 min\nmonthly", "1 min\nmonthly", "10 min\nseasonal", "1 min\nseasonal")
lagged_significant_days = c(91, 182, 274, 365, 456, 547, 639, 730)

######### LOAD DATA ####

#file names

file1minEW_path = paste(getwd(),"S_EW_1min.mat", sep ="/")
ts_1min_EW = load_matlab_as_tsibble(file1minEW_path)

file10minEW_path = paste(getwd(),"S_EW_10min.mat", sep ="/")
ts_10min_EW = load_matlab_as_tsibble(file10minEW_path)

file1minNS_path = paste(getwd(),"S_NS_1min.mat", sep ="/")
ts_1min_NS = load_matlab_as_tsibble(file1minNS_path)

file10minNS_path = paste(getwd(),"S_NS_10min.mat", sep ="/")
ts_10min_NS = load_matlab_as_tsibble(file10minNS_path)

########################### GRAPHS ####################################
######## Daily observations

# Data Frames
ts_days_main = ts_10min_EW %>%
  index_by(Datetime = ~ lubridate::floor_date(., unit = "day")) %>%
  mutate(Datetime = date(Datetime)) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  select(Datetime, nor, lap, log)

ts_days_comp = ts_1min_EW %>%
  index_by(Datetime = ~ lubridate::floor_date(., unit = "day")) %>%
  mutate(Datetime = date(Datetime)) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  select(Datetime, nor, lap, log)

ts_days_main_NS = ts_10min_NS %>%
  index_by(Datetime = ~ lubridate::floor_date(., unit = "day")) %>%
  mutate(Datetime = date(Datetime)) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  select(Datetime, nor, lap, log)

ts_days_comp_NS = ts_1min_NS %>%
  index_by(Datetime = ~ lubridate::floor_date(., unit = "day")) %>%
  mutate(Datetime = date(Datetime)) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  select(Datetime, nor, lap, log)

# Plots
# Comparison between 10min overlapped segments and 1min segments.
p1 <- ts_days_main %>% full_join( ts_days_comp, by = "Datetime" ) %>%
  filter(Datetime < "2021-01-01" ) %>%
  select(Datetime, nor.x, nor.y) %>%
  transmute(mnor = nor.x, cnor = nor.y ) %>%
  mutate( mmonor = ma_daily_month(mnor) ) %>%
  mutate( cmonor = ma_daily_month(cnor) ) %>%
  mutate( myenor = ma_daily_halfyear( mnor ) ) %>%
  mutate( cyenor = ma_daily_halfyear( cnor ) ) %>%
  gather("data_type", "value", 2:7) %>%
  mutate(data_type = factor(data_type, levels = my_elements) ) %>%
  index_by(Datetime) %>%
  ggplot(aes(x = Datetime, y = value, colour = data_type)) +
  geom_line() +
  scale_x_yearquarter(breaks = yearquarter("2016Q1") + seq(0,16, by = 4), date_minor_breaks = "3 months" ) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), limits = c(0.0, 1.0) ) +
  scale_colour_manual(values = my_color_scale,
                      name= NULL,
                      breaks = my_elements,
                      labels = my_labels
  ) +
  theme(legend.position = "top", 
        legend.box = "horizontal",
        plot.margin = unit( c(5.5,5.5,0,5.5), "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  ) +
  guides(colour = guide_legend(nrow=1)) + 
  labs(x = NULL, 
       y = "Gaussian Occurrence [%]"
  )

#p1 formatted to be saved
p1 + theme(legend.position = "bottom") + 
  scale_x_yearquarter(date_breaks = "1 year", 
                      date_minor_breaks = "3 months", 
                      name = "Date[day]")


p2 <- na_appraiser(ts_1min_NS, "day") %>% 
  filter(Datetime < "2021-01-01" ) %>%
  ggplot(aes(x = Datetime, y = 0, fill = na_percent)) +
  geom_tile() + 
  scale_fill_viridis(discrete = FALSE, name = NULL) +
  scale_x_yearquarter(breaks = yearquarter("2016Q1") + seq(0,16, by = 4), date_minor_breaks = "3 months" ) +
  scale_y_continuous(breaks = NULL, name = NULL) +
  labs(x = "Date [day]") + 
  theme(plot.margin = unit( c(0, 5.5, 5.5, 5.5), "pt"),
        legend.position = "none",
        #legend.position = "bottom", 
        #legend.box = "horizontal",
        #legend.key.height = unit(0.5, "cm"),
        #legend.title = NULL,
  )

plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(0.85,0.15) )

#Yearly moving average difference between 10 min and 1min
ts_days_main %>% full_join(ts_days_comp, by = "Datetime" ) %>%
  #filter(Datetime >= "2017-01-01", Datetime < "2018-01-01") %>%
  select(Datetime, nor.x, nor.y) %>%
  transmute(mnor = nor.x, cnor = nor.y) %>%
  mutate( myenor = ma_daily_halfyear(mnor) ) %>%
  mutate( cyenor = ma_daily_halfyear(cnor) ) %>%
  mutate(diff = myenor - cyenor) %>%
  mutate(average = (myenor + cyenor)/2) %>%
  transmute(value = diff/average) %>%
  #select(Datetime, myenor, cyenor, diff, average) %>%
  #gather("data_type", "value", 2:5) %>%
  ggplot( aes(x = Datetime, y = value) ) +
  geom_line()


#daily averaged, yearly wrapped % of distributions displayed as colour graph.
p1 <- ts_days_comp %>%
  filter(Datetime < "2021-01-01" ) %>%
  mutate( mis = 1 ) %>%
  mutate( s_mis = lap + log + nor ) %>%
  unite( "mis", s_mis, mis, sep = "/" ) %>%
  mutate( nor = lap + log + nor ) %>%
  mutate( s_nor = lap + log ) %>%
  unite( "nor", s_nor, nor, sep = "/") %>%
  mutate( lap = lap + log ) %>%
  mutate( s_lap = log ) %>%
  unite( "lap", s_lap, lap, sep = "/") %>%
  mutate( s_log = 0 ) %>%
  unite( "log", s_log, log, sep = "/") %>%
  gather("data_type", "value", 2:5) %>% 
  separate( value, c("s_val","e_val"), sep = "/" ) %>%
  mutate( s_val = as.numeric(s_val) ) %>%
  mutate( e_val = as.numeric(e_val) ) %>%
  mutate( t_year = year(Datetime) ) %>%
  ggplot( aes( x = Datetime, y = s_val, colour = data_type ) ) +
  facet_wrap(~ t_year, ncol = 1, scales = "free", strip.position = "right") +
  geom_segment( aes( xend = Datetime, yend = e_val ), size = 0.5 ) +
  scale_x_continuous(breaks = NULL,  name = NULL, expand = c(0,0) ) +
  scale_y_continuous(breaks = seq(0.25,1,0.25), limits = c(0, 1), name = NULL) +
  scale_colour_manual(values = c( "black","lightblue","yellow","darkgreen"),
                      name= NULL,
                      breaks = c("mis","nor","lap","log"),
                      labels = c("Unclassified","Gaussian","Laplacian","Logistic")
  ) + 
  theme(legend.position = "top", 
        legend.box = "horizontal",
        legend.box.margin = margin(-10,-10,-10,-10),
  )

#extra graph to place legend and x axis label at the bottom.
p2 <- ts_days_comp %>%
  filter(year(Datetime) ==  year("2019-01-01") ) %>%
  mutate( mis = 1 ) %>%
  mutate( s_mis = lap + log + nor ) %>%
  unite( "mis", s_mis, mis, sep = "/" ) %>%
  mutate( nor = lap + log + nor ) %>%
  mutate( s_nor = lap + log ) %>%
  unite( "nor", s_nor, nor, sep = "/") %>%
  mutate( lap = lap + log ) %>%
  mutate( s_lap = log ) %>%
  unite( "lap", s_lap, lap, sep = "/") %>%
  mutate( s_log = 0 ) %>%
  unite( "log", s_log, log, sep = "/") %>%
  gather("data_type", "value", 2:5) %>% 
  separate( value, c("s_val","e_val"), sep = "/" ) %>%
  mutate( s_val = as.numeric(s_val) ) %>%
  mutate( e_val = as.numeric(e_val) ) %>%
  mutate( t_year = year(Datetime) ) %>%
  ggplot( aes( x = Datetime, y = s_val, colour = data_type ) ) +
  facet_wrap(~ t_year, ncol = 1, scales = "free", strip.position = "right") +
  geom_segment( aes( xend = Datetime, yend = e_val ), size = 0.5) +
  scale_x_date(breaks = "month", date_labels = "%b", name = NULL, expand = c(0,0) ) +
  scale_y_continuous(breaks = seq(0.25,1,0.25), limits = c(0, 1), name = NULL) +
  scale_colour_manual(values = c( "black","lightblue","yellow","darkgreen"),
                      name= NULL,
                      breaks = c("mis","nor","lap","log"),
                      labels = c("Unclassified","Gaussian","Laplacian","Logistic")
  ) + 
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
  )

plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(0.95,0.05) )

#daily NA % wrapped yearly.
na_appraiser(ts_1min_EW, "day") %>% 
  mutate(t_year = year(Datetime)) %>%
  ggplot(aes(x = Datetime, y = 0, fill = na_percent)) +
  facet_wrap(~ t_year, ncol = 1, scales = "free", strip.position = "right") +
  geom_tile() + 
  scale_fill_viridis(discrete = FALSE) +
  scale_x_continuous(breaks = NULL,  name = NULL, expand = c(0,0) ) +
  scale_y_continuous(name = NULL) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
  )

p2 <- na_appraiser(ts_1min_EW, "day") %>% 
  filter(year(Datetime) == year("2019-01-01") ) %>%
  mutate(t_year = year(Datetime)) %>%
  ggplot(aes(x = Datetime, y = 0, fill = na_percent)) +
  facet_wrap(~ t_year, ncol = 1, scales = "free", strip.position = "right") +
  geom_tile() + 
  scale_fill_viridis(discrete = FALSE, name = NULL) +
  scale_x_date(breaks = "month", date_labels = "%b", name = NULL, expand = c(0,0) ) +
  scale_y_continuous(breaks = NULL, name = NULL) +
  theme(legend.position = "bottom", 
        legend.box = "horizontal",
        legend.key.height = unit(0.5, "cm"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
  )

plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(0.85,0.15) )

#ACF comparison between 10min and 1min
p1 <- ts_days_main_NS %>% full_join(ts_days_comp_NS, by = "Datetime" ) %>%
  filter(Datetime >= "2016-01-01", Datetime < "2019-01-01") %>%
  select(Datetime, nor.x, nor.y) %>%
  transmute(mmonor = ma_daily_month(nor.x), 
            cmonor = ma_daily_month(nor.y)
  ) %>%
  gather("data_type", "value", 2:3) %>%
  #filter (Datetime <= ("2019-12-31")) %>%
  mutate(data_type = factor(data_type, levels = my_elements[3:4] ) ) %>%
  ACF(value, lag_max = 365*2) %>% 
  ggplot(aes(x = lag, y = acf, colour = data_type)) +
  geom_line( size = 1.2, alpha = 0.7 ) +
  geom_line(y = 0, linetype = "dashed", colour = "black") +
  geom_line(y = 0.049, linetype = "dotdash", colour = "khaki4") +
  geom_line(y = - 0.049, linetype = "dotdash", colour = "khaki4") +
  scale_x_continuous( breaks = lagged_significant_days, name = "Lag [Day]" ) +
  scale_y_continuous( breaks = seq(-0.2, 1.0, by = 0.2), limits = c(-0.3, 1.0), name = "Auto Correlation Coefficient" ) + 
  scale_colour_manual(values = my_color_scale[3:4],
                      name= NULL,
                      breaks = my_elements[3:4],
                      labels = my_labels[3:4]
  ) + 
  theme(legend.position = "none")

p2 <- ts_days_main %>% full_join(ts_days_comp, by = "Datetime" ) %>%
  filter(Datetime >= "2016-01-01", Datetime < "2019-01-01") %>%
  select(Datetime, nor.x, nor.y) %>%
  transmute(mmonor = ma_daily_month(nor.x), 
            cmonor = ma_daily_month(nor.y)
  ) %>%
  gather("data_type", "value", 2:3) %>%
  #filter (Datetime <= ("2019-12-31")) %>%
  mutate(data_type = factor(data_type, levels = my_elements[3:4] ) ) %>%
  ACF(value, lag_max = 365*2) %>% 
  ggplot(aes(x = lag, y = acf, colour = data_type)) +
  geom_line( size = 1.2, alpha = 0.7 ) +
  geom_line(y = 0, linetype = "dashed", colour = "black") +
  geom_line(y = 0.049, linetype = "dotdash", colour = "khaki4") +
  geom_line(y = - 0.049, linetype = "dotdash", colour = "khaki4") +
  scale_x_continuous( breaks = lagged_significant_days, name = "Lag [Day]" ) +
  scale_y_continuous( breaks = seq(-0.2, 1.0, by = 0.2), limits = c(-0.3, 1.0), name = NULL ) + 
  scale_colour_manual(values = my_color_scale[3:4],
                      name= NULL,
                      breaks = my_elements[3:4],
                      labels = my_labels[3:4]
  )+
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
  )

prow <- plot_grid(p1, p2, align = "vh", nrow = 1)

legend <- get_legend(p1 + theme(
  legend.position="top",
  legend.box = "horizontal",
  legend.box.margin = margin(-10,-10,-10,-10)
)
)

plot_grid(legend, prow, nrow = 2, rel_heights = c(0.10,0.90))

######### Hourly observations

#Data Frames
ts_hours_main = ts_10min_EW %>%
  index_by(Datetime = ~ lubridate::ceiling_date(., unit = "hour")) %>%
  summarize( across(everything(),  ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  select(Datetime, nor, lap, log)

ts_hours_comp = ts_1min_EW %>%
  index_by(Datetime = ~ lubridate::ceiling_date(., unit = "hour")) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  select(Datetime, nor, lap, log)

ts_hours_comp2 = ts_1min_NS %>%
  index_by(Datetime = ~ lubridate::ceiling_date(., unit = "hour")) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  select(Datetime, nor, lap, log)

# indexed by hours with hour separated from day.
ts_day_hours_main = ts_1min_EW %>%
  index_by(time_h = ~ lubridate::floor_date(., unit = "hour")) %>%
  select(time_h, nor, lap, log) %>%
  summarize( across(everything(), mean), na.rm = TRUE ) %>%
  mutate(hour = lubridate::hour(time_h)) %>%
  as_tsibble(index = time_h, key = hour) %>%
  index_by(Datetime = ~ lubridate::floor_date(., unit = "day")) %>%
  mutate(Datetime = date(Datetime)) %>%
  group_by(Datetime, hour) %>%
  summarize( across(everything(), mean), na.rm = TRUE ) %>%
  select(Datetime, hour, nor, lap, log)

# Correlation Coefficients

corr_year_date = "2018-01-01"
ts_hours_comp_year = ts_hours_comp %>%
  filter( year(Datetime) == year(corr_year_date) )

ts_hours_comp2_year = ts_hours_comp2 %>%
  filter( year(Datetime) == year(corr_year_date) )

cor_nor = cor.test(ma_hourly_day(ts_hours_comp_year$nor), ma_hourly_day(ts_hours_comp2_year$nor))
cor_log = cor.test(ma_hourly_day(ts_hours_comp_year$log), ma_hourly_day(ts_hours_comp2_year$log))
cor_lap = cor.test(ma_hourly_day(ts_hours_comp_year$lap), ma_hourly_day(ts_hours_comp2_year$lap))


# Plots
# all data with trends wrapped by hours.
ts_day_hours_main %>%
  #filter(year(Datetime) == year("2018-01-01") | year(Datetime) == year("2019-01-01") ) %>%
  select(Datetime, hour, nor) %>%
  group_by(hour) %>%
  group_modify( ~ {
    .x %>% 
      mutate( monor = ma_daily_month( nor ) ) 
  } ) %>%
  group_modify( ~ {
    .x %>% 
      mutate( yenor = ma_daily_halfyear( nor ) ) 
  } ) %>%
  gather("data_type", "value", 3:5) %>%
  mutate(data_type = factor(data_type, levels = c("nor", "monor", "yenor" ) ) ) %>%
  ggplot(aes(x = Datetime, y = value, colour = data_type)) +
  geom_line() +
  facet_wrap(~ hour, nrow = 4) +
  #scale_x_yearquarter(breaks = NULL) + #date_breaks = "1 year", date_minor_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), limits = c(0.0, 1.0) ) +
  scale_colour_manual(values = c("red", "darkgreen", "blue"),
                      name= NULL,
                      breaks = c("nor", "monor", "yenor"),
                      labels = c("raw", "monthly MA", "season MA")
  ) +
  theme(legend.position = "top", legend.box = "horizontal") +
  guides(colour = guide_legend(nrow=1)) + 
  labs(x = NULL, 
       y = "Gaussian Occurrence [%]"
  )

# separate hours' yearly trends plotted together
p1 <- ts_day_hours_main %>%
  #filter(year(Datetime) == year("2018-01-01") | year(Datetime) == year("2019-01-01") ) %>%
  select(Datetime, hour, nor) %>%
  group_by(hour) %>%
  group_modify( ~ {
    .x %>% 
      mutate( monor = ma_daily_month( nor ) ) 
  } ) %>%
  group_modify( ~ {
    .x %>% 
      mutate( yenor = ma_daily_halfyear( nor ) ) 
  } ) %>%
  select( Datetime, hour, yenor ) %>%
  ggplot( aes( x = Datetime, y = yenor, colour = as.factor(hour) ) ) +
  geom_line() +
  scale_x_yearquarter(name = "Date [Day]", date_breaks = "1 year", date_minor_breaks = "3 months") +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), limits = c(0.0, 1.0), name = "Gaussian Occurrence [%]" ) +
  scale_colour_viridis(discrete = TRUE, name = NULL) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  guides(colour = guide_legend(nrow=2))

p2 <- ts_day_hours_main %>%
  #filter(year(Datetime) == year("2018-01-01") | year(Datetime) == year("2019-01-01") ) %>%
  as_tsibble(index = hour, key = Datetime) %>%
  index_by(hour) %>%
  summarize( yenor = mean_na_handling( nor, na_ratio_threshold ) ) %>%
  mutate( yenor = yenor/max(yenor) ) %>%
  #filter(hour >=12, hour < 18) %>%
  ggplot ( aes(x = 0, y = yenor, colour = as.factor(hour) ) ) +
  geom_point( aes(size = 3) ) +
  scale_colour_viridis(discrete = TRUE, name = NULL, breaks = NULL) +
  scale_size(name = NULL, breaks = NULL) +
  scale_y_continuous(name = NULL, breaks = NULL ) +
  scale_x_continuous(name = NULL, breaks = NULL)

grid.arrange(p1, p2, widths = c(0.95, 0.05), ncol = 2)

# Comparison between 10min overlapped segments and 1min segments.
p1 <- ts_hours_main %>% full_join( ts_hours_comp, by = "Datetime" ) %>%
  filter(yearmonth(Datetime) == yearmonth("2016-01-01") ) %>%
  select(Datetime, nor.x, nor.y) %>%
  transmute(mnor = nor.x, cnor = nor.y ) %>%
  mutate( mmonor = ma_hourly_month(mnor) ) %>%
  mutate( cmonor = ma_hourly_month(cnor) ) %>%
  mutate( myenor = ma_hourly_halfyear( mnor ) ) %>%
  mutate( cyenor = ma_hourly_halfyear( cnor ) ) %>%
  gather("data_type", "value", 2:7) %>%
  mutate(data_type = factor(data_type, levels = c("mnor", "cnor", "mmonor", "cmonor", "myenor", "cyenor"))) %>%
  index_by(Datetime) %>%
  ggplot(aes(x = Datetime, y = value, colour = data_type)) +
  geom_line() +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), limits = c(0.0, 1.0) ) +
  scale_colour_manual(values = my_color_scale,
                      name= NULL,
                      breaks = my_elements,
                      labels = my_labels
  ) +
  theme(legend.position = "top", legend.box = "horizontal") +
  guides(colour = guide_legend(nrow=1)) + 
  labs(y = "Gaussian Occurrence [%]")

p2 <- na_appraiser(ts_1min_NS, "hour") %>% 
  ggplot(aes(x = Datetime, y = 0, fill = na_percent)) +
  geom_tile() + 
  scale_fill_viridis(discrete = FALSE) +
  scale_x_yearquarter(date_breaks = "1 year", date_minor_breaks = "3 months") +
  scale_y_continuous(breaks = 0.5) +
  labs(x = "Datetime [hour]") + 
  theme(legend.position = "bottom", 
        legend.box = "horizontal",
  )

grid.arrange(p1, p2, heights = c(0.75, 0.25), nrow = 2)

#1month detail comparison
ts_hours_main %>% full_join( ts_hours_comp, by = "Datetime" ) %>%
  filter(yearmonth(Datetime) == yearmonth("2020-07-01") ) %>%
  select(Datetime, nor.x, nor.y) %>%
  transmute(mnor = nor.x, cnor = nor.y ) %>%
  gather("data_type", "value", 2:3) %>%
  mutate(data_type = factor(data_type, levels = c("mnor", "cnor", "mmonor", "cmonor", "myenor", "cyenor"))) %>%
  index_by(Datetime) %>%
  ggplot(aes(x = Datetime, y = value, colour = data_type)) +
  geom_line() +
  scale_x_datetime(date_breaks = "7 days", date_minor_breaks = "1 day", name = "Datetime [Hour]") +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), limits = c(0.0, 1.0), name = "Gaussian Occurrence [%]" ) +
  scale_colour_manual(values = my_color_scale,
                      name= NULL,
                      breaks = my_elements,
                      labels = my_labels
  ) + 
  theme(legend.position = "top", legend.box = "horizontal") +
  guides( colour = guide_legend(nrow=1) )


#hourly averaged, yearly wrapped % of distributions per month displayed as colour graph.
p1 <- ts_hours_comp %>%
  filter(Datetime < "2021-01-01" ) %>%
  filter(month(Datetime) == month("2021-11-01") ) %>%
  mutate( mis = 1 ) %>%
  mutate( s_mis = lap + log + nor ) %>%
  unite( "mis", s_mis, mis, sep = "/" ) %>%
  mutate( nor = lap + log + nor ) %>%
  mutate( s_nor = lap + log ) %>%
  unite( "nor", s_nor, nor, sep = "/") %>%
  mutate( lap = lap + log ) %>%
  mutate( s_lap = log ) %>%
  unite( "lap", s_lap, lap, sep = "/") %>%
  mutate( s_log = 0 ) %>%
  unite( "log", s_log, log, sep = "/") %>%
  gather("data_type", "value", 2:5) %>% 
  separate( value, c("s_val","e_val"), sep = "/" ) %>%
  mutate( s_val = as.numeric(s_val) ) %>%
  mutate( e_val = as.numeric(e_val) ) %>%
  mutate( t_year = year(Datetime) ) %>%
  ggplot( aes( x = Datetime, y = s_val, colour = data_type ) ) +
  facet_wrap(~ t_year, ncol = 1, scales = "free", strip.position = "right") +
  geom_segment( aes( xend = Datetime, yend = e_val ) ) +
  scale_x_continuous(breaks = NULL,  name = NULL, expand = c(0,0) ) +
  scale_y_continuous(breaks = seq(0.25,1,0.25), limits = c(0, 1), name = NULL) +
  scale_colour_manual(values = c( "black","lightblue","yellow","darkgreen"),
                      name= NULL,
                      breaks = c("mis","nor","lap","log"),
                      labels = c("Unclassified","Gaussian","Laplacian","Logistic")
  ) + 
  theme(legend.position = "top",
        legend.box = "horizontal",
        legend.box.margin = margin(-10,-10,-10,-10))

#extra graph to place legend and x axis label at the bottom.
p2 <- ts_hours_comp %>%
  filter(year(Datetime) == year("2019-01-01") ) %>%
  filter(month(Datetime) == month("2019-04-01") ) %>%
  mutate( mis = 1 ) %>%
  mutate( s_mis = lap + log + nor ) %>%
  unite( "mis", s_mis, mis, sep = "/" ) %>%
  mutate( nor = lap + log + nor ) %>%
  mutate( s_nor = lap + log ) %>%
  unite( "nor", s_nor, nor, sep = "/") %>%
  mutate( lap = lap + log ) %>%
  mutate( s_lap = log ) %>%
  unite( "lap", s_lap, lap, sep = "/") %>%
  mutate( s_log = 0 ) %>%
  unite( "log", s_log, log, sep = "/") %>%
  gather("data_type", "value", 2:5) %>% 
  separate( value, c("s_val","e_val"), sep = "/" ) %>%
  mutate( s_val = as.numeric(s_val) ) %>%
  mutate( e_val = as.numeric(e_val) ) %>%
  mutate( t_year = year(Datetime) ) %>%
  ggplot( aes( x = Datetime, y = s_val, colour = data_type ) ) +
  facet_wrap(~ t_year, ncol = 1, scales = "free", strip.position = "right") +
  geom_segment( aes( xend = Datetime, yend = e_val ) ) +
  scale_x_datetime(breaks = "day", date_labels = "%d", name = NULL, expand = c(0,0) ) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0, 1), name = NULL) +
  scale_colour_manual(values = c( "black","lightblue","yellow","darkgreen"),
                      name= NULL,
                      breaks = c("mis","nor","lap","log"),
                      labels = c("Unclassified","Gaussian","Laplacian","Logistic")
  ) + 
  theme(legend.position = "none", 
        legend.box = "horizontal",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
  )

plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(0.95,0.05) )

#hourly NA % for a month wrapped yearly.
na_appraiser(ts_1min_EW, "hour") %>% 
  filter(month(Datetime) == month("2020-10-01") ) %>%
  mutate(t_year = year(Datetime)) %>%
  ggplot(aes(x = Datetime, y = 0, fill = na_percent)) +
  facet_wrap(~ t_year, ncol = 1, scales = "free", strip.position = "right") +
  geom_tile() + 
  scale_fill_viridis(discrete = FALSE) +
  scale_x_continuous(breaks = NULL,  name = NULL, expand = c(0,0) ) +
  scale_y_continuous(name = NULL) +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),  
  )

p2 <- na_appraiser(ts_1min_EW, "hour") %>% 
  filter(year(Datetime) == month("2020-10-01") ) %>%
  filter(month(Datetime) == month("2020-10-01") ) %>%
  mutate(t_year = year(Datetime)) %>%
  ggplot(aes(x = Datetime, y = 0, fill = na_percent)) +
  facet_wrap(~ t_year, ncol = 1, scales = "free", strip.position = "right") +
  geom_tile() + 
  scale_fill_viridis(discrete = FALSE, name = NULL) +
  scale_x_datetime(breaks = "day", date_labels = "%d", name = NULL, expand = c(0,0) ) +
  scale_y_continuous(breaks = NULL, name = NULL) +
  theme(legend.position = "bottom", 
        legend.box = "horizontal",
        legend.key.height = unit(0.5, "cm"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
  )

plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(0.85,0.15) )


ACFfunction_hour <- ma_hourly_day

#ACF of months per year
p1 <- ts_hours_comp %>%
  mutate(t_year = year(Datetime), t_month = month(Datetime)) %>%
  filter(t_month >= 7, t_month < 8 ) %>%
  group_by(t_year) %>%
  group_modify( ~ {
    .x %>% 
      ACF( ma_hourly_day(nor), lag_max = 24 * 7)
  } ) %>%
  ggplot( aes(x = lag, y = acf, colour = as.factor(t_year) ) ) + 
  geom_line( size = 1.2, alpha = 0.7 ) +
  #x-axis and significance interval
  geom_line(y = 0, linetype = "dashed", colour = "black") +
  geom_line(y = 0.049, linetype = "dotdash", colour = "khaki4") +
  geom_line(y = - 0.049, linetype = "dotdash", colour = "khaki4") +
  scale_x_continuous( breaks = seq(12,24*7,12), name = "lag [Hour]" ) + 
  scale_y_continuous( breaks = seq(-0.5, 1.0, by = 0.25), limits = c(-0.6, 1.0), name = "Auto Correlation Coefficient" ) +
  scale_colour_manual(values = years_color_scale,
                      name= NULL,
                      breaks = c("2016","2017","2018","2019","2020"),
                      #                      labels = c("2016-2018","2017-2019","2018-2020")
  ) + 
  theme(
    legend.position="none",
  )

p2 <- ts_hours_comp %>%
  mutate(t_year = year(Datetime), t_month = month(Datetime)) %>%
  filter(t_month >= 11, t_month < 12 ) %>%
  group_by(t_year) %>%
  group_modify( ~ {
    .x %>% 
      ACF( ma_hourly_day(nor), lag_max = 24 * 7)
  } ) %>%
  ggplot( aes(x = lag, y = acf, colour = as.factor(t_year) ) ) + 
  geom_line( size = 1.2, alpha = 0.7 ) +
  #x-axis and significance interval
  geom_line(y = 0, linetype = "dashed", colour = "black") +
  geom_line(y = 0.049, linetype = "dotdash", colour = "khaki4") +
  geom_line(y = - 0.049, linetype = "dotdash", colour = "khaki4") +
  scale_x_continuous( breaks = seq(12,24*7,12), name = "lag [Hour]" ) + 
  scale_y_continuous( breaks = seq(-0.5, 1.0, by = 0.25), limits = c(-0.6, 1.0), name = NULL ) +
  scale_colour_manual(values = years_color_scale,
                      name= NULL,
                      breaks = c("2016","2017","2018","2019","2020"),
                      #                      labels = c("2016-2018","2017-2019","2018-2020")
  ) + 
  theme(
    legend.position="none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  )

prow <- plot_grid(p1, p2, align = "vh", nrow = 1)
prow
legend <- get_legend(p1 + theme(
  legend.position="top",
  legend.box = "horizontal"
)
)

plot_grid(legend, prow, nrow = 2, rel_heights = c(0.10,0.90))

######### Hourly averaged by month
# Data frame
ts_month_hours_EW = ts_1min_EW %>%
  index_by(time_h = ~ lubridate::ceiling_date(., "hour") ) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  index_by(hour = ~ lubridate::hour(.) ) %>%
  mutate(Datetime = lubridate::floor_date(time_h, "month")) %>%
  mutate(Datetime = yearmonth(Datetime)) %>%
  group_by(Datetime, hour) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  as_tsibble(key=Datetime, index=hour) %>%
  select(Datetime, hour, nor, lap, log) %>%
  filter(Datetime < yearmonth("2021-01-01") )

ts_month_hours_NS = ts_1min_NS %>%
  index_by(time_h = ~ lubridate::ceiling_date(., "hour") ) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  index_by(hour = ~ lubridate::hour(.) ) %>%
  mutate(Datetime = lubridate::floor_date(time_h, "month")) %>%
  mutate(Datetime = yearmonth(Datetime)) %>%
  group_by(Datetime, hour) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  as_tsibble(key=Datetime, index=hour) %>%
  select(Datetime, hour, nor, lap, log) %>%
  filter(Datetime < yearmonth("2021-01-01") )

# Plot
ts_month_hours_EW %>%
  filter(year(Datetime) == year("2020-01-01") ) %>%
  ggplot(aes(x = hour, y = nor ) ) +
  geom_line() +
  facet_wrap(~ as.Date(Datetime) ) + 
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.2), limits = c(0.0, 1.0) ) +
  scale_color_viridis(discrete = TRUE) + 
  scale_x_continuous(breaks = seq(0,23,4) )


######### Hourly averaged by year
ts_year_hours_EW = ts_month_hours_EW %>%
  index_by(year_t = lubridate::year(Datetime) ) %>%
  mutate(hour_t = hour) %>%
  group_by(year_t, hour_t) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  as_tsibble(key = year_t, index = hour_t) %>%
  select(hour_t,year_t,nor,lap,log)

ts_year_hours_NS = ts_month_hours_NS %>%
  index_by(year_t = lubridate::year(Datetime) ) %>%
  mutate(hour_t = hour) %>%
  group_by(year_t, hour_t) %>%
  summarize( across(everything(), ~ mean_na_handling(.x, na_ratio_threshold) ) ) %>%
  as_tsibble(key = year_t, index = hour_t) %>%
  select(hour_t,year_t,nor,lap,log)

# Plots
p1 <- ts_year_hours_EW %>% 
  ggplot(aes(x = hour_t, y = nor, color = as.factor(year_t) ) ) +
  geom_line( size = 1.2, alpha = 0.7 ) + 
  scale_y_continuous(breaks = seq(0.0, 0.7, 0.1), limits = c(0.0, 0.7), name = "Gaussian Occurrence [%]" ) +
  scale_x_continuous(breaks = seq(0,23,1), name = "Time [Hour]" ) +
  scale_colour_manual(values = years_color_scale,
                      name= NULL,
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10)
  )

p2 <- ts_year_hours_NS %>%
  ggplot(aes(x = hour_t, y = nor, color = as.factor(year_t) ) ) +
  geom_line( size = 1.2, alpha = 0.7 ) + 
  scale_y_continuous(breaks = seq(0.0, 0.7, 0.1), limits = c(0.0, 0.7), name = "Gaussian Occurrence [%]" ) +
  scale_x_continuous(breaks = seq(0,23,1), name = "Time [Hour]" ) +
  scale_colour_manual(values = years_color_scale,
                      name= NULL,
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10)
  ) 

prow <- plot_grid(p1, p2, align = "vh", nrow = 1)

legend_dist <- get_legend(p1 + theme(
  legend.position="top",
  legend.box = "horizontal"
)
)

plot_grid(legend_dist, prow, nrow = 2, rel_heights = c(0.1, 0.90))

#intensity and size of three major thunderstorm centers (from Schumann Resonance for Tyros - Nickolaenko and Hayakawa, 2014)
t_storms = tibble(
  hour_t = 0:23,
  Africa = c(0.45, 0.51, 0.43, 0.32, 0.24, 0.18, 0.13, 0.14, 0.22, 0.35, 0.55, 
             0.88, 1.38, 2, 2.4, 2.7, 2.47, 2.03, 1.63, 1.3, 1.02, 0.8, 0.6, 0.51),
  America = c(1.75, 1.38, 1.17, 0.91, 0.78, 0.65, 0.52, 0.43, 0.4, 0.35, 0.28, 
              0.22, 0.18, 0.17, 0.18, 0.28, 0.47, 0.87, 1.42, 1.9, 2.24, 2.45, 2.3, 2.08),
  Asia = c(0.08, 0.12, 0.15, 0.22, 0.38, 0.6, 0.88, 1.24, 1.46, 1.65, 1.58, 
           1.37, 1.14, 0.95, 0.7, 0.6, 0.5, 0.4, 0.3, 0.3, 0.2, 0.15, 0.1, 0.08)  
)

p3 <- t_storms %>% 
  full_join( ts_year_hours_EW %>% 
               index_by(hour_t) %>%
               summarize( nor_EW = mean(nor*4) ), by = "hour_t") %>%
  full_join( ts_year_hours_NS %>% 
               index_by(hour_t) %>%
               summarize( nor_NS = mean(nor*4) ), by = "hour_t") %>%
  pivot_longer(c(Africa, America, Asia, nor_EW, nor_NS), names_to = "key", values_to ="value") %>%
  ggplot( aes( x = hour_t, y = value, color = as.factor(key) ) ) +
  geom_line( size = 1.2, alpha = 0.7 ) +
  scale_colour_manual(values = c("springgreen", "darkolivegreen4", 
                                 "darkgreen", "lightskyblue", "darkblue"),
                      name= NULL,
                      breaks = c("Africa", "America", "Asia", "nor_EW", "nor_NS"),
                      labels = c("Africa (left)", "America (left)", "Asia (left)", 
                                 "EW average (right)", "NS average (right)")
  ) +
  scale_y_continuous(breaks = seq(0.0, 3, 0.5), limits = c(0.0, 3), name = "Intensity [AU]",
                     sec.axis = sec_axis(trans = ~./4, name = "Gaussian Ocurrence [%]")
  ) +
  scale_x_continuous(breaks = seq(0,23,1), name = "Time [Hour]") +
  theme(legend.position="top",
        legend.box = "horizontal",
        legend.spacing.x = unit(0.5, 'cm'),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)
  )

# spatial parameters of global thunderstorm centers (Schumann Resonance for Tyros).
t_sp = tibble(
  Af = c(10, 10, 5, 28, 8, 8),
  Am = c(4, 30, -15, -66, 12, -6),
  As = c(2.5, 10, 0, 120, 20, 0)
)

#Monthly coordinates of storm centers
lat_Af = mean( -( t_sp$Af[1] + t_sp$Af[2]*cos((1:12+3)*pi/6) + t_sp$Af[3]*cos((1:12+3)*pi/3) ) )
long_Af = mean( t_sp$Af[4] + t_sp$Af[5]*cos((1:12+3)*pi/6) + t_sp$Af[6]*cos((1:12+3)*pi/3) )

lat_Am = mean( -( t_sp$Am[1] + t_sp$Am[2]*cos((1:12+3)*pi/6) + t_sp$Am[3]*cos((1:12+3)*pi/3) ) )
long_Am = mean( t_sp$Am[4] + t_sp$Am[5]*cos((1:12+3)*pi/6) + t_sp$Am[6]*cos((1:12+3)*pi/3) )

lat_As = mean( -( t_sp$As[1] + t_sp$As[2]*cos((1:12+3)*pi/6) + t_sp$As[3]*cos((1:12+3)*pi/3) ) )
long_As = mean( t_sp$As[4] + t_sp$As[5]*cos((1:12+3)*pi/6) + t_sp$As[6]*cos((1:12+3)*pi/3) )

#with the average latitudes and longitudes of the GTC we can calculate the azimuth
#from our observatory (omnicalculator.com/other/azimuth) Coordinates (Lat, Lon)
# Observatory = (37.22, -2.55)
# African Azimuth = 11,537 km, 2.907 rad
# American Azimuth = 12,698 km, 3.548 rad
# Asian Azimuth = 15,326 km, 2.355 rad

# REAL African Azimuth = 6.141 km, 2.486 rad
# REAL American Azimuth = 7.980 km, 4.363 rad
# REAL Asian Azimuth = 13,013 km, 1.239 rad

# Coefficients for magnetic field intensity from azimuth as (NS, EW) (Nickolaenko 1997)

C_Af = c( abs( sin(2.486) ), abs( cos(2.486) ) )
C_Am = c( abs( sin(4.363) ), abs( cos(4.363) ) )
C_As = c( abs( sin(1.239) ), abs( cos(1.239) ) )

#Corrected intensity - North south comparison
p1 <- t_storms %>% 
  mutate(Africa = Africa * C_Af[1], America = America * C_Am[1], Asia = Asia * C_As[1] ) %>%
  full_join( ts_year_hours_NS %>% 
               index_by(hour_t) %>%
               summarize( nor_NS = mean(nor*4) ), by = "hour_t") %>%
  pivot_longer(c(Africa, America, Asia, nor_NS), names_to = "key", values_to ="value") %>%
  ggplot( aes( x = hour_t, y = value, color = as.factor(key) ) ) +
  geom_line( size = 1.2, alpha = 0.7 ) +
  scale_colour_manual(values = c("springgreen", "darkolivegreen4", 
                                 "darkgreen", "darkblue"),
                      name= NULL,
                      breaks = c("Africa", "America", "Asia", "nor_NS"),
                      labels = c("Africa (left)", "America (left)", "Asia (left)", 
                                 "NS average (right)")
  ) +
  scale_y_continuous(breaks = seq(0.0, 2.5, 0.5), limits = c(0.0, 2.5), name = NULL,
                     sec.axis = sec_axis(trans = ~./4, name = "Gaussian Ocurrence [%]")
  ) +
  scale_x_continuous(breaks = seq(0,23,1), name = "Time [Hour]") +
  theme(legend.position = "none",
        #legend.position="top",
        #legend.box = "horizontal",
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 10)
  )

#Corrected intensity - East West comparison
p2 <- t_storms %>% 
  mutate(Africa = Africa * C_Af[2], America = America * C_Am[2], Asia = Asia * C_As[2] * 1/2) %>%
  full_join( ts_year_hours_EW %>% 
               index_by(hour_t) %>%
               summarize( nor_EW = mean(nor*4) ), by = "hour_t") %>%
  pivot_longer(c(Africa, America, Asia, nor_EW), names_to = "key", values_to ="value") %>%
  ggplot( aes( x = hour_t, y = value, color = as.factor(key) ) ) +
  geom_line( size = 1.2, alpha = 0.7 ) +
  scale_colour_manual(values = c("springgreen", "darkolivegreen4", 
                                 "darkgreen", "lightskyblue"),
                      name= NULL,
                      breaks = c("Africa", "America", "Asia", "nor_EW"),
                      labels = c("Africa (left)", "America (left)", "Asia (left)", 
                                 "EW average (right)")
  ) +
  scale_y_continuous(breaks = seq(0.0, 2.5, 0.5), limits = c(0.0, 2.5), name = "Intensity [AU]",
                     sec.axis = sec_axis(trans = ~./4, name = NULL)
  ) +
  scale_x_continuous(breaks = seq(0,23,1), name = "Time [Hour]") +
  theme(legend.position = "none",
        #legend.position="top",
        #legend.box = "horizontal",
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 10)
  )

prow <- plot_grid(p2, p1, align = "vh", nrow = 1)

legend_dist <- get_legend(p3 + theme(
  legend.position="top",
  legend.box = "horizontal"
)
)

plot_grid(legend_dist, prow, nrow = 2, rel_heights = c(0.1, 0.90))

