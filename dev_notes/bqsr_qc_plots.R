
recal_data <- list.files(path="recal", pattern="recal_data.table1", full.names = TRUE) %>%
  #sample(10, replace=FALSE) %>%
  purrr::set_names() %>%
  purrr::map(function(x){
    print(x)
    #furrr::future_map(function(x){
    if(file.size(x) == 0 | file.size(x) %>% str_replace(".table1", ".table2") == 0) return(NULL)
    # Read in reports
    report1 <- x %>% 
      gsa.read.gatkreport() 
    report2 <- x %>% 
      str_replace(".table1", ".table2") %>%
      gsa.read.gatkreport()
    
    # Get recal tables
    recal1 <- report1$RecalTable1
    recal2 <- report2$RecalTable1
    
    if(!is.null(recal1) & !is.null(recal2)){
      recal_out <- recal1 %>%
        dplyr::select(ReadGroup, QualityScore, EmpiricalQuality, Observations) %>%
        mutate(type = "before") %>%
        bind_rows(recal2 %>%
                    dplyr::select(ReadGroup, QualityScore, EmpiricalQuality, Observations) %>%
                    mutate(type = "after"))
    } else{
      recal_out <- NULL
    }
    
    # Get covariates table
    cov1 <- report1$RecalTable2
    cov2 <- report2$RecalTable2
    
    if(!is.null(cov1) & !is.null(cov2)){
      cov_out <- cov1 %>%
        mutate(type = "before") %>%
        bind_rows(cov2 %>%
                    mutate(type = "after"))
    } else{
      cov_out <- NULL
    }
    out <- list(recal = recal_out,
                cov = cov_out)
    return(out)
  })

saveRDS(recal_data, "recal.rds")

#recal_data <- readRDS("recal.rds")

quals <- purrr::map(recal_data, "recal") %>%
  bind_rows(.id = "sample_id") %>%
  mutate(sample_id = sample_id %>%
           str_remove("^.*/")%>%
           str_remove("\\..*$"))  %>%
  mutate(Sequencer = case_when(
    str_detect(ReadGroup, "HLVKYDMXX") ~ "NovaSeq",
    TRUE ~ "HiSeq"
  )) %>%
  mutate(ReadGroup = paste0(Sequencer, "_", ReadGroup))


gg.bsqr_plot <- quals %>%
  group_by(ReadGroup, QualityScore, EmpiricalQuality, type) %>%
  summarise(Observations = sum(Observations)) %>%
  ggplot(aes(x=QualityScore ,y=EmpiricalQuality, alpha=log10(Observations), colour=type)) +
  geom_abline(intercept=0, slope=1, linetype=2) + 
  geom_point()+
  xlab("Reported Quality Score") +
  ylab("Empirical Quality Score") +
  scale_color_manual(values=c("before"="maroon1","after"="blue")) +
  #base_theme +
  facet_wrap(~ReadGroup) +
  coord_fixed()+
  theme(axis.text.x = element_text(angle=0, hjust = 0.5),
        legend.position = "right") +
  labs(colour = "BSQR",
       alpha = "Log10(Bases)")

gg.bsqr_plot


gg.bsqr_grouped <- quals %>%
  group_by(ReadGroup, QualityScore, EmpiricalQuality, type) %>%
  summarise(Observations = sum(Observations)) %>%
  mutate(Sequencer = case_when(
    str_detect(ReadGroup, "HLVKYDMXX") ~ "NovaSeq",
    TRUE ~ "HiSeq"
  )) %>%
  mutate(ReadGroup = paste0(Sequencer, "_", ReadGroup)) %>%
  ggplot(aes(x=QualityScore ,y=EmpiricalQuality, alpha=log10(Observations), colour=type)) +
  geom_abline(intercept=0, slope=1, linetype=2) + 
  geom_point()+
  xlab("Reported Quality Score") +
  ylab("Empirical Quality Score") +
  scale_color_manual(values=c("before"="maroon1","after"="blue")) +
  #base_theme +
  facet_wrap(~Sequencer) +
  coord_fixed()+
  theme_half_open()+
  theme(axis.text.x = element_text(angle=0, hjust = 0.5),
        legend.position = "mone") +
  labs(colour = "BSQR",
       alpha = "Log10(Bases)")

gg.bsqr_grouped

# Distribution of quality scores
gg.qual_distribution <- quals %>%
  mutate(Sequencer = case_when(
    str_detect(ReadGroup, "HLVKYDMXX") ~ "NovaSeq",
    TRUE ~ "HiSeq"
  )) %>%
  group_by(Sequencer, QualityScore, type) %>%
  summarise(Observations = sum(Observations)) %>%
  mutate(type = factor(type, levels=c("before", "after"))) %>%
  ggplot(aes(x=QualityScore, y= Observations, fill=type)) +
  geom_col(position="identity", alpha=0.7)+
  scale_fill_manual(values=c("before"="maroon1","after"="blue")) +
  theme_half_open()+
  scale_y_log10()+
  facet_grid(~Sequencer, scale="free_y") +
  theme(legend.position = "none") +
  labs(x = "Reported Quality Score",
       y = "Log10(Bases)")

gg.qual_distribution

# Which samples have the highest residuals?
library(tidymodels)
gg.sample_bqsr <- quals %>%
  #filter(type == "after") %>%
  group_by(sample_id , type,ReadGroup) %>%
  rmse(truth = EmpiricalQuality, estimate = QualityScore) %>%
  left_join(
    read_csv("sample_data/sample_info.csv" ) %>% 
      mutate(sample_id  = str_replace(sample_id, pattern=" ", replacement="")) %>%
      filter(!duplicated(sample_id))) %>%
  ggplot(aes(x = sample_id , y = .estimate, colour=ReadGroup)) +
  geom_point()+
  #geom_boxplot()+
  facet_grid(type~species, scales="free_x", space = "free") #+ 
# base_theme

gg.sample_bqsr

# Load covariate data
covs <- purrr::map(recal_data, "cov") %>%
  bind_rows(.id = "sample_id") %>%
  mutate(sample_id = sample_id %>%
           str_remove("^.*/")%>%
           str_remove("\\..*$")) %>%
  mutate(Sequencer = case_when(
    str_detect(ReadGroup, "HLVKYDMXX") ~ "NovaSeq",
    TRUE ~ "HiSeq"
  )) %>%
  mutate(ReadGroup = paste0(Sequencer, "_", ReadGroup))


# Residuals by machine cycle
gg.rmse_cycle <- covs %>%
  filter(CovariateName == "Cycle") %>%
  mutate(CovariateValue = as.numeric(CovariateValue)) %>%
  group_by(CovariateValue, type, ReadGroup, Sequencer) %>%
  rmse(truth = EmpiricalQuality, estimate = QualityScore) %>%
  ggplot(aes(x = CovariateValue, y =.estimate, colour=type))+
  geom_point() + 
  scale_colour_manual(values=c("before"="maroon1","after"="blue")) +
  theme_half_open()+
  theme(legend.position = "bottom")+
  labs(x= "Read position",
       y = "RMSE",
       colour="BQSR")

gg.rmse_cycle

# Residuals by mutation type
gg.rmse_mutations <- covs %>%
  filter(CovariateName == "Context") %>%
  group_by(sample_id, EventType, type) %>%
  rmse(truth = EmpiricalQuality, estimate = QualityScore) %>%
  ggplot(aes(x = EventType, y =.estimate))+
  geom_jitter(width=0.1, height = 0) +
  facet_grid(type~.)

gg.rmse_mutations

# Residuals for dinucleotide mutation type
gg.rmse_dinucleotide <- covs %>%
  filter(CovariateName == "Context", EventType == "M") %>%
  mutate(CovariateValue = CovariateValue %>%
           str_replace("^A", "A>")%>%
           str_replace("^C", "C>")%>%
           str_replace("^G", "G>")%>%
           str_replace("^T", "T>"))%>%
  group_by(CovariateValue, type) %>%
  rmse(truth = EmpiricalQuality, estimate = QualityScore) %>%
  ggplot(aes(x = CovariateValue, y =.estimate, colour=type))+
  geom_point()  + 
  theme_half_open()+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none")+
  #coord_flip()+
  scale_colour_manual(values=c("before"="maroon1","after"="blue")) +
  labs(x= "Mutation",
       y = "RMSE",
       colour="BQSR")

gg.rmse_dinucleotide

# Residuals for INDEL mutation type
#gg.rmse_indel <- covs %>%
#  filter(CovariateName == "Context", !EventType == "M") %>%
#  group_by(CovariateValue, EventType, type) %>%
#  rmse(truth = EmpiricalQuality, estimate = QualityScore) %>%
#  ggplot(aes(x = CovariateValue, y =.estimate, colour=EventType))+
#  geom_point() +
#  facet_grid(type~.)
#
#gg.rmse_indel

gg.bqsr_fig <- (gg.bsqr_grouped | gg.rmse_dinucleotide) / gg.qual_distribution / gg.rmse_cycle

gg.bqsr_fig

# Writ out supplementary figure
pdf("figs/supplementary/bsqr_plots.pdf", width=11, height=8,paper = "a4r")
gg.bqsr_fig
try(dev.off(), silent=TRUE)
