# UK Kidney Donor Risk Index (UKKDRI) and Recipient Risk Index (UKKRRI)

library(readxl)
library(dplyr)
library(ggplot2)
library(transplantr)

out_dir <- "output/data_preparation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Donor risk index (deceased donors) -------------------------------------
donors <- read_excel("data/donors.xlsx") %>%
  mutate(
    donor_htn = case_when(
      donor_htn %in% c(1, "1", TRUE, "Yes", "YES") ~ 1,
      donor_htn %in% c(0, "0", FALSE, "No", "NO")  ~ 0,
      TRUE ~ NA_real_
    ),
    donor_sex_num = if_else(sex_FM == "F", 1, 0),
    donor_cmv = case_when(
      donor_cmv %in% c(1, "1", "Positive", "POS", "Pos") ~ 1,
      donor_cmv %in% c(0, "0", "Negative", "NEG", "Neg") ~ 0,
      TRUE ~ NA_real_
    ),
    across(c(donor_age, donor_height, donor_egfr_adm, donor_los_days), as.numeric),
    UKKDRI = ukkdri(
      age    = donor_age,
      height = donor_height,
      htn    = donor_htn,
      sex    = donor_sex_num,   # 1 = female
      cmv    = donor_cmv,
      gfr    = donor_egfr_adm,
      hdays  = donor_los_days
    ),
    UKKDRI_quartile = ukkdri_q(UKKDRI, prefix = TRUE, fct = TRUE)
  )

donors %>%
  select(donor_id, donor_age, donor_egfr_adm, UKKDRI, UKKDRI_quartile) %>%
  head()
write.csv(donors, file.path(out_dir, "donors_ukkdri.csv"), row.names = FALSE)

# ---- Recipient risk index ---------------------------------------------------
recipients <- read.csv("data/recipients.csv", stringsAsFactors = FALSE) %>%
  mutate(
    UKKRRI = round(ukkrri(
      age  = Recipient.Age,
      dx   = Recipient.OnDialysis,
      wait = Recipient.Waittime,
      dm   = Recipient.Diabetes
    ), 2),
    UKKRRI_Quartile = ukkrri_q(UKKRRI, prefix = TRUE)
  )
write.csv(recipients, file.path(out_dir, "recipients_ukkrri.csv"), row.names = FALSE)

# ---- Recipient vs donor risk quartiles --------------------------------------
indices <- read.csv("data/recipient_donor_indices.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(UKKRRI_Quartile), !is.na(UKKDRI_Quartile)) %>%
  mutate(
    UKKRRI_Quartile = factor(UKKRRI_Quartile, levels = c("R1", "R2", "R3", "R4")),
    UKKDRI_Q = ntile(as.numeric(UKKDRI_Quartile), 4),
    UKKDRI_Q = factor(paste0("D", UKKDRI_Q), levels = c("D1", "D2", "D3", "D4"))
  )

p <- ggplot(count(indices, UKKRRI_Quartile, UKKDRI_Q, name = "n"),
            aes(x = UKKRRI_Quartile, y = n, fill = UKKDRI_Q)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("#e8f1fa", "#b7d4ea", "#6fa8dc", "#1f4e79"),
                    name = "Donor Risk Index") +
  labs(x = "Recipient Risk Index (R1-R4)", y = "Number of transplant recipients") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(0.9, "lines")
  )
p

ggsave(file.path(out_dir, "recipient_by_donor_risk.png"), plot = p,
       width = 7, height = 5, dpi = 300)
