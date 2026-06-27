#!/usr/bin/env Rscript
# ============================================================
# REVISION re-analysis pipeline (Psychology Reports submission)
# Beyond Income: SSC, Psychological Resources & SQOL in Korea
# ------------------------------------------------------------
# Addresses ARQOL desk-reject (CMV) + AI peer-review comments.
# KEY CORRECTION: financial_worry was mis-mapped to B3 (주거형태,
#   housing type) in analysis.R. Correct item is C8 (생활비 조달
#   걱정 정도, 0-10). This script uses C8.
# Output: results/r_revision_results.json
# ============================================================

suppressMessages({
  library(lavaan)
  library(psych)
  library(jsonlite)
})

# Portable paths: run from the repository root (Rscript analysis.R).
# Place the KOSSDA data file (see README, Data availability) at ./data/korean_happiness_qol_2021.csv
PROJ <- tryCatch({ a <- commandArgs(trailingOnly=FALSE); d <- dirname(sub("--file=","",a[grep("--file=",a)])); if (length(d) && nzchar(d)) normalizePath(d) else getwd() }, error=function(e) getwd())
data_file   <- file.path(PROJ, "data", "korean_happiness_qol_2021.csv")
results_dir <- file.path(PROJ, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

cat(strrep("=", 60), "\nREVISION pipeline | lavaan", as.character(packageVersion("lavaan")), "\n", strrep("=", 60), "\n\n")

# ============================================================
# 1. LOAD & PREPROCESS
# ============================================================
df <- read.csv(data_file, fileEncoding = "EUC-KR")

vars_map <- list(
  cantril = "q1", satisfaction = "C1", yesterday_happy = "C2", eudaimonia = "C4",
  yesterday_depress = "C3", autonomy = "C5", health_physical = "C6", health_mental = "C7",
  covid_economy = "C10", subjective_class = "C13", household_income = "B1_1",
  financial_worry = "C8",                 # <- CORRECTED (was B3 = housing type)
  gender = "A1_2", birth_year = "A1_3", education = "A1_4",
  marital_status = "A1_5", employment = "A1_8",
  region = "reg3", urban_rural = "reg2", weight = "wt"
)
df2 <- data.frame(row.names = 1:nrow(df))
for (nm in names(vars_map)) if (vars_map[[nm]] %in% names(df)) df2[[nm]] <- df[[vars_map[[nm]]]]

df2$female   <- ifelse(df2$gender == 2, 1, 0)
df2$age      <- 2021 - df2$birth_year + 1
df2$married  <- ifelse(df2$marital_status == 1, 1, 0)
df2$employed <- ifelse(df2$employment %in% 1:6, 1, 0)   # CORRECTED: all in employment (was <=2 = wage workers only)
df2$yesterday_depress_r <- 10 - df2$yesterday_depress
df2$financial_worry_r   <- 10 - df2$financial_worry          # C8 reverse: higher = less worry
df2$household_income[df2$household_income == 9999] <- NA
df2$log_income <- log(df2$household_income)
df2$log_income[is.infinite(df2$log_income)] <- NA

analysis_vars <- c("cantril","satisfaction","yesterday_happy","eudaimonia",
                   "log_income","education","health_physical","health_mental","autonomy",
                   "female","age","married","employed","subjective_class","covid_economy",
                   "yesterday_depress_r","financial_worry_r","region","urban_rural","weight")

df_clean <- df2[complete.cases(df2[, analysis_vars]), ]
df_clean$age_c    <- df_clean$age - mean(df_clean$age)
df_clean$age_sq_c <- df_clean$age_c^2
N <- nrow(df_clean)
cat(sprintf("N (listwise) = %d | age range %d-%d\n", N, min(df_clean$age), max(df_clean$age)))

# Objective SES composite (income + education + employment), z-averaged
df_clean$z_income   <- as.numeric(scale(df_clean$log_income))
df_clean$z_edu      <- as.numeric(scale(df_clean$education))
df_clean$z_employed <- as.numeric(scale(df_clean$employed))
df_clean$obj_ses    <- (df_clean$z_income + df_clean$z_edu + df_clean$z_employed) / 3
alpha_obj <- suppressWarnings(psych::alpha(df_clean[, c("z_income","z_edu","z_employed")], warnings = FALSE)$total$raw_alpha)

CTRL     <- "health_mental + health_physical + autonomy + female + age_c + age_sq_c + married + employed + yesterday_depress_r + financial_worry_r + covid_economy"
CTRL_NOEMP <- "health_mental + health_physical + autonomy + female + age_c + age_sq_c + married + yesterday_depress_r + financial_worry_r + covid_economy"
MEAS     <- "happiness =~ cantril + satisfaction + yesterday_happy + eudaimonia"

getcoef <- function(fit, rhs) {
  p <- parameterEstimates(fit, standardized = TRUE)
  r <- p[p$op == "~" & p$rhs == rhs, ]
  if (nrow(r) == 0) return(list(estimate = NA, se = NA, z = NA, p = NA, std = NA))
  list(estimate = round(r$est,4), se = round(r$se,4), z = round(r$z,2),
       p = round(r$pvalue,6), std = round(r$std.all,4))
}
fitlist <- function(fit) {
  fi <- fitMeasures(fit, c("cfi","tli","rmsea","srmr","aic","bic"))
  list(CFI=round(fi["cfi"],3), TLI=round(fi["tli"],3), RMSEA=round(fi["rmsea"],3),
       SRMR=round(fi["srmr"],3), AIC=round(fi["aic"],1), BIC=round(fi["bic"],1))
}
r2of <- function(fit) round(as.numeric(inspect(fit, "r2")["happiness"]), 4)

# ============================================================
# 2. CFA + RELIABILITY
# ============================================================
cat("\n[CFA + reliability]\n")
fit_cfa <- cfa(MEAS, data = df_clean, estimator = "MLR")
fi_cfa  <- fitMeasures(fit_cfa, c("cfi","tli","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr","chisq","df","pvalue"))
ss      <- standardizedSolution(fit_cfa)
lam     <- ss[ss$op == "=~", "est.std"]
alpha_h <- suppressWarnings(psych::alpha(df_clean[, c("cantril","satisfaction","yesterday_happy","eudaimonia")])$total$raw_alpha)
ave_h   <- mean(lam^2)
omega_h <- sum(lam)^2 / (sum(lam)^2 + sum(1 - lam^2))   # McDonald's omega / CR (congeneric)
cat(sprintf("  alpha=%.3f omega/CR=%.3f AVE=%.3f | CFI=%.3f TLI=%.3f RMSEA=%.3f SRMR=%.3f\n",
            alpha_h, omega_h, ave_h, fi_cfa["cfi"], fi_cfa["tli"], fi_cfa["rmsea"], fi_cfa["srmr"]))

# ============================================================
# 3. ICC (one-way ANOVA based, no lme4 needed)
# ============================================================
cat("\n[ICC by 3-region]\n")
icc1_fun <- function(y, g) {
  a <- summary(aov(y ~ factor(g)))[[1]]
  MSB <- a["factor(g)","Mean Sq"]; MSW <- a["Residuals","Mean Sq"]
  k <- mean(table(g))
  max(0, (MSB - MSW) / (MSB + (k - 1) * MSW))
}
icc <- sapply(c("cantril","satisfaction","yesterday_happy","eudaimonia"),
              function(v) round(icc1_fun(df_clean[[v]], df_clean$region), 4))
print(icc)

# ============================================================
# 4. MULTICOLLINEARITY: correlation matrix + VIF
# ============================================================
cat("\n[VIF / correlations]\n")
pred <- c("subjective_class","log_income","health_mental","health_physical","autonomy",
          "yesterday_depress_r","financial_worry_r","covid_economy",
          "female","age_c","age_sq_c","married","employed")
vif <- sapply(pred, function(v) {
  f <- as.formula(paste(v, "~", paste(setdiff(pred, v), collapse = "+")))
  round(1 / (1 - summary(lm(f, data = df_clean))$r.squared), 3)
})
print(vif)
cormat <- round(cor(df_clean[, c("cantril","satisfaction","yesterday_happy","eudaimonia",pred)]), 3)

# ============================================================
# 5. COMPETING MODELS M1/M2/M3 (corrected)
# ============================================================
cat("\n[Competing models]\n")
fit_m1 <- sem(paste(MEAS, "\nhappiness ~ log_income +", CTRL), data=df_clean, estimator="MLR")
fit_m2 <- sem(paste(MEAS, "\nhappiness ~ subjective_class +", CTRL), data=df_clean, estimator="MLR")
fit_m3 <- sem(paste(MEAS, "\nhappiness ~ log_income + subjective_class +", CTRL), data=df_clean, estimator="MLR")
cat(sprintf("  M1 R2=%.4f | M2 R2=%.4f | M3 R2=%.4f\n", r2of(fit_m1), r2of(fit_m2), r2of(fit_m3)))
cat(sprintf("  income M1 std=%.3f p=%.4f -> M3 std=%.3f p=%.4f\n",
            getcoef(fit_m1,"log_income")$std, getcoef(fit_m1,"log_income")$p,
            getcoef(fit_m3,"log_income")$std, getcoef(fit_m3,"log_income")$p))
cat(sprintf("  SSC M2 std=%.3f -> M3 std=%.3f\n", getcoef(fit_m2,"subjective_class")$std, getcoef(fit_m3,"subjective_class")$std))

# full M3 coefficient table
p_m3 <- parameterEstimates(fit_m3, standardized = TRUE)
struct <- p_m3[p_m3$op == "~" & p_m3$lhs == "happiness", ]
struct <- struct[order(-abs(struct$std.all)), ]
m3_coefs <- lapply(1:nrow(struct), function(i) list(
  variable=struct$rhs[i], estimate=round(struct$est[i],4), se=round(struct$se[i],4),
  z=round(struct$z[i],2), p=round(struct$pvalue[i],6), std=round(struct$std.all[i],4),
  significant=struct$pvalue[i] < 0.05))

# ============================================================
# 6. CMV: ULMC (narrow + wide) + confirmatory 1-factor CFA
# ============================================================
cat("\n[CMV / ULMC]\n")
ulmc_w <- paste0(MEAS, "
  method =~ a*cantril + a*satisfaction + a*yesterday_happy + a*eudaimonia +
            a*subjective_class + a*health_mental + a*health_physical + a*autonomy +
            a*yesterday_depress_r + a*financial_worry_r + a*covid_economy
  happiness ~~ 0*method
  method ~~ 1*method
  happiness ~ log_income + subjective_class + ", CTRL)
ulmc_n <- paste0(MEAS, "
  method =~ a*cantril + a*satisfaction + a*yesterday_happy + a*eudaimonia
  happiness ~~ 0*method
  method ~~ 1*method
  happiness ~ log_income + subjective_class + ", CTRL)
fit_un <- tryCatch(sem(ulmc_n, data=df_clean, estimator="MLR"), error=function(e) NULL)
fit_uw <- tryCatch(sem(ulmc_w, data=df_clean, estimator="MLR"), error=function(e) NULL)
all_lik <- c("cantril","satisfaction","yesterday_happy","eudaimonia","subjective_class",
             "health_mental","health_physical","autonomy","yesterday_depress_r",
             "financial_worry_r","covid_economy")
fit_1f <- tryCatch(cfa(paste("F1 =~", paste(all_lik, collapse="+")), data=df_clean, estimator="MLR"), error=function(e) NULL)
if (!is.null(fit_uw)) cat(sprintf("  ULMC-wide: SSC std=%.3f p=%.4f\n", getcoef(fit_uw,"subjective_class")$std, getcoef(fit_uw,"subjective_class")$p))
if (!is.null(fit_1f)) { f1 <- fitMeasures(fit_1f, c("cfi","rmsea","srmr")); cat(sprintf("  1-factor CFA CFI=%.3f RMSEA=%.3f SRMR=%.3f\n", f1["cfi"], f1["rmsea"], f1["srmr"])) }

# ============================================================
# 7. OBJECTIVE-SES COMPOSITE (fair head-to-head)
# ============================================================
cat("\n[Objective SES composite]\n")
fit_m4a <- sem(paste(MEAS, "\nhappiness ~ obj_ses +", CTRL_NOEMP), data=df_clean, estimator="MLR")
fit_m4b <- sem(paste(MEAS, "\nhappiness ~ obj_ses + subjective_class +", CTRL_NOEMP), data=df_clean, estimator="MLR")
fit_m4c <- tryCatch(sem(paste0(MEAS, "
  method =~ a*cantril + a*satisfaction + a*yesterday_happy + a*eudaimonia +
            a*subjective_class + a*health_mental + a*health_physical + a*autonomy +
            a*yesterday_depress_r + a*financial_worry_r + a*covid_economy
  happiness ~~ 0*method
  method ~~ 1*method
  happiness ~ obj_ses + subjective_class + ", CTRL_NOEMP), data=df_clean, estimator="MLR"), error=function(e) NULL)
cat(sprintf("  M4a ObjSES std=%.3f | M4b ObjSES std=%.3f SSC std=%.3f\n",
            getcoef(fit_m4a,"obj_ses")$std, getcoef(fit_m4b,"obj_ses")$std, getcoef(fit_m4b,"subjective_class")$std))

# ============================================================
# 8. CONSTRUCT-OVERLAP CHECK: exclude psychological-distress predictors
#    (drop mental health, autonomy, depression) to test R2 inflation
# ============================================================
cat("\n[Psych-resource-excluded models]\n")
CTRL_NP <- "health_physical + female + age_c + age_sq_c + married + employed + financial_worry_r + covid_economy"
fit_np1 <- sem(paste(MEAS, "\nhappiness ~ log_income +", CTRL_NP), data=df_clean, estimator="MLR")
fit_np2 <- sem(paste(MEAS, "\nhappiness ~ subjective_class +", CTRL_NP), data=df_clean, estimator="MLR")
fit_np3 <- sem(paste(MEAS, "\nhappiness ~ log_income + subjective_class +", CTRL_NP), data=df_clean, estimator="MLR")
cat(sprintf("  no-psych: M1 R2=%.4f M2 R2=%.4f M3 R2=%.4f | SSC(M3) std=%.3f income(M3) std=%.3f p=%.4f\n",
            r2of(fit_np1), r2of(fit_np2), r2of(fit_np3),
            getcoef(fit_np3,"subjective_class")$std, getcoef(fit_np3,"log_income")$std, getcoef(fit_np3,"log_income")$p))

# ============================================================
# 9. BLOCK-WISE INCREMENTAL R2
# ============================================================
cat("\n[Hierarchical R2 blocks]\n")
blkA <- "female + age_c + age_sq_c + married + employed"
blkB <- paste(blkA, "+ log_income + subjective_class")
blkC <- paste(blkB, "+ health_mental + health_physical + autonomy + yesterday_depress_r + financial_worry_r + covid_economy")
rA <- r2of(sem(paste(MEAS,"\nhappiness ~",blkA), data=df_clean, estimator="MLR"))
rB <- r2of(sem(paste(MEAS,"\nhappiness ~",blkB), data=df_clean, estimator="MLR"))
rC <- r2of(sem(paste(MEAS,"\nhappiness ~",blkC), data=df_clean, estimator="MLR"))
cat(sprintf("  R2: demog=%.4f  +econ=%.4f (Δ=%.4f)  +psych=%.4f (Δ=%.4f)\n", rA, rB, rB-rA, rC, rC-rB))

# ============================================================
# 10. FIML vs listwise (missing data)
# ============================================================
cat("\n[FIML]\n")
df_fiml <- df2   # keep NA rows; recenter age on available
df_fiml$age_c    <- df_fiml$age - mean(df_fiml$age, na.rm=TRUE)
df_fiml$age_sq_c <- df_fiml$age_c^2
fit_fiml <- tryCatch(sem(paste(MEAS, "\nhappiness ~ log_income + subjective_class +", CTRL),
                         data=df_fiml, estimator="MLR", missing="fiml"), error=function(e) NULL)
n_fiml <- if (!is.null(fit_fiml)) lavInspect(fit_fiml,"nobs") else NA
if (!is.null(fit_fiml)) cat(sprintf("  FIML N=%d | SSC std=%.3f income std=%.3f p=%.4f\n",
       n_fiml, getcoef(fit_fiml,"subjective_class")$std, getcoef(fit_fiml,"log_income")$std, getcoef(fit_fiml,"log_income")$p))

# ============================================================
# 11. WEIGHTED SEM (survey sampling weights)
# ============================================================
cat("\n[Weighted SEM]\n")
fit_w <- tryCatch(sem(paste(MEAS, "\nhappiness ~ log_income + subjective_class +", CTRL),
                      data=df_clean, estimator="MLR", sampling.weights="weight"), error=function(e){cat(" wfail:",conditionMessage(e),"\n");NULL})
if (!is.null(fit_w)) cat(sprintf("  weighted: SSC std=%.3f p=%.4f income std=%.3f p=%.4f\n",
       getcoef(fit_w,"subjective_class")$std, getcoef(fit_w,"subjective_class")$p,
       getcoef(fit_w,"log_income")$std, getcoef(fit_w,"log_income")$p))

# ============================================================
# 12. MEDIATION income -> SSC -> happiness
# ============================================================
cat("\n[Mediation income->SSC->happiness]\n")
med_mod <- paste0(MEAS, "
  subjective_class ~ aa*log_income + health_mental + health_physical + autonomy + female + age_c + age_sq_c + married + employed + yesterday_depress_r + financial_worry_r + covid_economy
  happiness ~ cp*log_income + bb*subjective_class + health_mental + health_physical + autonomy + female + age_c + age_sq_c + married + employed + yesterday_depress_r + financial_worry_r + covid_economy
  indirect := aa*bb
  total := cp + aa*bb")
fit_med <- tryCatch(sem(med_mod, data=df_clean, estimator="MLR"), error=function(e) NULL)
if (!is.null(fit_med)) {
  pm <- parameterEstimates(fit_med, standardized=TRUE)
  ind <- pm[pm$label=="indirect",]; tot <- pm[pm$label=="total",]; dir <- pm[pm$label=="cp",]
  cat(sprintf("  a*b indirect=%.4f p=%.4f | direct(cp)=%.4f p=%.4f | total=%.4f\n",
              ind$est, ind$pvalue, dir$est, dir$pvalue, tot$est))
}

# ============================================================
# 13. PER-FACET robustness (each happiness indicator separately)
# ============================================================
cat("\n[Per-facet OLS]\n")
facet <- list()
for (y in c("cantril","satisfaction","yesterday_happy","eudaimonia")) {
  f <- as.formula(paste(y, "~ log_income + subjective_class +", CTRL))
  m <- lm(f, data=df_clean); cc <- summary(m)$coefficients
  facet[[y]] <- list(
    ssc_b=round(cc["subjective_class","Estimate"],4), ssc_p=round(cc["subjective_class","Pr(>|t|)"],6),
    inc_b=round(cc["log_income","Estimate"],4), inc_p=round(cc["log_income","Pr(>|t|)"],6),
    r2=round(summary(m)$r.squared,4))
  cat(sprintf("  %s: SSC b=%.3f p=%.4f | income b=%.3f p=%.4f\n",
              y, facet[[y]]$ssc_b, facet[[y]]$ssc_p, facet[[y]]$inc_b, facet[[y]]$inc_p))
}

# ============================================================
# 14. DESCRIPTIVES (recompute with corrected financial worry)
# ============================================================
desc_vars <- c("cantril","satisfaction","yesterday_happy","eudaimonia","subjective_class",
               "log_income","health_mental","health_physical","autonomy","covid_economy",
               "yesterday_depress_r","financial_worry_r","female","age","married","employed")
dd <- psych::describe(df_clean[, desc_vars])
desc <- lapply(desc_vars, function(v) list(mean=round(dd[v,"mean"],3), sd=round(dd[v,"sd"],3),
                                           min=round(dd[v,"min"],3), max=round(dd[v,"max"],3)))
names(desc) <- desc_vars

# ============================================================
# SAVE
# ============================================================
out <- list(
  note = "Revision re-analysis. financial_worry corrected B3->C8. GFI dropped.",
  software = list(R=R.version.string, lavaan=as.character(packageVersion("lavaan"))),
  data = list(n_total=nrow(df2), n_analysis=N, age_min=min(df_clean$age), age_max=max(df_clean$age)),
  measurement = list(
    note="Predictors mental health(C7), physical health(C6), autonomy(C5), depression(C3), financial worry(C8), COVID economy(C10), SSC(C13) are all SINGLE survey items; SSC=1-5, others 0-10.",
    happiness_alpha=round(alpha_h,4), happiness_omega=round(omega_h,4), happiness_ave=round(ave_h,4),
    cfa_fit=list(CFI=round(fi_cfa["cfi"],3), TLI=round(fi_cfa["tli"],3),
                 RMSEA=round(fi_cfa["rmsea"],3), RMSEA_lo=round(fi_cfa["rmsea.ci.lower"],3),
                 RMSEA_hi=round(fi_cfa["rmsea.ci.upper"],3), SRMR=round(fi_cfa["srmr"],3),
                 chisq=round(fi_cfa["chisq"],3), df=as.integer(fi_cfa["df"]), p=as.numeric(fi_cfa["pvalue"])),
    std_loadings=round(lam,3)),
  obj_ses = list(components=c("z_income","z_edu","z_employed"), alpha=round(alpha_obj,4)),
  icc = as.list(icc),
  vif = as.list(vif),
  correlations = cormat,
  competing = list(
    M1 = list(r2=r2of(fit_m1), fit=fitlist(fit_m1), income=getcoef(fit_m1,"log_income")),
    M2 = list(r2=r2of(fit_m2), fit=fitlist(fit_m2), subjective_class=getcoef(fit_m2,"subjective_class")),
    M3 = list(r2=r2of(fit_m3), fit=fitlist(fit_m3), income=getcoef(fit_m3,"log_income"),
              subjective_class=getcoef(fit_m3,"subjective_class")),
    delta_r2 = list(M1_to_M3=round(r2of(fit_m3)-r2of(fit_m1),4), M2_to_M3=round(r2of(fit_m3)-r2of(fit_m2),4)),
    dAIC_M1_M2 = round(fitlist(fit_m2)$AIC - fitlist(fit_m1)$AIC,1)),
  full_model_coefs = m3_coefs,
  cmv = list(
    ulmc_narrow = if(!is.null(fit_un)) list(ssc=getcoef(fit_un,"subjective_class"), income=getcoef(fit_un,"log_income")) else NA,
    ulmc_wide   = if(!is.null(fit_uw)) list(ssc=getcoef(fit_uw,"subjective_class"), income=getcoef(fit_uw,"log_income")) else NA,
    one_factor_cfa = if(!is.null(fit_1f)) { f1<-fitMeasures(fit_1f,c("cfi","tli","rmsea","srmr")); list(CFI=round(f1["cfi"],3),TLI=round(f1["tli"],3),RMSEA=round(f1["rmsea"],3),SRMR=round(f1["srmr"],3)) } else NA),
  objective_ses_models = list(
    M4a = list(r2=r2of(fit_m4a), obj_ses=getcoef(fit_m4a,"obj_ses")),
    M4b = list(r2=r2of(fit_m4b), obj_ses=getcoef(fit_m4b,"obj_ses"), subjective_class=getcoef(fit_m4b,"subjective_class")),
    M4c = if(!is.null(fit_m4c)) list(r2=r2of(fit_m4c), obj_ses=getcoef(fit_m4c,"obj_ses"), subjective_class=getcoef(fit_m4c,"subjective_class")) else NA),
  overlap_excluded = list(
    M1=list(r2=r2of(fit_np1), income=getcoef(fit_np1,"log_income")),
    M2=list(r2=r2of(fit_np2), subjective_class=getcoef(fit_np2,"subjective_class")),
    M3=list(r2=r2of(fit_np3), income=getcoef(fit_np3,"log_income"), subjective_class=getcoef(fit_np3,"subjective_class"))),
  hierarchical_r2 = list(demographics=rA, plus_economic=rB, plus_psych=rC,
                         d_econ=round(rB-rA,4), d_psych=round(rC-rB,4)),
  fiml = if(!is.null(fit_fiml)) list(n=n_fiml, subjective_class=getcoef(fit_fiml,"subjective_class"), income=getcoef(fit_fiml,"log_income")) else NA,
  weighted = if(!is.null(fit_w)) list(subjective_class=getcoef(fit_w,"subjective_class"), income=getcoef(fit_w,"log_income")) else NA,
  mediation = if(!is.null(fit_med)) { pm<-parameterEstimates(fit_med,standardized=TRUE); list(
      a_income_to_ssc=round(pm[pm$label=="aa","est"],4), a_p=round(pm[pm$label=="aa","pvalue"],6),
      b_ssc_to_happy=round(pm[pm$label=="bb","est"],4), b_p=round(pm[pm$label=="bb","pvalue"],6),
      direct_cp=round(pm[pm$label=="cp","est"],4), direct_p=round(pm[pm$label=="cp","pvalue"],6),
      indirect=round(pm[pm$label=="indirect","est"],4), indirect_p=round(pm[pm$label=="indirect","pvalue"],6),
      total=round(pm[pm$label=="total","est"],4)) } else NA,
  per_facet = facet,
  descriptives = desc
)
write_json(out, file.path(results_dir, "r_revision_results.json"), pretty=TRUE, auto_unbox=TRUE, digits=6)
cat("\nSaved results/r_revision_results.json\n=== DONE ===\n")
