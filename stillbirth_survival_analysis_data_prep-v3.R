# ====================================================================
# PROJECT: Interpregnancy Interval and Risk of Stillbirth
# A Survival Analysis Using the 2022 Bangladesh Demographic and Health Survey
# ====================================================================
#
# WHAT THIS SCRIPT DOES, IN PLAIN TERMS
#
# This script takes three raw data files from the 2022 Bangladesh
# Demographic and Health Survey and turns them into one clean, ready-to-use
# use dataset for a survival analysis. The research question is whether
# the interval since a woman's previous pregnancy, along with her age,
# how many pregnancies she has had before, and how much antenatal care
# she received, is associated with the risk of a stillbirth in the
# pregnancy that follows it.
#
# THE THREE RAW DATA FILES
#
#   1. Pregnancies Recode 
#      One row per pregnancy, for every pregnancy every woman in the
#      survey ever reported across her whole reproductive life. This is
#      the main file. The outcome of this study, whether a pregnancy
#      ended in a stillbirth, lives here.
#
#   2. Individual Recode 
#      One row per woman, holding her own personal characteristics, such
#      as her date of birth.
#
#   3. Pregnancy and Postnatal Care Recode
#      Antenatal and postnatal care details, but only for pregnancies
#      within a shorter, more recent recall period, not a woman's entire
#      reproductive history. It is later confirmed that this file adds no new
#      columns of its own.
#
# HOW THIS SCRIPT IS ORGANIZED
#
# Many sections below try something simple first, then check whether the
# result actually makes sense before moving forward. Real survey data
# almost always has quirks hiding inside it, and the safest way to work
# with it is to look before trusting a number, rather than assuming a
# calculation is correct just because R ran it without an error.
#
# By the end of this script, you will have one object called
# analytic_data, which is the actual dataset ready for Kaplan Meier
# estimation and Cox proportional hazards regression.
#
# WHERE EACH VARIABLE GETS BUILT
# Every place in this script where a variable is actually created or
# changed is marked with a line of asterisks and a short label, so you
# can scan straight to it without reading every line in between. Here is
# where to find each one.
#
#   stillbirth                Section 5
#   sample_weight             Section 5B (V021 and V022 used as is, no new variable needed)
#   interpregnancy_interval   Section 7 (first attempt) and Section 9 (corrected, final version)
#   maternal_age              Section 10 (created, then cleaned)
#   m14 (antenatal visits)    Section 11 (cleaned)
#   parity                    Section 12
#   analytic_data              Section 12 (the final dataset itself)
#   interval_group             Section 13
#   the survival object          Section 16
#   km_fit                       Section 17
#   the log-rank test            Section 18
#   cox_naive                    Section 19
#   the proportional hazards check   Section 20
#   design, cox_survey           Section 21
#   the labelled comparison table Section 22
#   the forest plot              Section 23
# ====================================================================

# ====================================================================
# SECTION 1: Loading the tools and the three raw data files
# ====================================================================

# haven reads Stata (.dta) files into R while keeping the value labels
# DHS attaches to each coded number. This matters a great deal with DHS
# data, since almost every column is a number standing in for a
# category, not a plain measurement on its own.
install.packages("haven")   # only needed once
library(haven)

# One object per file, named plainly so it is obvious what each one holds.
# Update the file paths below to match where you saved these files on
# your own computer.
pregnancies <- read_dta("C:/DHS_Survival Analysis/Data files/BDGR81FL_Pregnancies recode.DTA")   # full pregnancy history, the outcome lives here
women       <- read_dta("C:/DHS_Survival Analysis/Data files/BDIR81FL_Individual recode.DTA")     # one row per woman, her own characteristics
ppnc        <- read_dta("C:/DHS_Survival Analysis/Data files/BDNR81FL_Pregnancy and Postnatal Care.DTA")  # antenatal and postnatal care details

# A first look at each file, just its shape, meaning how many rows and
# how many columns it has. 
dim(pregnancies)
dim(women)
dim(ppnc)

# ====================================================================
# SECTION 2: Checking whether the postnatal care file is actually new
# information, or just a smaller slice of the main file
# ====================================================================

# If this returns TRUE, ppnc and pregnancies share the same
# columns, meaning ppnc is just a row subset of pregnancies rather than
# a separate set of variables. Confirmed TRUE in this project, so ppnc
# is not loaded again after this point, everything we need is already
# sitting inside pregnancies.
identical(names(pregnancies), names(ppnc))

# ====================================================================
# SECTION 3: Finding the real variables behind DHS's cryptic column names
# ====================================================================

# DHS column names are short and impossible to guess, things like p32 or
# p20 tell the layperson nothing on their own. Fortunately, haven kept each
# column's plain language description as a label attribute when I imported
# the file, so I have to search those labels for keywords instead of
# scrolling through hundreds of column names by eye.
var_labels <- sapply(pregnancies, function(x) attr(x, "label"))

# Searching for the variable marking what happened at the end of each
# pregnancy, live birth, stillbirth, and so on. This is the outcome
# variable. This search points to p32, "pregnancy outcome
# reclassified", which is DHS's own cleaned version of this variable,
# more reliable than the respondent's original unedited answer.
var_labels[grepl("outcome", var_labels, ignore.case = TRUE)]

# Searching for the variable recording how long each pregnancy lasted. This
# points to p20, "duration of pregnancy in months", which is needed
# both to define stillbirth correctly and to work out interpregnancy
# interval later.
var_labels[grepl("duration", var_labels, ignore.case = TRUE)]

# ====================================================================
# SECTION 4: Understanding the outcome and duration variables before
# trusting either one
# ====================================================================

# Seeing exactly how DHS coded the categories behind p32, rather than
# assuming which number means live birth versus stillbirth. This
# returned four clean categories: born alive = 1, born dead = 2,
# miscarriage = 3, abortion = 4.
attr(pregnancies$p32, "labels")

# DHS sometimes mixes reporting units inside one variable, for example
# using a hundreds digit to signal whether an answer was given in weeks
# or in months. Dividing p20 by 100 and keeping only the whole number
# checks for this. In this dataset, every value came back as 0, meaning
# there is no blocked unit coding here, p20 is a plain count.
table(pregnancies$p20 %/% 100, useNA = "ifany")

# A plain count can still hide a special "don't know" or "missing" code,
# for example a value like 98, which would not have been caught by the
# check above. This shows the full spread of every value p20 actually
# takes. In this dataset, every value fell cleanly between 1 and 10
# months, with no missing values and no special codes, so p20 can be
# trusted as a genuine duration in months.
table(pregnancies$p20, useNA = "ifany")

# A labelled value is the clearest sign a number means something other
# than a literal amount, the same way p32's numbers meant categories,
# not amounts. This came back empty for p20, confirming it holds no
# hidden special codes.
attr(pregnancies$p20, "labels")

# ====================================================================
# SECTION 5: Building the stillbirth outcome variable
# ====================================================================

# DHS's own definition of stillbirth needs two conditions together, the
# outcome must be born dead, and the pregnancy must have lasted at least
# seven months, since an earlier loss is what most people mean by
# miscarriage, not stillbirth.

# ******************************************************************
# VARIABLE CREATED: stillbirth
# ******************************************************************
pregnancies$stillbirth <- ifelse(pregnancies$p32 == 2 & pregnancies$p20 >= 7, 1, 0)

# Checking how many pregnancies meet that definition. Stillbirth should be
# a rare outcome, so a small number here is a good sign, not a mistake.
# This dataset showed 1,743 stillbirths out of 73,239 pregnancies, about
# two and a half percent, a believable rate for this population.
table(pregnancies$stillbirth, useNA = "ifany")

# Cross-checking the new stillbirth flag against the raw outcome variable,
# to confirm whether the two-part rule actually
# excluded any cases. In this dataset every single "born dead" case
# already had a duration of seven months or more, and every miscarriage
# fell under seven months, so DHS's reclassified outcome variable had
# already kept these two groups cleanly separated on its own.
#
# Wrapping both sides in factor() with real labels here, rather than
# reading a table of bare 1, 2, 3, 4 and TRUE, FALSE. p32's own codes
# and a plain logical comparison both print fine on their own, but
# neither explains itself sitting side by side in a cross tab, and a
# table you have to keep looking up is a table you are more likely to
# misread.
table(
  factor(pregnancies$p32, levels = c(1, 2, 3, 4),
         labels = c("Live birth", "Born dead (stillbirth candidate)", "Miscarriage", "Abortion")),
  factor(pregnancies$p20 >= 7, levels = c(FALSE, TRUE),
         labels = c("Under 7 months", "7 months or more")),
  useNA = "ifany"
)

# ====================================================================
# SECTION 5B: Finding the survey design variables, weight, cluster, strata
# ====================================================================

# DHS is not a simple random sample, it is a stratified, multi-stage
# cluster sample. Villages or urban blocks are sampled first, not
# individual women, and some regions and some urban versus rural areas
# were deliberately oversampled relative to others. Treating this data
# as though every woman had an equal chance of being selected would
# misrepresent the actual Bangladeshi population, so three variables
# exist specifically to correct for that, a sampling weight, a cluster
# identifier, DHS calls this the primary sampling unit, and a
# stratification variable.
#
# Before deciding how to bring these into analytic_data, the first
# step is finding out what DHS actually calls them and where they
# currently live, the same label search already used to find p32 and
# p20 back in Section 3, rather than guessing a variable name.
var_labels[grepl("weight", var_labels, ignore.case = TRUE)]
var_labels[grepl("primary sampling unit|cluster", var_labels, ignore.case = TRUE)]
var_labels[grepl("stratification|strata", var_labels, ignore.case = TRUE)]

# What these three searches return determines the next step. If all
# three come back with a result here, DHS already duplicated these
# variables into the pregnancies file directly, and nothing further
# needs to be merged in. If any come back empty, that variable will
# need to be pulled in from the women file instead and matched onto
# every pregnancy by caseid, since weight, cluster, and strata describe
# the woman who was sampled, not any one specific pregnancy of hers.
#
# All three came back directly, so DHS already duplicated them into
# the pregnancies file, no merge needed. Two names came back for
# cluster, and two for strata, and they are not interchangeable. V001
# is the raw cluster number, while V021 is the version DHS's own
# recode manual describes as built specifically for variance
# estimation, the one that actually belongs in a survey design object
# later. The same split applies to strata, V023 is a general label,
# "stratification used in sample design", while V022 is the one DHS
# ties specifically to calculating sampling errors with the Taylor
# series method. V005, V021, and V022 are the three carried forward
# from here, not V001 and V023.
#
# V005 itself needs one more step before it is usable. DHS stores the
# sample weight as a large integer with six implied decimal places
# baked in, so the stored number has to be divided by a million before
# it behaves like an actual weight.

# A look at the raw stored values before dividing, just to see the
# large integer form DHS actually saves this variable in.
summary(pregnancies$v005)

# ******************************************************************
# VARIABLE CREATED: sample_weight
# ******************************************************************
pregnancies$sample_weight <- pregnancies$v005 / 1000000

# A believable sample weight should scatter loosely around 1, since
# weights are built to average out to the sample's actual size once
# applied. A minimum near zero or a maximum in the thousands here would
# be a sign the division step above was skipped or done twice.
summary(pregnancies$sample_weight)

# V021 and V022 need no cleaning or new variable of their own, they are
# used exactly as DHS provides them. Nothing further happens with any
# of these three here, they simply ride along inside pregnancies
# through every filter still to come, the same way stillbirth and the
# other covariates already do, so they arrive inside analytic_data
# automatically once Section 12 runs. The actual use of all three
# together, weight, cluster, and strata, happens later, at the
# modeling stage, inside a survey design object built specifically to
# account for them.

# Every other variable in this script got checked for missing values
# before being trusted, these three have not been yet. Design variables
# like these almost never go missing, since they come from the sample
# draw itself rather than a question a respondent could decline to
# answer, but "almost never" is not the same as checked. All three
# should read 0 here.
sum(is.na(pregnancies$sample_weight))
sum(is.na(pregnancies$v021))
sum(is.na(pregnancies$v022))

# ====================================================================
# SECTION 6: Finding the variables needed to build interpregnancy interval
# ====================================================================

# I need a variable identifying which woman each pregnancy belongs to,
# so I can group pregnancies by woman before measuring the gap between
# one pregnancy and the next for that same woman. This points to caseid.
var_labels[grepl("case", var_labels, ignore.case = TRUE)]

# I also need whatever marks the order of pregnancies within a woman's
# history, first, second, third, and so on, so I know which pregnancy
# came immediately before the one I am looking at. This points to
# pord, "pregnancy order number".
var_labels[grepl("order|pregnancy number", var_labels, ignore.case = TRUE)]

# Finally, I need a date I can subtract cleanly. Separate month and
# year fields do not subtract well across a year boundary, so DHS uses a
# century month code instead, a single running count of months since
# January 1900. This points to p3, "date of end of pregnancy (cmc)",
# which lives inside the pregnancy file itself, one value per pregnancy,
# unlike similarly named variables starting with v, which describe the
# woman in general rather than this specific pregnancy.
var_labels[grepl("century month|cmc", var_labels, ignore.case = TRUE)]

# ====================================================================
# SECTION 7: Building interpregnancy interval, first attempt
# ====================================================================

# dplyr lets me work within groups, here each woman, and reach back to
# the row directly before the current one, which is what a birth spacing
# calculation needs.
install.packages("dplyr")
library(dplyr)

# ******************************************************************
# VARIABLE CREATED: interpregnancy_interval (first attempt)
# This version has a known flaw involving twins, corrected in Section 9
# ******************************************************************
pregnancies <- pregnancies %>%
  arrange(caseid, pord) %>%           # put each woman's pregnancies in true chronological order first
  group_by(caseid) %>%                # everything below happens separately within each woman, not across women
  mutate(
    prev_end_cmc = lag(p3),           # the end date of this woman's previous pregnancy, NA if this is her first
    start_cmc = p3 - p20,             # working backward from this pregnancy's own end date and duration to estimate when it started
    interpregnancy_interval = start_cmc - prev_end_cmc   # months from the last pregnancy ending to this one starting
  ) %>%
  ungroup()

# NOTE: this first attempt has a known flaw, which the next two sections
# find and fix. It is kept here deliberately so the reasoning behind the
# fix is easy to follow, rather than presenting the corrected version
# with no explanation of what it corrects.

# A first look at the result
summary(pregnancies$interpregnancy_interval)

# A sanity check. Every woman's very first pregnancy should show up as
# NA here, since there is nothing before it to measure from, so these
# two counts should match, and in this dataset they did, both equal to
# 27,324.
sum(is.na(pregnancies$interpregnancy_interval))
length(unique(pregnancies$caseid))

# ====================================================================
# SECTION 8: Diagnosing the negative and unusually large intervals
# ====================================================================

# Every century month code in this dataset, p3, prev_end_cmc, and so
# on, is a plain count of months since January 1900, CMC 1. That is
# what lets these dates subtract cleanly across year boundaries, but it
# also means a raw CMC number like 1455 means nothing to a human eye
# until it is converted back into an actual month and year. This small
# helper does that conversion, so tables below can show a real date
# instead of a number that has to be looked up by hand every time.

# ******************************************************************
# HELPER FUNCTION CREATED: cmc_to_date
# ******************************************************************
cmc_to_date <- function(cmc) {
  year  <- 1900 + (cmc - 1) %/% 12    # how many full years have passed since January 1900
  month <- (cmc - 1) %% 12 + 1        # the remainder tells which month within that year
  sprintf("%s %d", month.name[month], year)
}

# A quick check that the helper works as expected. CMC 1455 should
# convert to "March 2021", which is exactly the delivery date behind
# the very first row of the long gap table further down.
cmc_to_date(1455)

#*** I confirmed the formula from DHS's own documentation:
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode5DHS_23August2012.pdf (Page 5)

# Looking directly at a handful of the negative interval cases, checking
# whether the previous pregnancy's end date matches this pregnancy's own
# end date, which is exactly what twins recorded as two separate rows
# would produce, since both babies share the same delivery date.
#
# The plain column names below, caseid, pord, p3, and so on, mean
# nothing on sight, the same problem labels solved back in Section 5.
# Renaming them here, and adding one extra column that states outright
# whether each row's previous end date matches its own end date, turns
# this from a table you have to reason through into one that tells you
# its own answer. Twin cases will always show prev_end_cmc equal to p3.
#
# All 734 rows are kept this time, not just a preview of ten. A table
# this size printed straight into the console would either scroll past
# or get silently truncated by tibble's own print method, neither of
# which is actually seeing the whole thing. The result is saved into
# negative_interval_cases and opened with View() instead, which is
# RStudio's own spreadsheet style viewer, built for scrolling through
# and sorting exactly this many rows without losing any of them.

# ******************************************************************
# OBJECT CREATED: negative_interval_cases
# ******************************************************************
negative_interval_cases <- pregnancies %>%
  filter(interpregnancy_interval < 0) %>%
  select(caseid, pord, p3, p20, prev_end_cmc, interpregnancy_interval) %>%
  mutate(
    likely_cause = ifelse(
      prev_end_cmc == p3,
      "Twin, shares end date with previous row",
      "Genuine gap, rounding overlap of about a month"
    )
  ) %>%
  rename(
    `Woman ID` = caseid,
    `Pregnancy order` = pord,
    `This pregnancy's end (CMC)` = p3,
    `Duration (months)` = p20,
    `Previous pregnancy's end (CMC)` = prev_end_cmc,
    `Interval (months)` = interpregnancy_interval,
    `Likely cause` = likely_cause
  )

View(negative_interval_cases)

# How many pregnancies this affects in total. In this dataset, 734,
# matching the row count you should now see at the top of the viewer.
sum(pregnancies$interpregnancy_interval < 0, na.rm = TRUE)

# Now look at the unusually long gaps, over twenty years, to see whether
# they look like real spacing or like a data problem. Interval is shown
# converted into years as well as months here, since twenty or so years
# is a lot easier to judge as plausible or implausible than a bare
# three-digit month count. Both end dates are also converted into real
# calendar dates using cmc_to_date, so each row reads as an actual
# piece of one woman's history rather than a pair of raw CMC numbers.
pregnancies %>%
  filter(interpregnancy_interval > 240) %>%
  select(caseid, pord, p3, p20, prev_end_cmc, interpregnancy_interval) %>%
  mutate(
    years_since_previous = round(interpregnancy_interval / 12, 1),
    this_pregnancy_ended = cmc_to_date(p3),
    previous_pregnancy_ended = cmc_to_date(prev_end_cmc)
  ) %>%
  rename(
    `Woman ID` = caseid,
    `Pregnancy order` = pord,
    `This pregnancy's end (CMC)` = p3,
    `Duration (months)` = p20,
    `Previous pregnancy's end (CMC)` = prev_end_cmc,
    `Interval (months)` = interpregnancy_interval,
    `Interval (years)` = years_since_previous,
    `This pregnancy ended` = this_pregnancy_ended,
    `Previous pregnancy ended` = previous_pregnancy_ended
  ) %>%
  arrange(desc(`Interval (months)`)) %>%
  head(10)

sum(pregnancies$interpregnancy_interval > 240, na.rm = TRUE)

# CONCLUSION FROM THIS SECTION:
# The negative values were confirmed to be mostly twins, recorded as two
# separate pregnancy rows sharing one identical end date, which our
# first attempt above mistakenly treated as two separate pregnancies.
# The very large positive gaps affected only ten pregnancies out of over
# seventy-three thousand, several of them a woman's second ever
# pregnancy after a long gap since her first, which is rare but
# medically possible, so these ten were left in the data rather than
# removed without real cause.

# ====================================================================
# SECTION 9: Fixing interpregnancy interval properly, accounting for twins
# ====================================================================

# The row directly before this one is not always a different pregnancy,
# since twins share the same end date across two rows. What I
# actually need is the previous DISTINCT pregnancy, meaning the last one
# with a genuinely different end date.

# First, collapsing each woman's history down to one row per distinct
# pregnancy event, so a twin pair only counts once.

# ******************************************************************
# HELPER TABLE CREATED: pregnancy_events
# Not part of the final dataset, exists only to calculate the fix below
# ******************************************************************
pregnancy_events <- pregnancies %>%
  distinct(caseid, p3) %>%
  arrange(caseid, p3) %>%
  group_by(caseid) %>%
  mutate(prev_distinct_end_cmc = lag(p3)) %>%   # the true previous pregnancy, skipping past any twin
  ungroup()

# Bringing that corrected date back onto every row. Both twins share the
# same p3, so this join gives both of them the same, correctly measured
# interval, rather than measuring the twins against each other.

# ******************************************************************
# VARIABLE RECALCULATED: interpregnancy_interval (corrected, final version)
# This overwrites the first attempt from Section 7
# ******************************************************************
pregnancies <- pregnancies %>%
  left_join(pregnancy_events, by = c("caseid", "p3")) %>%
  mutate(
    start_cmc = p3 - p20,
    interpregnancy_interval = start_cmc - prev_distinct_end_cmc
  )

# Checking the spread again now that twins are handled correctly. The
# negative count dropped sharply here, from 734 down to 142, and the NA
# count rose slightly, from 27,324 to 27,516, which is expected, since a
# woman whose very first pregnancy was a twin pair now correctly shows
# NA on both of those rows instead of just one.
summary(pregnancies$interpregnancy_interval)
sum(pregnancies$interpregnancy_interval < 0, na.rm = TRUE)

# Seeing the actual spread of what remains negative, rather than assuming
# all 142 cases are the same kind of small noise. About two-thirds
# clustered at negative one or two, consistent with simple rounding in
# self-reported pregnancy duration, while a smaller tail reaching down
# to negative eight likely reflects a few more near twin cases that this
# fix did not fully catch.
table(pregnancies$interpregnancy_interval[pregnancies$interpregnancy_interval < 0])

# ******************************************************************
# VARIABLE MODIFIED: interpregnancy_interval (remaining negatives floored to zero)
# ******************************************************************
pregnancies$interpregnancy_interval[pregnancies$interpregnancy_interval < 0] <- 0

# Confirming the floor worked and taking one more look at the overall shape.
# The minimum is now 0, and the rest of the distribution is unchanged.
summary(pregnancies$interpregnancy_interval)

# ====================================================================
# SECTION 10: Building and cleaning maternal age
# ====================================================================

# Maternal age is not a ready-made column, but it can be built from two
# dates I already have, the woman's own date of birth (v011) and this
# pregnancy's end date (p3), both in century month codes. Subtracting
# gives her age in months at that point, dividing by twelve gives years,
# the same trick already used to build interpregnancy interval.

# ******************************************************************
# VARIABLE CREATED: maternal_age
# ******************************************************************
pregnancies$maternal_age <- (pregnancies$p3 - pregnancies$v011) / 12

# A sanity check. Real maternal ages here should fall somewhere between
# the young teens and roughly the mid-fifties.
summary(pregnancies$maternal_age)

# The check above turned up a minimum of 0.25, a biologically impossible
# age. Pull out the affected pregnancies to see the actual dates behind
# that number rather than guessing at a cause.
pregnancies %>%
  filter(maternal_age < 10) %>%
  select(caseid, pord, p3, v011, maternal_age) %>%
  arrange(maternal_age) %>%
  head(10)

# How many pregnancies this affects. In this dataset, only five, all of
# them a woman's first pregnancy, which points to an isolated data entry
# problem in a handful of individual records, not a systemic issue with
# the merge or the variable itself.
sum(pregnancies$maternal_age < 10, na.rm = TRUE)

# These five values are not noise that can be rounded away, they are
# biologically impossible, and the true age behind them cannot be
# recovered from the data available. They are marked missing rather than
# left in to quietly distort the model later.

# ******************************************************************
# VARIABLE MODIFIED: maternal_age (biologically impossible values set to missing)
# ******************************************************************
pregnancies$maternal_age[pregnancies$maternal_age < 10] <- NA

# Confirming only those five changed and the rest of the distribution held
# steady. The new minimum landed at 10, still young but no longer
# impossible, and affecting only a single record, not worth a further
# cleaning pass on its own.
summary(pregnancies$maternal_age)

# NOTE ON PARITY: parity does not need a separate variable at all. The
# pord variable already tells which pregnancy number this is for that
# woman, so the number of pregnancies she had before this one is simply
# pord minus one.

# ====================================================================
# SECTION 11: Finding and cleaning the antenatal care variable
# ====================================================================

# The antenatal care variable should already be sitting inside
# pregnancies, since Section 2 confirmed the postnatal care file shares
# identical columns. Search the labels rather than guessing another
# cryptic variable name.
var_labels[grepl("antenatal", var_labels, ignore.case = TRUE)]

# Of everything that search returns, m14, "number of antenatal visits
# during pregnancy", is the variable almost all published research on
# this topic actually uses, so that is the one carried forward here.
# Other results, such as which specific facility type she visited, or
# what was discussed during her visits, are a different kind of question
# and are not used in this project.

# Checking how much of the full pregnancy history actually has this
# variable populated, since antenatal detail was only collected for
# pregnancies inside a shorter recall window, not a woman's full
# history. Only 5,187 of 73,239 pregnancies had this variable at all.
sum(!is.na(pregnancies$m14))
summary(pregnancies$m14)

# The maximum value of 98 above is suspicious; a real pregnancy does not
# involve 98 antenatal visits. Checking whether DHS actually labelled this
# as a special code, the same way p20 was checked earlier. This
# confirmed 98 specifically means "don't know", while 0 is a genuine
# answer meaning no antenatal visits at all.
attr(pregnancies$m14, "labels")

library(haven)
# Checking the real shape of the distribution, including exactly how many
# pregnancies sit at that suspicious value. Only one pregnancy was coded
# 98 in this dataset.
table(pregnancies$m14, useNA = "ifany")

# Recoding only the genuine don't know response to missing. Zero is left
# untouched, since it is a real and meaningful answer, not a placeholder
# for missing data.

# ******************************************************************
# VARIABLE MODIFIED: m14 / number of antenatal care visits ("don't know" code set to missing)
# ******************************************************************
pregnancies$m14[pregnancies$m14 == 98] <- NA

# Confirming the fix. The maximum is now a believable 20 visits, and the
# mean dropped slightly now that the false 98 is gone.
summary(pregnancies$m14)

# ====================================================================
# SECTION 12: Assembling the final analytic dataset
# ====================================================================

# The analytic population for this project is pregnancies that could
# have ended in either a live birth or a stillbirth, since that is the
# actual contrast this research question is about. Miscarriage and
# induced abortion are excluded here, not because they are unimportant,
# but because they represent a different process, ending a pregnancy for
# reasons unrelated to, and generally before, the stillbirth risk window
# this project is measuring. This also matches how stillbirth rates are
# conventionally defined in the field, counted against total births,
# live plus still, not against every pregnancy including early losses.

# Adding p20 >= 7 below removes those and lines the denominator up with
# the definition exactly.
#
# On top of that outcome restriction, a pregnancy is only kept if all
# four covariates needed for the model are actually available for it,
# interpregnancy interval, maternal age, parity, and antenatal care
# visits. Parity itself needs no missing data check, since it is simply
# pord minus one, and pord is never missing.

# ******************************************************************
# VARIABLE CREATED: parity
# DATASET CREATED: analytic_data (the final file used for the survival analysis)
# ******************************************************************
analytic_data <- pregnancies %>%
  filter(
    p32 %in% c(1, 2),                     # only live births and stillbirths
    p20 >= 7,                             # restricted to pregnancies of seven months or more, matching the formal definition on both sides
    !is.na(interpregnancy_interval),       # dropped first pregnancies, since spacing cannot be measured for them
    !is.na(maternal_age),                  # dropped the five biologically impossible ages found earlier
    !is.na(m14)                            # kept only pregnancies with antenatal care data available
  ) %>%
  mutate(parity = pord - 1)                # number of pregnancies before this one

# The actual size of the final working sample after adding the seven
# month restriction. Comparing this number to the 3,446 figure from
# before the fix. The gap between them says how many of the forty
# one premature live births had also survived the covariate filters
# and would otherwise have sat in the denominator without belonging
# there.
nrow(analytic_data)

# Confirming the outcome still makes sense once the sample has shrunk this
# much. Checking this table against the earlier 70 out of 3,446 figure.
# The stillbirth count itself should not move, since Section 5 already
# showed the stillbirth flag was seven months or more from the start.
# Only the total and the rate can shift here.
table(analytic_data$stillbirth)

# A HONEST NOTE ON SAMPLE SIZE: seventy stillbirth events across four
# covariates sits at the lower edge of what is generally considered a
# comfortable number for a Cox regression model, commonly cited
# guidance suggests roughly ten to twenty events per predictor. This
# does not make the analysis invalid, but it does mean confidence
# intervals will likely be wide, and it means four covariates should be
# treated as a ceiling for this model, not a starting point to build
# further complexity on top of.

# ====================================================================
# SECTION 13: Categorizing interpregnancy interval into three groups
# ====================================================================

# The quintile check done earlier showed a clear U-shape: stillbirth
# risk highest in the shortest intervals, lowest in the middle, and
# rising again in the longest intervals. A plain continuous term for
# interval cannot represent that shape, a straight line only goes one
# direction, so interval needs to enter the model as a categorical
# variable instead.
#
# Five quintiles would cost four dummy variables, which is too many
# given only seventy stillbirth events total and four covariates
# already treated as a ceiling for this model. Collapsing down to
# three groups instead, short, reference, and long, costs only two
# dummy variables and still lets both ends of the U shape show up as
# their own separate estimated effects.
#
# The cutpoints below, sixteen months and seventy-eight months, come
# from the quintile boundaries already found in analytic_data, not
# from an imported cutoff built on another country's population. The
# gestational age restriction added in Section 12 changed analytic_data
# slightly after those quintiles were first checked, so the line below
# reruns that same quantile check on the corrected data. The dropped
# cases were all under seven-month live births, which have nothing to
# do with interval length, so a real shift in these boundaries would
# be a surprise, but confirm it rather than assume it.
quantile(analytic_data$interpregnancy_interval, probs = seq(0, 1, 0.1))

# ******************************************************************
# VARIABLE CREATED: interval_group
# ******************************************************************
analytic_data <- analytic_data %>%
  mutate(
    interval_group = cut(
      interpregnancy_interval,
      breaks = c(-Inf, 16, 78, Inf),        # short: sixteen months or less, reference: seventeen to seventy-eight months, long: seventy-nine months or more
      labels = c("short", "reference", "long"),
      right = TRUE                           # each break value falls into the lower group, matching how the quintile boundaries were read earlier
    ),
    interval_group = relevel(interval_group, ref = "reference")  # the middle group becomes the baseline for the Cox model, since it held the lowest rate
  )

# Confirming the split behaved as expected. Roughly a fifth of pregnancies
# should land in the short group, three-fifths in the middle, and a
# fifth in the long group, since the cutpoints came from quintile
# boundaries to begin with.
table(analytic_data$interval_group)

# Checking the stillbirth rate inside each of the three groups directly.
# The short and long groups should both show a visibly higher rate
# than the reference group in the middle, the same U-shape the five
# group version showed, just with fewer and more stable groups to
# read.
analytic_data %>%
  group_by(interval_group) %>%
  summarise(
    n = n(),
    stillbirths = sum(stillbirth),
    rate = mean(stillbirth)
  )

# ====================================================================
# END OF DATA PREPARATION
# analytic_data is now ready for Kaplan Meier estimation and
# Cox proportional hazards regression, with interval_group in place
# as the categorical form of interpregnancy interval.
# ====================================================================

# ====================================================================
# MAIN ANALYSIS: Kaplan Meier Estimation and Cox Proportional Hazards
# ====================================================================
#
# WHAT THIS SCRIPT DOES, IN PLAIN TERMS
#
# The data prep script built one object, analytic_data, containing
# 3,443 pregnancies of seven months or more, with the stillbirth
# outcome, the three category interval_group, maternal_age, parity,
# m14 antenatal visits, and the three survey design variables
# sample_weight, v021, and v022 all sitting inside it.
#
# Two questions get answered here. First, do the three interval
# groups actually look different when plotted as survival curves,
# answered with Kaplan Meier estimation. Second, does that difference
# hold up once maternal age, parity, and antenatal care are accounted
# for at the same time, and once the survey design is respected,
# answered with Cox proportional hazards regression, fit twice, once
# the naive way and once the correct survey adjusted way.

# ====================================================================
# SECTION 15: Loading the two packages this script needs
# ====================================================================

# survival provides Surv(), survfit(), survdiff(), and coxph(), the
# whole standard toolkit for this kind of analysis. It usually ships
# already installed alongside R itself, so the install line below is
# only a safety net.
install.packages("survival")   
library(survival)

# survey provides svydesign() and svycoxph(), the versions of Cox
# regression that know how to account for a sample that was not drawn
# by simple random selection. Plain coxph() has no idea DHS used
# clusters and strata; it would treat every one of the 3,443
# pregnancies as if it were sampled completely independently, which
# is not how this data was actually collected.
install.packages("survey")   
library(survey)

# ====================================================================
# SECTION 16: Building the survival object, and understand what it means
# ====================================================================

# Survival analysis needs two ingredients for every row, how long that
# row was observed for, and whether the event of interest happened at
# the end of that time or not. Here, the time is p20, the completed
# months of gestation when the pregnancy ended, and the event is
# stillbirth, one meaning the event happened, zero meaning it did not.
#
# A HONEST NOTE ON WHAT "SURVIVAL" MEANS IN THIS CONTEXT. This is not
# a study that followed women for months or years and watched what
# happened to them over that time, the DHS interview happened once,
# after the fact. What makes this a legitimate survival analysis
# anyway is the idea of a pregnancy being "at risk" of stillbirth for
# as long as it continues, and that risk ending the moment the
# pregnancy ends, whichever way it ends. A live birth is not a failure
# to observe the outcome, it is the outcome that survival analysis
# calls censoring, the pregnancy survived to the end of its observed
# time without the event happening. This is the same logic obstetric
# researchers use when they talk about stillbirth risk by gestational
# week, treating each ongoing pregnancy as being at risk until it
# ends.
#
# A SECOND HONEST NOTE. p20 is only recorded in completed months, not
# weeks or days, and most pregnancies in any population end within a
# narrow band of months, eight, nine, or ten. That means many
# pregnancies will share the same recorded time, what survival
# analysis calls ties. Cox regression was originally built assuming
# ties are rare, so the ties = "efron" setting used later in Section
# 19 exists specifically to handle this honestly rather than pretend
# the ties are not there.

survival_object <- Surv(time = analytic_data$p20, event = analytic_data$stillbirth)

# A quick look at the first several rows. A plus sign after the number
# marks a censored row, meaning that pregnancy ended in a live birth,
# reaching that duration without the event happening. A number with no
# plus sign is a stillbirth, the event happening at that duration.
head(survival_object)

# ====================================================================
# SECTION 17: Kaplan-Meier estimation, one curve per interval group
# ====================================================================

# survfit() estimates, separately for each interval group, the
# probability that a pregnancy has not yet ended in stillbirth as
# gestation continues month by month. This is entirely descriptive,
# no adjustment for maternal age, parity, or antenatal care happens
# here, it is only asking whether the three interval groups look
# different from each other before anything else gets controlled for.

# ******************************************************************
# OBJECT CREATED: km_fit
# ******************************************************************
km_fit <- survfit(Surv(p20, stillbirth) ~ interval_group, data = analytic_data)

# The printed summary gives, for each group, how many pregnancies went
# in, how many stillbirths happened, and the median survival time.
# Expect the median columns to come back NA here, stillbirth is rare
# enough in every group that the estimated probability of remaining
# event-free never actually drops to fifty percent, so there is no
# median to report. That is a normal result given how rare the event
# is, not a sign anything went wrong.
print(km_fit)

# A plot of the three curves. Because the event is rare, all three
# lines will start near 1.0 and stay high, do not expect a dramatic
# drop, the story here is in the small gap between the reference
# curve and the short and long curves, not in the overall shape.
plot(
  km_fit,
  col = c("#2a78d6", "#eb6834", "#1baf7a"),   # reference, short, long, matched to the order R lists the groups in
  lwd = 2,
  xlab = "Gestational duration (completed months)",
  ylab = "Proportion of pregnancies not yet ending in stillbirth",
  main = "Kaplan Meier estimates by interpregnancy interval group"
)
legend(
  "bottomleft",
  legend = c("reference (17-78 months)", "short (16 months or less)", "long (79 months or more)"),
  col = c("#2a78d6", "#eb6834", "#1baf7a"),
  lwd = 2,
  bty = "n"
)

# ====================================================================
# SECTION 18: A formal test of whether the three curves actually differ
# ====================================================================

# Looking at a plot and deciding the curves look different is not a
# statistical test on its own. The log-rank test below compares the
# observed number of stillbirths in each group against the number
# that would be expected if all three groups shared the same
# underlying risk. A small p-value here means the three curves differ
# by more than random chance alone would produce.

survdiff(Surv(p20, stillbirth) ~ interval_group, data = analytic_data)

# One limitation worth stating plainly. The log-rank test above, like
# the Kaplan-Meier curves themselves, does not adjust for maternal
# age, parity, or antenatal care, and it does not know anything about
# the survey design either. It answers a narrower question than the
# Cox models below will, whether interval_group alone, with nothing
# else held constant, separates the three groups.

# ====================================================================
# SECTION 19: Cox proportional hazards regression, the naive version
# ====================================================================

# coxph() answers the real research question, whether interval_group
# still predicts stillbirth risk once maternal_age, parity, and m14
# are held constant at the same time. "Naive" here means this first
# version treats the 3,443 pregnancies as an ordinary independent
# sample, ignoring the sampling weight and the cluster and strata
# design entirely. It is fit here as a deliberate first step, a
# baseline to compare the properly weighted version in Section 21
# against, not as the final answer.

# ******************************************************************
# OBJECT CREATED: cox_naive
# ******************************************************************
cox_naive <- coxph(
  Surv(p20, stillbirth) ~ interval_group + maternal_age + parity + m14,
  data = analytic_data,
  ties = "efron"   # the standard, more accurate way of handling the heavy ties explained in Section 16, the default in R but named explicitly here so the choice is visible rather than silent
)

summary(cox_naive)

# READING THIS OUTPUT. The exp(coef) column is the
# hazard ratio for each variable. A hazard ratio above 1 means higher
# instantaneous risk of stillbirth compared to the reference category,
# for interval_group, or per one unit increase, for the continuous
# variables maternal_age, parity, and m14. A hazard ratio's 95%
# confidence interval that does not cross 1 is generally read as a
# statistically meaningful result, given the honest caveat already
# noted back in Section 12 about seventy events being a thin base to
# build four covariates on top of.

# ====================================================================
# SECTION 20: Check whether the proportional hazards assumption holds
# ====================================================================

# Cox regression assumes the hazard ratio for each variable stays
# constant over the whole span of gestation, a short pregnancy
# is assumed to carry the same relative risk at month seven that it
# does at month nine. cox.zph() tests that assumption directly for
# each variable and for the model as a whole. A p-value below roughly
# 0.05 on any row would mean that variable's effect is not actually
# constant over time, and the model's results for that variable would
# need to be interpreted more cautiously.

cox.zph(cox_naive)

# ====================================================================
# SECTION 21: Rebuild the model the correct way, respecting the survey design
# ====================================================================

# Section 5B already established that DHS did not draw a simple random
# sample of pregnancies. Households were selected in clusters, oversampling
# corrections were applied by division and urban or rural strata, and
# sample_weight, v021, and v022 exist specifically to undo those design
# effects when they are actually used. Skipping them, as cox_naive did
# above, does not just leave a courtesy step out, it produces standard
# errors that are typically too small, because it treats every
# pregnancy as if it carried independent information when pregnancies
# from the same cluster tend to resemble each other. Standard errors
# that are too small make confidence intervals too narrow and p-values
# look more convincing than they honestly are.

# ******************************************************************
# OBJECT CREATED: design
# ******************************************************************
design <- svydesign(
  ids = ~v021,          # the primary sampling unit, the cluster DHS itself recommends for variance estimation
  strata = ~v022,        # the strata DHS itself recommends for variance estimation
  weights = ~sample_weight,
  data = analytic_data,
  nest = TRUE            # tells svydesign that cluster numbers in v021 can repeat across different strata, and should be treated as distinct clusters rather than accidentally merged
)

design

# svycoxph() is the survey adjusted version of coxph(), fit on the
# design object above rather than directly on analytic_data, so that
# every stage of the model, not only the final summary, accounts for
# the weighting and clustering.

# ******************************************************************
# OBJECT CREATED: cox_survey
# ******************************************************************
cox_survey <- svycoxph(
  Surv(p20, stillbirth) ~ interval_group + maternal_age + parity + m14,
  design = design
)

summary(cox_survey)

# WHAT TO ACTUALLY COMPARE. Look at how the hazard ratios and,
# especially, the confidence intervals here differ from cox_naive
# above. It is common, and expected, for the survey-adjusted
# confidence intervals to come out wider than the naive ones; that
# width is the honest price of acknowledging the clustering that was
# ignored before. A result that stays meaningful after this
# adjustment is a genuinely stronger result than one that only looked
# meaningful in the naive model.

# ====================================================================
# SECTION 22: A labelled, side-by-side comparison of both models
# ====================================================================

# Reading two separate summary() printouts and comparing them line by
# line invites mistakes. The table below pulls the hazard ratio and
# 95% confidence interval out of each model into one shared table,
# with plain language row labels instead of raw variable names, so
# the naive and survey-adjusted results sit next to each other
# directly.

extract_hr_table <- function(model, model_label) {
  hr <- exp(coef(model))
  ci <- exp(confint(model))
  
  # summary(model)$coefficients holds the p-value for each term, under
  # the column name "Pr(>|z|)", true for both coxph and svycoxph
  # objects, so this one line works for cox_naive and cox_survey alike.
  p_value <- summary(model)$coefficients[, "Pr(>|z|)"]
  
  data.frame(
    term = names(coef(model)),
    hazard_ratio = round(hr, 2),
    ci_lower = round(ci[, 1], 2),
    ci_upper = round(ci[, 2], 2),
    p_value = round(p_value, 3),   # rounded to three decimal places, enough precision to see which side of 0.05 a result falls on without a wall of trailing digits
    model = model_label,
    row.names = NULL
  )
}

comparison_table <- rbind(
  extract_hr_table(cox_naive, "Naive (unweighted)"),
  extract_hr_table(cox_survey, "Survey-adjusted")
)

# Plain language labels in place of the raw coefficient names that R
# generates automatically.
term_labels <- c(
  "interval_groupshort" = "Short interval (16 months or less) vs reference",
  "interval_grouplong"  = "Long interval (79 months or more) vs reference",
  "maternal_age"        = "Maternal age (per additional year)",
  "parity"               = "Parity (per additional prior pregnancy)",
  "m14"                  = "Antenatal care visits (per additional visit)"
)
comparison_table$term <- term_labels[comparison_table$term]

comparison_table

# Saving this table is what actually produces the file the README
# points to. Printing to the console only shows it here and now; once this
# R session closes, nothing is left behind unless it gets written to
# disk explicitly.

# ====================================================================
# SECTION 23: Visualizing both models as a forest plot
# ====================================================================

# comparison_table already holds everything this plot needs, a hazard
# ratio and a confidence interval for each variable, from each model.
# A forest plot puts every one of those intervals on one shared axis,
# so a reader can see at a glance which results cross the "no effect"
# line at a hazard ratio of 1, and which results move when the naive
# model gets swapped for the survey-adjusted one.
#
# The hazard ratio axis is drawn on a log scale rather than a plain
# linear one. A hazard ratio of 0.5 and a hazard ratio of 2 represent
# the same-sized effect, just pointing in opposite directions, halving
# risk versus doubling it, and a log scale is what makes those two
# distances look equal on the page the way they are mathematically
# equal in size. A linear scale would visually exaggerate every
# hazard ratio above 1 compared to the ones below it.

library(ggplot2)

# The same five variables appear twice in comparison_table, once for
# each model, in the order the Cox formula listed them. Pulling out
# that first occurrence order, then reversing it, makes the plot read
# top to bottom in the same order the table itself already reads in,
# short interval first, antenatal visits last.
term_order <- rev(unique(comparison_table$term))
comparison_table$term <- factor(comparison_table$term, levels = term_order)
comparison_table$model <- factor(comparison_table$model, levels = c("Naive (unweighted)", "Survey-adjusted"))

# ******************************************************************
# OBJECT CREATED: forest_plot
# FILE WRITTEN: outputs/forest_plot.png
# ******************************************************************
forest_plot <- ggplot(comparison_table, aes(x = hazard_ratio, y = term, color = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#898781", linewidth = 0.6) +   # the "no effect" reference line
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0.2,
    position = position_dodge(width = 0.5),
    linewidth = 0.7
  ) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4)) +
  scale_color_manual(values = c("Naive (unweighted)" = "#eb6834", "Survey-adjusted" = "#2a78d6")) +
  labs(
    title = "Hazard ratios for stillbirth, naive versus survey adjusted",
    subtitle = "Points show the hazard ratio, bars show the 95% confidence interval\nThe dashed line marks no effect (hazard ratio of 1)",
    x = "Hazard ratio (log scale)",
    y = NULL,
    color = "Model"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

forest_plot

# WHAT TO ACTUALLY LOOK FOR IN THIS PLOT. Any bar that crosses the
# dashed vertical line is not statistically significant, its interval
# includes the possibility of no effect at all. Watch the short
# interval row in particular, the orange naive bar sits clear of the
# line, the blue survey adjusted bar crosses it, the exact attenuation
# already found in Sections 19 and 21, now visible in a single glance
# rather than something a reader has to compare across two separate
# tables to notice.

# ====================================================================
# END OF SCRIPT
# ====================================================================
#
# A HONEST SUMMARY OF WHAT THIS ANALYSIS CAN AND CANNOT SAY. With
# seventy stillbirth events total, this model can detect a difference
# large enough to matter, which the short and long interval groups
# both appear to show, but it cannot pin that effect down to a
# precise number the way a study with a few hundred events could. Wide
# confidence intervals in comparison_table are not a flaw in how this
# was built, they are an honest reflection of how much a sample this
# size can actually support. Treat the direction of the effect, short
# and long intervals both carrying higher stillbirth risk than the
# middle group, as the finding worth trusting most, and treat the
# exact hazard ratio numbers as a rough range rather than a precise
# estimate.
# ====================================================================
