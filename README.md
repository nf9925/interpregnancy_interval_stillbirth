Interpregnancy_interval_stillbirth

A survival analysis project built on the 2022 Bangladesh Demographic and Health Survey, asking whether the length of time since a woman's previous pregnancy, alongside her age, how many pregnancies she has had before, and how much antenatal care she received, is associated with her risk of stillbirth in the pregnancy that follows.

Why this question

The United Nations Inter-Agency Group for Child Mortality Estimation (UN IGME) reports that Bangladesh records more than 63,000 stillbirths every year. Roughly one baby is stillborn for every 41 births, the highest rate anywhere in South Asia. Stillbirth also remains one of the more understudied outcomes in global maternal and child health relative to its actual scale, which is part of why this project treats it as the outcome worth building a full analysis around, rather than a side note to live birth outcomes.

The specific research design here follows Stephansson, Dickman, and Cnattingius, "The influence of interpregnancy interval on the subsequent risk of stillbirth and early neonatal death," published in Obstetrics and Gynaecology in 2003 (PMID 12850614). That paper used Swedish national birth registry data. The project asks the same underlying question in a very different setting: a nationally representative household survey in Bangladesh, using DHS's own sampling weights, clusters, and strata to keep that national claim honest.

Data source

This project uses three files from the 2022 Bangladesh Demographic and Health Survey, the Pregnancies Recode (BDGR81FL), the Individual Recode (BDIR81FL), and the Pregnancy and Postnatal Care Recode (BDNR81FL), all distributed by the DHS Program.

None of that data lives in this repository. DHS microdata is only available after registering directly with the DHS Program at dhsprogram.com, and its terms of use do not permit redistributing the files, so anyone wanting to run this script needs to register and download their own copies of the three files listed above, then update the file paths near the top of scripts/analysis.R. When citing DHS data in any output, use the exact citation DHS provides at the point of download, since that wording is set by DHS itself, not reproduced here.

Methods, in brief

The first half of scripts/analysis.R builds one analytic dataset from the three raw files. Stillbirth is defined as a pregnancy ending in fetal death at seven months of gestation or more, matching the same definition on both the stillbirth and live birth sides of the comparison, the same convention Stephansson's paper used. Interpregnancy interval is calculated from the gap between consecutive pregnancy end dates for the same woman, with a dedicated fix for twin pregnancies that would otherwise be double-counted. Because only 70 stillbirth events exist in the final sample, the interval gets collapsed into three categories, short, reference, and long, rather than the finer five group split a larger sample could support.

The second half of the same script builds Kaplan-Meier curves and a log-rank test as a first, unadjusted look at whether the three interval groups differ, then fits a Cox proportional hazards model two ways, once treating the sample as an ordinary independent sample, and once properly accounting for DHS's actual cluster and stratified sampling design through the survey package. Both versions are compared directly, in a labelled table and a forest plot, rather than only reporting whichever version looked more convincing.

Key finding

The naive, unweighted model suggested a real effect, pregnancies following a short interval showed a statistically significant 91 percent higher hazard of stillbirth. That result did not survive proper adjustment for DHS's survey design, the same comparison in the correctly weighted model came back with a hazard ratio of 1.45 and a p-value of 0.23, no longer significant. The one result that held up under both versions of the model was antenatal care, more visits were associated with meaningfully lower stillbirth risk in both the naive and the survey-adjusted model.

The full numbers behind that comparison are in outputs/comparison_table.csv, and the same comparison is drawn out visually in outputs/forest_plot.png.

Repository structure
.
├── README.md
├── LICENSE
├── .gitignore
├── scripts/
│   └── analysis.R                    Builds the analytic dataset, then runs the full survival analysis
└── outputs/
    ├── comparison_table.csv             Hazard ratios, 95% CIs, and p-values, naive vs survey-adjusted
    ├── kaplan_meier_group_summary.csv   Sample size and stillbirth rate by interval group
    └── forest_plot.png                  Visual comparison of both Cox models

Reproducing this analysis

None of the raw DHS files is included in this repository; see Data source above. Register with the DHS Program and download the three files named above. Update the file paths near the top of scripts/analysis.R to point to your own copies, then run the script top to bottom in a single R session; the data prep half and the analysis half both live in this one file and run in order automatically.

Packages needed: haven, dplyr, ggplot2, survival, and survey.

Honest limitations

Seventy stillbirth events is a thin base for a model carrying four covariates, commonly cited guidance suggests something closer to ten to twenty events per predictor. Confidence intervals throughout this project are correspondingly wide, and the exact hazard ratio values are best read as a rough range rather than a precise estimate. Gestational duration is only recorded in completed months in this data, not weeks or days, which produces heavy ties in the survival model, handled here with the Efron approximation rather than ignored. And DHS's antenatal visit count only has non-missing values for a recent recall subset of pregnancies, which is the main reason the final analytic sample is far smaller than the raw number of stillbirths in the full dataset, a property of how DHS collected this data, not a choice made in this analysis.
