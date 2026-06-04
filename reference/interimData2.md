# Interim enrollment and event data after enrollment completion

A data frame with 300 rows and 7 columns:

- `trialsdt`:

  The trial start date

- `usubjid`:

  The unique subject ID

- `randdt`:

  The randomization date

- `treatment`:

  The treatment group number

- `treatment_description`:

  Description of the treatment group

- `time`:

  The day of event or censoring since randomization

- `event`:

  The event indicator: 1 for event, 0 for non-event

- `dropout`:

  The dropout indicator: 1 for dropout, 0 for non-dropout

- `cutoffdt`:

  The cutoff date

For ongoing subjects, both `event` and `dropout` are equal to 0.

## Usage

``` r
interimData2
```

## Format

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 300
rows and 9 columns.
