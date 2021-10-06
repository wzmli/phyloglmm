# phyloglmm

It looks like

```
make fit.$platform.$numsite.$size.$seed.Rout
```

is supposed to generate a fit using:

- `$platform`: pez phyr lme4 brms gls glmmTMB phylolm ?
- `$numsite`: `ss` (single), `ms` (multiple)
- `$size`: `small`, `medium`, `large`
- `$seed`: integer value

Fit uncorrelated intercept slope glmms using pez and lme4

## Todo

- single site simulation (problematic, not sure why it doesn't work)
- multivariate simulation to get correlation between slopes and intercepts
- extend example with MCMCglmm, glmmTMB, gls, PGLMM, brms
- what is caper? there is a function call pgls
