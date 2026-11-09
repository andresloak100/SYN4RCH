# D-101 power

Produced by `python3 analysis/d101_programmes.py power` before the test was run: real gene
programmes, populations and injection sites with synthetic wiring only.

```
populations 73; eligible within-region pairs 86
programme guidance_adhesion (741 genes); a pass needs p < 0.0167 (the first Holm step of three)
  iid              mean adjusted r -0.008   decision 1 0.015   transfer 0.007
  gene      0.00   mean adjusted r -0.004   decision 1 0.013   transfer 0.007
  gene      0.50   mean adjusted r +0.152   decision 1 0.163   transfer 0.072
  gene      0.75   mean adjusted r +0.272   decision 1 0.438   transfer 0.282
  gene      1.00   mean adjusted r +0.384   decision 1 0.748   transfer 0.583
  position  0.60   mean adjusted r +0.019   decision 1 0.025   transfer 0.018
  position  1.00   mean adjusted r +0.067   decision 1 0.045   transfer 0.033
  cells     1.00   mean adjusted r -0.011   decision 1 0.015   transfer 0.013
programme transcription_factors (1326 genes); a pass needs p < 0.0167 (the first Holm step of three)
  iid              mean adjusted r -0.003   decision 1 0.018   transfer 0.007
  gene      0.00   mean adjusted r -0.014   decision 1 0.020   transfer 0.005
  gene      0.50   mean adjusted r +0.165   decision 1 0.170   transfer 0.098
  gene      0.75   mean adjusted r +0.300   decision 1 0.445   transfer 0.263
  gene      1.00   mean adjusted r +0.447   decision 1 0.807   transfer 0.608
  position  0.60   mean adjusted r +0.026   decision 1 0.018   transfer 0.010
  position  1.00   mean adjusted r +0.064   decision 1 0.033   transfer 0.022
  cells     1.00   mean adjusted r -0.051   decision 1 0.003   transfer 0.000
python3 analysis/d101_programmes.py power  54.54s user 9.70s system 96% cpu 1:06.31 total
```
