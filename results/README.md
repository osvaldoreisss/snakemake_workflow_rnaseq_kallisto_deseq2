This folder is meant to hold all files generated during analyses conducted by the workflow.
Best practice is to have a separate folder for intermediate files of each analysis step and to use separate subfolders for:

* `logs/` of analysis steps, optimally with a subfolder per rule (especially when rules are run a lot times in parallel)
* `plots/` generated by workflow rules
* `tables/` generated by workflow rules

Further, this folder could contain an `tmp/` subfolder if you fear that the system's `tmp` might overflow from your analyses.