# CTPerfusionPipeline

This code replicates the CT perfusion preprocessing pipeline from Teghipco, A., Kim, H., Rorden, C., Newman-Norlund, R., Sharif M., Sikorski D., Hillis, A.E. (in preperation). Excellence is a habit: Enhancing predictions of language impairment by identifying stable features in clinical perfusion scans. It assumes you are working with perfusion measure images output by RAPID and CTA head/neck images. See documentation for full pipeline; the two sets of images will be co-registered, then normalized. Using the resulting transforms, an atlas (e.g., JHU atlas) will be brought into native space. You can then use scalar_atlas.m to extract perfusion measures within each ROI of the atlas. 

Please see documentation in ct_rgb.m for the script that does the majority of the work. Additionally, see runPipeline.m for more information on setting up the pipeline. Note, you will need SPM12 or later to use.
