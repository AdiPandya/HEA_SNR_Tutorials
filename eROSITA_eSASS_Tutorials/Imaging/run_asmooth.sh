#!/bin/bash
source /science/InitScripts/iaat-xmmsas.sh
input_image=../Data/Images/merged_image_200_2300_masked.fits
cheesemask=../Data/Source_cat/cheesemask_PS_1arcmin.fits
masked_image=../Data/Images/merged_image_200_2300_masked_smoothed.fits
input_expmap=../Data/Images/merged_expmap_200_2300_masked_masked.fits
masked_expmap=../Data/Images/merged_expmap_200_2300_masked_masked.fits
output_smooth_image=../Data/Images/merged_image_200_2300_masked_smoothed.fits
desiredsnr=30

farith $input_image $cheesemask $masked_image MUL clobber=yes
farith $input_expmap $cheesemask $masked_expmap MUL clobber=yes

asmooth inset=$masked_image 
    outset=$output_smooth_image
    weightset=$masked_expmap 
    withweightset=yes
    withexpimageset=yes
    expimageset=$masked_expmap
    desiredsnr=$desiredsnr
