#!/bin/bash
source /science/InitScripts/iaat-xmmsas.sh # Source the xmmsas environment here
input_image=../Data/Images/merged_image_200_2300.fits
cheesemask=../Data/Source_cat/cheesemask_PS_1.0arcmin.fits
masked_image=../Data/Images/merged_image_200_2300_masked.fits
input_expmap=../Data/Images/merged_expmap_200_2300.fits
masked_expmap=../Data/Images/merged_expmap_200_2300_masked.fits
output_smooth_image=../Data/Images/merged_image_200_2300_asmooth.fits
desiredsnr=30

farith $input_image $cheesemask $masked_image MUL clobber=yes
farith $input_expmap $cheesemask $masked_expmap MUL clobber=yes

asmooth inset=$masked_image     outset=$output_smooth_image     weightset=$masked_expmap     withweightset=yes     withexpimageset=yes     expimageset=$masked_expmap     desiredsnr=$desiredsnr 