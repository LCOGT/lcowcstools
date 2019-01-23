# lcowcstools
wcstools for lco

This tool is meant to derive WCS SIP (or as appropiate) solutions for LCO BANZAI e91 processed images that already 
have a valid linear WCS solution. 

## Targeted workflow per image

1. Read e91 file: load CD matrix, CRVAL, CRPIX, load source catalog extension
2. Fetch reference catalog (gaia) from area around image
3. Sort reference and source catalog by magnitude

Iterate n times:
 4. transform reference catalog to pixel coordiantes (inverse TAN / CD matrix matrix applied)
 5. Transform source catalog to undistorted pixel coordiantes (in itereation 0, SIP is identity transformation)
 6. Match both catalogs by smalles distance, with configurable max allowable distance
 7. Fit WCS SIP solution 
 8. Decide on n iterations / convergance if to iterate from 4.
 
9. If desired, create graphical output of fit: x/y errors vs x/y coordiante / radius. Bar Diagram of SIP coefficiencts. 
10. Save resulting WCS into sqlite database, together with meta information: seeing, alt, az, exptime, telfocus


## Targeted workflow for evaluation:
1. Read WCS from sqlite database, by camera identifier
2. Plot SIP coefficents vs time, alt, az, seeing, telfocus.
3. Averrage SIP coefficients on TBD criterion from input images and generate average SIP for inclusion into site software processing. 



## Notes:
Parts of the star matching are already implemented in the photzp tools. 
 
 

