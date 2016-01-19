# generate final pics from PNGs in the PUBLICATIONS folder ...
# pdflatex Fig1 # Fig1 is a scheme, nothing to overlay here ...
pdflatex Fig2
pdflatex Fig3
pdflatex Fig4
pdflatex Fig5
pdflatex Fig6
pdflatex Fig7
rm *.aux *.log

# let's make PNG for our manuscripts and compressed TIF for submission ...
convert -density 600 ../figure1.png -quality 100 Fig1.png
convert -density 600 ../figure1.png -quality 100 -compress lzw Fig1.tiff
# let's make PNG for our manuscripts and compressed TIF for submission ...
convert -density 600 Fig2.pdf -quality 100 Fig2.png
convert -density 600 Fig2.pdf -quality 100 -compress lzw Fig2.tiff
# let's make PNG for our manuscripts and compressed TIF for submission ...
convert -density 600 Fig3.pdf -quality 100 Fig3.png
convert -density 600 Fig3.pdf -quality 100 -compress lzw Fig3.tiff
# let's make PNG for our manuscripts and compressed TIF for submission ...
convert -density 600 Fig4.pdf -quality 100 Fig4.png
convert -density 600 Fig4.pdf -quality 100 -compress lzw Fig4.tiff
# let's make PNG for our manuscripts and compressed TIF for submission ...
convert -density 600 Fig5.pdf -quality 100 Fig5.png
convert -density 600 Fig5.pdf -quality 100 -compress lzw Fig5.tiff
# let's make PNG for our manuscripts and compressed TIF for submission ...
convert -density 600 Fig6.pdf -quality 100 Fig6.png
convert -density 600 Fig6.pdf -quality 100 -compress lzw Fig6.tiff
# let's make PNG for our manuscripts and compressed TIF for submission ...
convert -density 600 Fig7.pdf -quality 100 Fig7.png
convert -density 600 Fig7.pdf -quality 100 -compress lzw Fig7.tiff


# move PNGs to the article folder ...
mv *.png ../Archive


# that's it ....
# tiffs are only for submission ...




