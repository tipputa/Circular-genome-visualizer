# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 14:21:32 2017

@author: tipputa
"""


header="""
# circos.conf

karyotype = data/karyotype.txt
chromosomes_units = 1
#chromosomes_display_default = no
#chromosomes                 = chr1


<ideogram>

<spacing>
default = 0.00r
</spacing>

# Ideogram position, fill and outline

radius           = 0.8r
thickness        = 25p
fill             = yes
#stroke_color     = white #dgrey
#stroke_thickness = 2p

# Minimum definition for ideogram labels.

show_label       = no
# see etc/fonts.conf for list of font names
label_font       = default
label_radius     = 1r + 75p
label_size       = 30
label_parallel   = yes

show_bands = yes
fill_bands = yes
# band_label = yes
# band_transparency = 5

show_ticks          = no
show_tick_labels    = no

# chromosomes_color   = chr1=red
</ideogram>
#<<include data/ticks.conf>>
<highlights> 
"""


tail="""
</highlights>



################################################################
# The remaining content is standard and required. It is imported
# from default files in the Circos distribution.
#
# These should be present in every Circos configuration file and
# overridden as required. To see the content of these files,
# look in etc/ in the Circos distribution.

<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>
"""